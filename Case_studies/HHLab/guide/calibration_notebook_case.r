rm(list = ls())
graphics.off()

# Set directory
dir_workspace <- here::here()

# Load functions
# function_list <- list.files("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/Functions/", full.names = TRUE)
# for (i in function_list) {
#     source(i)
# }



#############################################
# Load Functions
#############################################

# Interpolation function passing by some specific nodes:
interpolation_specific_points <- function(total_points = 100,
                                          specific_nodes) {
    if (total_points <= length(specific_nodes)) stop("Total points defined is lower or equal to the specific points required. Please either increase the number of total_points or not used interpolation function")
    # Calculate the intervals between the specific points
    intervals <- diff(specific_nodes)
    # Distribute the number of points proportionally to each interval
    total_local_points <- total_points - 1 # last value during diff function
    points_per_interval_float <- (intervals / sum(intervals)) * total_local_points

    points_per_interval <- round(points_per_interval_float)

    # Generate sub-sequences for each interval, excluding the duplicate endpoint
    specific_seq_temp <- unlist(mapply(
        function(start, end, n) seq(start, end, length.out = n + 1)[-(n + 1)], # Exclude the endpoint
        specific_nodes[-length(specific_nodes)], # Start of each interval
        specific_nodes[-1], # End of each interval
        points_per_interval # Number of points per interval
    ))

    # Add the final endpoint manually
    specific_seq <- sort(
        c(
            specific_seq_temp,
            specific_nodes[
                -c(1)
            ]
        )
    )

    if (!all(specific_nodes %in% specific_seq)) stop("Distance so close, increase number of interpolate KP to get all points required")

    return(specific_seq)
}

# Function to assign reach to a grid
assign_reach_from_a_grid <- function(reach_KP_boundaries, grid) {
    reaches <- rep(NA, length(grid))

    if (nrow(reach_KP_boundaries) == 1) {
        all_set <- dplyr::between(grid, reach_KP_boundaries$KP_start, reach_KP_boundaries$KP_end)
        if (!all(all_set)) stop("reach_KP_boundaries must contain all grid")
        reaches[which(all_set)] <- reach_KP_boundaries$reach[1]
    } else {
        duplicated_numbers <- grid[duplicated(grid)]
        if ((nrow(reach_KP_boundaries) - 1) != length(duplicated_numbers)) stop("Number of rows of specific_points_Model minus 1 must coincide with the number of repeated values in the grid")

        for (i in 1:nrow(reach_KP_boundaries)) {
            idx_position_reaches <- which(dplyr::between(grid, reach_KP_boundaries$KP_start[i], reach_KP_boundaries$KP_end[i]))

            idx_duplicated_values <- which(grid[idx_position_reaches] %in% duplicated_numbers)

            if (length(idx_duplicated_values) != 0) {
                if (length(idx_duplicated_values) != 2 & length(idx_duplicated_values) != 4) stop("Duplicated values are not found, the number of duplicated values must be two or four")

                # First duplicated values
                if (i == 1 & length(idx_duplicated_values) == 2) {
                    reaches[idx_position_reaches] <- reach_KP_boundaries$reach[i]

                    reaches[idx_position_reaches[idx_duplicated_values[2]]] <- reach_KP_boundaries$reach[i + 1]
                } else if (length(idx_duplicated_values) == 4) {
                    idx_position_reaches <- idx_position_reaches[-idx_duplicated_values[c(1, 2)]]

                    reaches[idx_position_reaches] <- reach_KP_boundaries$reach[i]

                    reaches[idx_position_reaches[(idx_duplicated_values[4] - 2)]] <- reach_KP_boundaries$reach[i + 1] # -2 because two values have been removed previously
                } else if (length(idx_duplicated_values) == 2) {
                    idx_position_reaches <- idx_position_reaches[-idx_duplicated_values[c(1, 2)]]
                    reaches[idx_position_reaches] <- reach_KP_boundaries$reach[i]
                }
            } else {
                stop("It should never get here. Something is wrong")
                reaches[idx_position_reaches] <- reach_KP_boundaries$reach[i]
            }
        }
    }

    if (any(is.na(reaches))) {
        return(stop("Some values in the vector reaches have not been assigned"))
    }

    return(reaches)
}

# Function to get covariate in piecewise function
getCovariate_piecewise <- function(KP_grid, shiftPoints) {
    if (any(diff(KP_grid) < 0)) {
        return(stop("KP grid must be increasing all the time"))
    }
    if (KP_grid[1] == min(shiftPoints)) {
        return(stop("Shift point must be different to the first KP_grid"))
    }
    if (last(KP_grid) == max(shiftPoints)) {
        return(stop("Shift point must be different to the last KP_grid"))
    }
    if (!any(KP_grid < min(shiftPoints))) {
        return(stop(paste0("Shift point ", min(shiftPoints), "is out of the KP grid")))
    }
    if (!any(KP_grid > max(shiftPoints))) {
        return(stop(paste0("Shift point ", max(shiftPoints), " is out of the KP grid")))
    }
    if (any(duplicated(shiftPoints))) {
        return(stop("Non duplicated values in shiftPoints"))
    }

    all_shiftPoints <- sort(c(
        shiftPoints,
        KP_grid[1],
        last(KP_grid)
    ))

    # Assign each value in KP_grid to an interval
    intervals <- cut(KP_grid, breaks = all_shiftPoints, include.lowest = TRUE, right = TRUE)

    # Create binary matrix version 0:
    covariate_matrix_v0 <- stats::model.matrix(~ intervals - 1)
    covariate_matrix_v1 <- cbind(KP_grid, covariate_matrix_v0)
    # Get duplicated values from grid respresenting switch of reach
    id_shifts_reaches <- which(duplicated(KP_grid) | duplicated(KP_grid, fromLast = TRUE))

    if (length(id_shifts_reaches) != 0) {
        # Identify which of these are the second (or later) occurrences
        is_second_occurrence <- duplicated(KP_grid[id_shifts_reaches])
        # If TRUE, matrix should be corrected to consider shift point coincide with switch of reach
        # Correct the matrix for the second (or later) occurrences
        for (i in which(is_second_occurrence)) {
            id <- id_shifts_reaches[i]
            # Find the column where the current 1 is set
            current_col <- which(covariate_matrix_v1[id, ] == 1)
            # If the current column is not the last one, set the next column to 1 and the current to 0
            if (current_col < ncol(covariate_matrix_v1)) {
                covariate_matrix_v1[id, current_col] <- 0
                covariate_matrix_v1[id, current_col + 1] <- 1
            }
        }
    }
    covariate_matrix <- covariate_matrix_v1[, -c(1)]
    colnames(covariate_matrix) <- NULL
    return(covariate_matrix)
}

# Function to estimate covariate's values of a polynomial degree fixed
getlegendre <- function(
    polynomial_degree,
    covariate_discretization) {
    if (length(unique(covariate_discretization)) == 1) {
        stop("covariate_discretization must contain more than one unique value")
    }
    normalized_values <- 2 * (covariate_discretization - min(covariate_discretization)) / (max(covariate_discretization) - min(covariate_discretization)) - 1

    # Write co variant matrix for main channel
    # Create an empty data frame and fill it with polynomial values
    legendre_df_covariate <- data.frame(x = normalized_values) # Start with x values
    legendre_df <- legendre_df_covariate

    if (polynomial_degree < 0) stop("The polynomial degree must be positive")
    if (any(!dplyr::between(normalized_values, -1, 1))) stop("Range of normalized_values should be between [-1,1]")

    if (polynomial_degree == 0) {
        return(rep(1, length(normalized_values)))
    }
    if (polynomial_degree == 1) {
        return(normalized_values)
    }

    # P_0(normalized_values) = 1
    Pn_1 <- rep(1, length(normalized_values))
    # P_1(normalized_values) = normalized_values
    Pn <- normalized_values

    for (k in 1:(polynomial_degree - 1)) {
        P_next <- ((2 * k + 1) * normalized_values * Pn - k * Pn_1) / (k + 1)
        Pn_1 <- Pn
        Pn <- P_next
    }
    return(Pn)
}

# Function to get covariate in Legendre polynomial
getCovariate_Legendre <- function(max_polynomial_degree, ...) {
    if (!(max_polynomial_degree >= 0 && max_polynomial_degree == floor(max_polynomial_degree))) stop("max_polynomial_degree must be a positive integer")

    # Capture the first argument passed via ...
    args <- list(...)
    # Check if required arguments are present
    if (!("covariate_discretization" %in% names(args))) {
        stop("covariate_discretization must be provided as a named argument")
    }
    # Extract covariate_discretization
    covariate_discretization <- args$covariate_discretization

    Z <- matrix(
        data = 0,
        nrow = length(covariate_discretization),
        ncol = max_polynomial_degree + 1 # Because degrees start at 0
    )
    for (polynomial_degree in seq(0, max_polynomial_degree)) {
        Z[, polynomial_degree + 1] <- getlegendre(
            polynomial_degree = polynomial_degree,
            covariate_discretization = covariate_discretization
        )
    }

    return(Z)
}

# Block-diagonal binding of matrices
# source: https://stackoverflow.com/questions/17495841/block-diagonal-binding-of-matrices
block_diagonal_matrix <- function(...) {
    d <- list(...)
    nrows <- sum(sapply(d, NROW))
    ncols <- sum(sapply(d, NCOL))
    ans <- matrix(0, nrows, ncols)
    i1 <- 1
    j1 <- 1
    for (m in d) {
        if (!is.matrix(m)) {
            return(stop("All arguments must be matrices"))
        }
        i2 <- i1 + NROW(m) - 1
        j2 <- j1 + NCOL(m) - 1
        ans[i1:i2, j1:j2] <- m
        i1 <- i2 + 1
        j1 <- j2 + 1
    }
    return(ans)
}

RUGFile_constructor_SU_correction <- function(data_layer_1) {
    # Calculate the length of RUGFile using data_layer_1
    length_RUGFile <- sum(
        sapply(data_layer_1, function(id_data_layer_1) {
            length(id_data_layer_1$reach) - length(unique(id_data_layer_1$reach))
        })
    )
    # Initialize RUGFile_first_part as a data frame with 5 columns
    RUGFile_first_part <- data.frame(
        matrix(
            data = 0,
            nrow = length_RUGFile,
            ncol = 5 # RUGFILE columns : reach, kp_start, kp_end, kmin,kmoy
        )
    )
    colnames(RUGFile_first_part) <- c("reach", "KP_start", "KP_end", "Kmin", "Kmoy")
    # Correct information from data_layer_1 to adapt it to RUGFile, because size of both must be equal!
    data_layer_1_adapted <- vector(mode = "list", length = length(data_layer_1))
    i1 <- 1
    counter_SU <- 1
    for (id_data_layer_1 in data_layer_1) {
        id_reaches <- id_data_layer_1$reach
        id_KP <- id_data_layer_1$KP

        for (i in unique(id_reaches)) {
            mask <- which(id_reaches == i)
            i2 <- i1 + length(mask) - 2 # start in i1, end in length of mask (all data of a specific reach) - 2 (first data already count in mask and remove last value because the count is about intervals)
            RUGFile_first_part[i1:i2, c(1, 2, 3)] <- data.frame(
                id_reach = id_reaches[mask][-c(1)],
                KP_start = id_KP[mask][-c(length(id_KP[mask]))],
                KP_end = id_KP[mask][-c(1)]
            )
            i1 <- i2 + 1
        }

        # Correct SU removing second data of interpolated grid in KP, Reach and Z by id_reach
        SU_before_adaptation <- data.frame(
            KP = id_data_layer_1$KP,
            reach = id_data_layer_1$reach,
            Z = id_data_layer_1$Z
        )

        if (length(unique(id_reaches)) >= 1) {
            SU_adapted_sliced <- data.frame(SU_before_adaptation %>%
                group_by(reach) %>%
                slice(-2) %>% # Remove the second row of each group
                ungroup())

            data_layer_1_adapted[[counter_SU]]$KP <- SU_adapted_sliced$KP
            data_layer_1_adapted[[counter_SU]]$reach <- SU_adapted_sliced$reach
            data_layer_1_adapted[[counter_SU]]$Z <- as.matrix(SU_adapted_sliced[, -c(1, 2)])
            data_layer_1_adapted[[counter_SU]]$prior <- id_data_layer_1$prior
        } else {
            return(stop("id_reaches must be higher than 0"))
        }
        counter_SU <- counter_SU + 1
    }
    return(list(
        RUGFile = RUGFile_first_part,
        data_layer_1_adapted = data_layer_1_adapted
    ))
}

CalData_plot <- function(
    data,
    scales_free = "free_y",
    y_label,
    title_label,
    col_label = NULL,
    plot_water_depth = TRUE,
    wrap = TRUE) {
    if (plot_water_depth) {
        plot_CalData <-
            ggplot(data, aes(
                x = x,
                color = ID,
                y = h_mean
            ))

        if (any(data$Yu != 0)) {
            plot_CalData <- plot_CalData +
                geom_errorbar(aes(
                    ymin = h_mean - 1.96 * Yu,
                    ymax = h_mean + 1.96 * Yu
                ))
        }
    } else {
        plot_CalData <-
            ggplot(data = data, aes(
                x = x,
                y = z_mean,
                col = ID
            )) +
            geom_line(
                aes(
                    y = z_riverbed,
                    col = "riverbed"
                )
            )

        if (any(data$Yu != 0)) {
            plot_CalData <- plot_CalData +
                geom_errorbar(aes(
                    ymin = z_mean - 1.96 * Yu,
                    ymax = z_mean + 1.96 * Yu
                ))
        }
    }
    plot_CalData <- plot_CalData +
        geom_point() +
        labs(
            x = "Streamwise position (meters)",
            y = y_label,
            title = title_label,
            col = col_label
        ) +
        theme_bw() +
        theme(
            plot.title = element_text(hjust = 0.5),
            legend.title = element_text(hjust = 0.5)
        )

    if (wrap) {
        plot_CalData <- plot_CalData +
            facet_wrap(~ID, scales = scales_free, ncol = 1)
    }
    return(plot_CalData)
}
get_specific_SR_points <- function(SR_properties) {
    if (identical(SR_properties$function_SR, getCovariate_piecewise)) {
        SR_shifts_points <- SR_properties$shiftPoints
        return(SR_shifts_points)
    }
}
`%out%` <- function(x, y) {
    x[!x %in% y]
}


Get_KP_Reach_SR <- function(SR_key_HM,
                            data_layer_1_ID_MC_TRIB) {
    # Initialize KP and reaches as lists of lists
    KP <- reaches <- vector(mode = "list", length = length(SR_key_HM))
    names(KP) <- names(reaches) <- names(SR_key_HM)
    # Loop through each main channel and tributary
    for (ID_MC_TRIB in seq_along(SR_key_HM)) {
        # Get KP boundaries of each SR
        SR_KP_boundaries_list <- SR_key_HM[[ID_MC_TRIB]]
        data_KP <- data_layer_1_ID_MC_TRIB[[ID_MC_TRIB]]$KP
        data_reaches <- data_layer_1_ID_MC_TRIB[[ID_MC_TRIB]]$reach
        # Check if KP boundaries of each SR respect the KP boundaries of ID
        if (min(data_KP) != min(unlist(SR_KP_boundaries_list)) ||
            max(data_KP) != max(unlist(SR_KP_boundaries_list))) {
            return(stop(sprintf(
                "KP boundary mismatch: data min/max (%s, %s) vs HM min/max (%s, %s)",
                min(data_KP), max(data_KP),
                min(unlist(SR_KP_boundaries_list)), max(unlist(SR_KP_boundaries_list))
            )))
        }

        # Initialize sublists for KP and reaches
        KP[[ID_MC_TRIB]] <- reaches[[ID_MC_TRIB]] <- vector("list", length(SR_KP_boundaries_list))
        names(KP[[ID_MC_TRIB]]) <- names(reaches[[ID_MC_TRIB]]) <- names(SR_KP_boundaries_list)
        # Loop through each SR
        for (i in seq_along(SR_KP_boundaries_list)) {
            boundaries <- SR_KP_boundaries_list[[i]]

            # Check if boundaries are within data_KP range
            if (any(!dplyr::between(boundaries, min(data_KP), max(data_KP)))) {
                return(stop("Values in SR_key_HM are outside the range of Input_ID_key_main_channel_tributaries"))
            }

            # SR traitement : Extract KP and reaches within boundaries
            position_match <- which(data_KP >= boundaries[1] & data_KP <= boundaries[2])
            SR_KP_temp <- data_KP[position_match]
            SR_reach_temp <- data_reaches[position_match]

            # Handle duplicated KP values
            # Check if duplicated values exist indicating switch in MAGE reach
            if (any(duplicated(SR_KP_temp))) {
                duplicated_values <- SR_KP_temp[duplicated(SR_KP_temp)]
                # Check  if boundaries conditions match with duplicated data (switch of MAGE reach) to treat separately
                if (any(duplicated_values %in% boundaries)) {
                    # Correct position_match if boundaries are duplicated
                    id_boundaries_duplicated_values <- which(duplicated_values == boundaries)

                    # First check: detect if a intermedial node, duplicated values are equal to both boundaries
                    # Solution, remove first and last value
                    if (length(id_boundaries_duplicated_values) == 2) {
                        position_match <- position_match[-c(1, length(position_match))]
                    } else if (id_boundaries_duplicated_values == 1) {
                        # Second check: detect if duplicated value is first boundaries (left node of the reach, commonly case of final MAGE reach)
                        # Solution, remove first value
                        position_match <- position_match[-1]
                    } else {
                        # Otherwise: duplicated value is the second element of boundaries (right node of the reach, commonly case of initial MAGE reach)
                        # Solution, remove last value
                        position_match <- position_match[-length(position_match)]
                    }
                    KP[[ID_MC_TRIB]][[i]] <- data_KP[position_match]
                    reaches[[ID_MC_TRIB]][[i]] <- data_reaches[position_match]
                } else { # Duplicated nodes are not in boundary conditions
                    KP[[ID_MC_TRIB]][[i]] <- SR_KP_temp
                    reaches[[ID_MC_TRIB]][[i]] <- SR_reach_temp
                }
            } else { # No duplicated values
                KP[[ID_MC_TRIB]][[i]] <- SR_KP_temp
                reaches[[ID_MC_TRIB]][[i]] <- SR_reach_temp
            }
        }
    }
    return(list(KP = KP, reaches = reaches))
}

############################################
# End source functions
############################################


############################################
# Calibration settings
############################################

############################################
# Construction of SU, RUGFile and Z_file (spatialisation)
############################################

# Input information about hydraulic model (HM)
Input_reachesHM <- data.frame(
    reach = c(1, 2, 3, 4),
    KP_start = c(0, 5, 10, 25),
    KP_end = c(10, 15, 20, 35)
)

############################################
# SU in Kmin and Kflood
############################################
# First layer: ID aggregated reaches: main channel or tributary
Input_ID_key_main_channel_tributaries <- list(
    main_channel = c(1, 3),
    tributary_1 = 2,
    tributary_2 = 4
)


# Second layer: Kmin
# Key to relate KP_cut by spatialisation reach (SR) to KP discretized
Input_Kmin_SR_key_reachesHM <- list(
    main_channel = list(
        SU1 = c(0, 20) # only a SU going from KP = 0 to 20, representing start KP of reach 1 and last KP of reach 3
    ),
    tributary_1 = list(
        SU1 = c(5, 15)
    ),
    tributary_2 = list(
        SU1 = c(25, 35)
    )
)

# Check if size is respected between ID and Kmin
if (length(Input_ID_key_main_channel_tributaries) != length(Input_Kmin_SR_key_reachesHM)) stop("Size must be equal between ID key list and Kmin SR")

# Second layer: Main channel
# Third layer: SR_1
Kmin_Z_properties_MC_SR_1 <- list(
    function_SR = getCovariate_piecewise,
    shiftPoints = c(5, 10, 14, 18) # KP to defined the shift point of constant friction
)
# Second layer: Tributary 1
# Third layer: SR_1
Kmin_Z_properties_TB1_SR_1 <- list(
    function_SR = getCovariate_Legendre,
    max_polynomial_degree = 2
)
# Second layer: Tributary 2
# Third layer: SR_1
Kmin_Z_properties_TB2_SR_1 <- list(
    function_SR = getCovariate_Legendre,
    max_polynomial_degree = 1
)
# Second layer : Kflood
# Key to relate KP_cut by spatialisation reach (SR) to KP discretized
Input_Kflood_SR_key_reachesHM <- list(
    main_channel = list(
        SU1 = c(0, 10),
        SU2 = c(10, 20)
    ),
    tributary_1 = list(
        SU1 = c(5, 15)
    ),
    tributary_2 = list(
        SU1 = c(25, 35)
    )
)
# Check if size is respected between ID and Kmin
if (length(Input_ID_key_main_channel_tributaries) != length(Input_Kflood_SR_key_reachesHM)) stop("Size must be equal between ID key list and Kflood SR")
# Second layer: Main channel
# Third layer: SR_1
Kflood_Z_properties_MC_SR_1 <- list(
    function_SR = getCovariate_piecewise,
    shiftPoints = c(2.5, 7)
)
# Second layer: Main channel
# Third layer: SR_2
Kflood_Z_properties_MC_SR_2 <- list(
    function_SR = getCovariate_Legendre,
    max_polynomial_degree = 2
)
# Second layer: Tributary 1
# Third layer: SR_1
Kflood_Z_properties_TB1_SR_1 <- list(
    function_SR = getCovariate_piecewise,
    shiftPoints = c(8)
)
# Second layer: Tributary 2
# Third layer: SR_1
Kflood_Z_properties_TB2_SR_1 <- list(
    function_SR = getCovariate_Legendre,
    max_polynomial_degree = 0
)

# Create the link between Z_properties_SR either in Kmin or Kmoy to Input_ID_key_main_channel_tributaries
Z_properties_ID_1 <- list(
    Kmin_Z_properties_MC_SR_1,
    Kflood_Z_properties_MC_SR_1,
    Kflood_Z_properties_MC_SR_2
)
Z_properties_ID_2 <- list(
    Kmin_Z_properties_TB1_SR_1,
    Kflood_Z_properties_TB1_SR_1
)
Z_properties_ID_3 <- list(
    Kmin_Z_properties_TB2_SR_1,
    Kflood_Z_properties_TB2_SR_1
)

Z_properties_all <- list(
    Z_properties_ID_1,
    Z_properties_ID_2,
    Z_properties_ID_3
)
if (length(Input_ID_key_main_channel_tributaries) != length(Z_properties_all)) stop(paste0("Number of components in Input_ID_key_main_channel_tributaries (", length(Input_ID_key_main_channel_tributaries), ") must be equal to the number of components in Z_properties_all (", length(Z_properties_all), ")"))

# Get nodes needed for doing the interpolation-> boundaries KP of each MAGE reach from each Input_ID_key_main_channel_tributaries
data_layer_1_ID_MC_TRIB <- vector(mode = "list", length = length(Input_ID_key_main_channel_tributaries))
for (i in seq_along(data_layer_1_ID_MC_TRIB)) {
    # Get hydraulic model information by SR
    data_layer_1_ID_MC_TRIB[[i]]$ID_MC_TRIB_HM_SR <- Input_reachesHM[which(Input_reachesHM$reach %in% Input_ID_key_main_channel_tributaries[[i]]), ]

    # Get nodes of MAGE reaches
    MAGE_nodes <- unique(c(
        data_layer_1_ID_MC_TRIB[[i]]$ID_MC_TRIB_HM_SR$KP_start,
        data_layer_1_ID_MC_TRIB[[i]]$ID_MC_TRIB_HM_SR$KP_end
    ))
    # Add some specific points related to spatialisation Spatialisation reach (SR)
    SR_nodes <- sort(unlist(sapply(Z_properties_all[[i]], function(SR_properties) {
        get_specific_SR_points(SR_properties = SR_properties)
    })))

    data_layer_1_ID_MC_TRIB[[i]]$MAGE_nodes <- MAGE_nodes
    data_layer_1_ID_MC_TRIB[[i]]$SR_nodes <- SR_nodes

    # Interpolated values passing by MAGE_nodes and duplicating them to respect the switch between reaches
    discretization_KP_MAGE <- interpolation_specific_points(
        total_points = 100,
        specific_nodes = MAGE_nodes
    )

    # Remove SR_nodes already presented in MAGE_nodes
    extra_nodes <- SR_nodes %out% MAGE_nodes

    # Find values in extra_nodes that are not in discretization_KP_MAGE
    new_nodes_to_add <- extra_nodes[!extra_nodes %in% discretization_KP_MAGE]

    # KP by ID_MC_TRIB
    data_layer_1_ID_MC_TRIB[[i]]$KP <- sort(c(discretization_KP_MAGE, new_nodes_to_add))
    # 2. ID reach HM
    data_layer_1_ID_MC_TRIB[[i]]$reach <- assign_reach_from_a_grid(
        reach_KP_boundaries = data_layer_1_ID_MC_TRIB[[i]]$ID_MC_TRIB_HM_SR,
        grid = data_layer_1_ID_MC_TRIB[[i]]$KP
    )
}


############################################
# Kmin view
############################################

# First of all, declare where SU will be stored

Kmin_SR <- vector(mode = "list", length = length(Input_ID_key_main_channel_tributaries))
names(Kmin_SR) <- names(Input_ID_key_main_channel_tributaries)


# Get KP and reach of all SR (third layer) in Kmin (second layer)
Kmin_KP_reach_list <- Get_KP_Reach_SR(
    SR_key_HM = Input_Kmin_SR_key_reachesHM,
    data_layer_1_ID_MC_TRIB = data_layer_1_ID_MC_TRIB
)

Kmin_KP <- Kmin_KP_reach_list$KP
Kmin_reach <- Kmin_KP_reach_list$reach

# Get covariate matrix  and prior and the SU by ID (main channel or tributary)
# ID = Main channel
# SR = 1
KP_MC_SR_1 <- Kmin_KP[[1]][[1]]
reach_MC_SR_1 <- Kmin_reach[[1]][[1]]

Z <- Kmin_Z_properties_MC_SR_1[[1]](
    KP_grid = KP_MC_SR_1,
    shiftPoints = Kmin_Z_properties_MC_SR_1[[2]]
)

prior <- list(
    RBaM::parameter(name = "alpha0_1", init = 35),
    RBaM::parameter(name = "alpha0_2", init = 55),
    RBaM::parameter(name = "alpha0_3", init = 30),
    RBaM::parameter(name = "alpha0_4", init = 20),
    RBaM::parameter(name = "alpha0_5", init = 25)
)

if (dim(Z)[2] != length(prior)) stop(paste0("Number of columns of Z (", dim(Z)[2], ") must be equal to the length(prior) = ", length(prior)))


Kmin_SR[[1]]$SU1 <- list(
    KP = KP_MC_SR_1,
    reach = reach_MC_SR_1,
    Z = Z,
    prior = prior
)

# Get covariate matrix  and prior and the SU by ID (main channel or tributary)
# ID = Tributary 1
# SR = 1
KP_TB1_SR_1 <- Kmin_KP[[2]][[1]]
reach_TB1_SR_1 <- Kmin_reach[[2]][[1]]

Z <- Kmin_Z_properties_TB1_SR_1[[1]](
    max_polynomial_degree = Kmin_Z_properties_TB1_SR_1[[2]],
    covariate_discretization = KP_TB1_SR_1
)

prior <- list(
    RBaM::parameter(name = "alpha0", init = 30),
    RBaM::parameter(name = "alpha1", init = 0),
    RBaM::parameter(name = "alpha2", init = 0)
)

if (dim(Z)[2] != length(prior)) stop(paste0("Number of columns of Z (", dim(Z)[2], ") must be equal to the length(prior) = ", length(prior)))

Kmin_SR[[2]]$SU1 <- list(
    KP = KP_TB1_SR_1,
    reach = reach_TB1_SR_1,
    Z = Z,
    prior = prior
)

# Get covariate matrix  and prior and the SU by ID (main channel or tributary)
# ID = Tributary 2
# SR = 1
KP_TB2_SR_1 <- Kmin_KP[[3]][[1]]
reach_TB2_SR_1 <- Kmin_reach[[3]][[1]]

Z <- Kmin_Z_properties_TB2_SR_1[[1]](
    max_polynomial_degree = Kmin_Z_properties_TB2_SR_1[[2]],
    covariate_discretization = KP_TB2_SR_1
)

prior <- list(
    RBaM::parameter(name = "alpha0", init = 30),
    RBaM::parameter(name = "alpha1", init = 0)
)

if (dim(Z)[2] != length(prior)) stop(paste0("Number of columns of Z (", dim(Z)[2], ") must be equal to the length(prior) = ", length(prior)))

Kmin_SR[[3]]$SU1 <- list(
    KP = KP_TB2_SR_1,
    reach = reach_TB2_SR_1,
    Z = Z,
    prior = prior
)
############################################
# End Kmin view
############################################

############################################
# Kflood view
############################################

# First of all, declare where SU will be stored

Kflood_SR <- vector(mode = "list", length = length(Input_ID_key_main_channel_tributaries))
names(Kflood_SR) <- names(Input_ID_key_main_channel_tributaries)


# Get KP and reach of all SR (third layer) in Kflood (second layer)
Kflood_KP_reach_list <- Get_KP_Reach_SR(
    SR_key_HM = Input_Kflood_SR_key_reachesHM,
    data_layer_1_ID_MC_TRIB = data_layer_1_ID_MC_TRIB
)

Kflood_KP <- Kflood_KP_reach_list$KP
Kflood_reach <- Kflood_KP_reach_list$reach

# Get covariate matrix  and prior and the SU by ID (main channel or tributary)
# ID = Main channel
# SR = 1
KP_MC_SR_1 <- Kflood_KP[[1]][[1]]
reach_MC_SR_1 <- Kflood_reach[[1]][[1]]

Z <- Kflood_Z_properties_MC_SR_1[[1]](
    KP_grid = KP_MC_SR_1,
    shiftPoints = Kflood_Z_properties_MC_SR_1[[2]]
)

prior <- list(
    RBaM::parameter(name = "alpha0_1", init = 35),
    RBaM::parameter(name = "alpha0_2", init = 55),
    RBaM::parameter(name = "alpha0_3", init = 30)
)

if (dim(Z)[2] != length(prior)) stop(paste0("Number of columns of Z (", dim(Z)[2], ") must be equal to the length(prior) = ", length(prior)))


Kflood_SR[[1]]$SU1 <- list(
    KP = KP_MC_SR_1,
    reach = reach_MC_SR_1,
    Z = Z,
    prior = prior
)

Kflood_KP <- Kflood_KP_reach_list$KP
Kflood_reach <- Kflood_KP_reach_list$reach

# Get covariate matrix  and prior and the SU by ID (main channel or tributary)
# ID = Main channel
# SR = 2
KP_MC_SR_2 <- Kflood_KP[[1]][[2]]
reach_MC_SR_2 <- Kflood_reach[[1]][[2]]

Z <- Kflood_Z_properties_MC_SR_2[[1]](
    max_polynomial_degree = Kflood_Z_properties_MC_SR_2[[2]],
    covariate_discretization = KP_MC_SR_2
)

prior <- list(
    RBaM::parameter(name = "alpha0", init = 35),
    RBaM::parameter(name = "alpha1", init = 55),
    RBaM::parameter(name = "alpha2", init = 30)
)

if (dim(Z)[2] != length(prior)) stop(paste0("Number of columns of Z (", dim(Z)[2], ") must be equal to the length(prior) = ", length(prior)))


Kflood_SR[[1]]$SU2 <- list(
    KP = KP_MC_SR_2,
    reach = reach_MC_SR_2,
    Z = Z,
    prior = prior
)

# Get covariate matrix  and prior and the SU by ID (main channel or tributary)
# ID = Tributary 1
# SR = 1
KP_TB1_SR_1 <- Kflood_KP[[2]][[1]]
reach_TB1_SR_1 <- Kflood_reach[[2]][[1]]

Z <- Kflood_Z_properties_TB1_SR_1[[1]](
    shiftPoints = Kflood_Z_properties_TB1_SR_1[[2]],
    KP_grid = KP_TB1_SR_1
)

prior <- list(
    RBaM::parameter(name = "alpha0_1", init = 30),
    RBaM::parameter(name = "alpha0_2", init = 20)
)

if (dim(Z)[2] != length(prior)) stop(paste0("Number of columns of Z (", dim(Z)[2], ") must be equal to the length(prior) = ", length(prior)))

Kflood_SR[[2]]$SU1 <- list(
    KP = KP_TB1_SR_1,
    reach = reach_TB1_SR_1,
    Z = Z,
    prior = prior
)

# Get covariate matrix  and prior and the SU by ID (main channel or tributary)
# ID = Tributary 2
# SR = 1
KP_TB2_SR_1 <- Kflood_KP[[3]][[1]]
reach_TB2_SR_1 <- Kflood_reach[[3]][[1]]

Z <- Kflood_Z_properties_TB2_SR_1[[1]](
    max_polynomial_degree = Kflood_Z_properties_TB2_SR_1[[2]],
    covariate_discretization = KP_TB2_SR_1
)

prior <- list(
    RBaM::parameter(name = "alpha0", init = 30)
)

if (dim(Z)[2] != length(prior)) stop(paste0("Number of columns of Z (", dim(Z)[2], ") must be equal to the length(prior) = ", length(prior)))

Kflood_SR[[3]]$SU1 <- list(
    KP = KP_TB2_SR_1,
    reach = reach_TB2_SR_1,
    Z = Z,
    prior = prior
)
########################
# End Flood View
########################


########################
# RUGFile
########################

# Create RUGFile data and adapt SR for correcting grid and respect restriction of sizes between RUGFile and Z_file (spatialisaiton)
RUGFile_and_SU_corrected <- RUGFile_constructor_SU_correction(SR = Kmin_SR)

RUGFile <- RUGFile_and_SU_corrected$RUGFile
Kmin_SR_adapted <- RUGFile_and_SU_corrected$SR_adapted

# Get all Z and reach blocks
Kmin_Z_matrix_block <- block_diagonal_matrix(
    Kmin_SR_adapted[[1]]$Z,
    Kmin_SR_adapted[[2]]$Z,
    Kmin_SR_adapted[[3]]$Z
)

if (nrow(RUGFile) != nrow(Kmin_Z_matrix_block)) stop("RugFile must have the same size as Z file (spatialisation)")

#######################################
# End Kmin and RUGFile
#######################################


############################################
# Kflood settings
############################################

Input_Kflood_SR_key_reachesHM <- list(
    1,
    3,
    2,
    4
)

# SR_1
Kflood_Z_properties_MC_SR_1 <- list(
    function_SR = getCovariate_piecewise,
    shiftPoints = c(2, 7)
)

Kflood_SR <- Kflood_HM_SR <- KP <- reaches <- vector(mode = "list", length = length(Input_Kflood_SR_key_reachesHM))

for (i in seq_along(Input_Kflood_SR_key_reachesHM)) {
    # Get hydraulic model information by SR
    Kflood_HM_SR[[i]] <- Input_reachesHM[which(Input_reachesHM$reach %in% Input_Kflood_SR_key_reachesHM[[i]]), ]
    # SR1
    KP[[i]] <- interpolation_specific_points(
        total_points = 100,
        specific_nodes = nodes_to_interpolate
    )

    # 2. ID reach HM
    reaches[[i]] <- assign_reach_from_a_grid(
        reach_KP_boundaries = Kflood_HM_SR[[i]],
        grid = KP[[i]]
    )
}


Z <- Kflood_Z_properties_MC_SR_1[[1]](
    KP_grid = KP[[1]],
    shiftPoints = Kflood_Z_properties_MC_SR_1[[2]]
)

prior <- list(
    RBaM::parameter(name = "alpha0_1", init = 35),
    RBaM::parameter(name = "alpha0_2", init = 55),
    RBaM::parameter(name = "alpha0_3", init = 30)
)

if (dim(Z)[2] != length(prior)) stop(paste0("Number of columns of Z (", dim(Z_matrix_block)[2], ") must be equal to the length(prior) = ", length(prior)))

Kflood_SR[[1]] <- list(
    KP = KP[[1]],
    reach = reaches[[1]],
    Z = Z,
    prior = prior
)

# SR2
Kflood_Z_properties_MC_SR_2 <- list(
    function_SR = getCovariate_Legendre,
    max_polynomial_degree = 2
)

Z <- Kflood_Z_properties_MC_SR_2[[1]](
    max_polynomial_degree = Kflood_Z_properties_MC_SR_2[[2]],
    covariate_discretization = KP[[2]]
)

prior <- list(
    RBaM::parameter(name = "alpha0", init = 30),
    RBaM::parameter(name = "alpha1", init = 0),
    RBaM::parameter(name = "alpha2", init = 0)
)

if (dim(Z)[2] != length(prior)) stop(paste0("Number of columns of Z (", dim(Z_matrix_block)[2], ") must be equal to the length(prior) = ", length(prior)))

Kflood_SR[[2]] <- list(
    KP = KP[[2]],
    reach = reaches[[2]],
    Z = Z,
    prior = prior
)

# SR3
Kflood_Z_properties_TB1_SR_1 <- list(
    function_SR = getCovariate_Legendre,
    max_polynomial_degree = 1
)

Z <- Kflood_Z_properties_TB1_SR_1[[1]](
    max_polynomial_degree = Kflood_Z_properties_TB1_SR_1[[2]],
    covariate_discretization = KP[[3]]
)

prior <- list(
    RBaM::parameter(name = "alpha0", init = 30),
    RBaM::parameter(name = "alpha1", init = 0)
)

if (dim(Z)[2] != length(prior)) stop(paste0("Number of columns of Z (", dim(Z_matrix_block)[2], ") must be equal to the length(prior) = ", length(prior)))

Kflood_SR[[3]] <- list(
    KP = KP[[3]],
    reach = reaches[[3]],
    Z = Z,
    prior = prior
)

# SR4
Kflood_Z_properties_TB2_SR_1 <- list(
    function_SR = getCovariate_Legendre,
    max_polynomial_degree = 0
)

Z <- Kflood_Z_properties_TB2_SR_1[[1]](
    max_polynomial_degree = Kflood_Z_properties_TB2_SR_1[[2]],
    covariate_discretization = KP[[3]]
)

prior <- list(
    RBaM::parameter(name = "alpha0", init = 30)
)

if (dim(Z)[2] != length(prior)) stop(paste0("Number of columns of Z (", dim(Z_matrix_block)[2], ") must be equal to the length(prior) = ", length(prior)))

Kflood_SR[[4]] <- list(
    KP = KP[[4]],
    reach = reaches[[4]],
    Z = Z,
    prior = prior
)
#######################################
# End Kflood
#######################################

#######################################
# Structural error settings
#######################################


remant_error_list <- list(
    remnantErrorModel(
        fname = "Config_RemnantSigma.txt",
        funk = "Constant",
        par = list(parameter(
            name = "intercept",
            init = 0.0005,
            prior.dist = "FlatPrior",
        ))
    ),
    remnantErrorModel(
        fname = "Config_RemnantSigma2.txt",
        funk = "Constant",
        par = list(parameter(
            name = "intercept",
            init = 10,
            prior.dist = "LogNormal",
            prior.par = c(log(10), 0.2)
        ))
    ),
    remnantErrorModel(
        fname = "Config_RemnantSigma3.txt",
        funk = "Constant",
        par = list(parameter(
            name = "intercept",
            init = 10,
            prior.dist = "LogNormal",
            prior.par = c(log(10), 0.2)
        ))
    ),
    remnantErrorModel(
        fname = "Config_RemnantSigma4.txt",
        funk = "Constant",
        par = list(parameter(
            name = "intercept",
            init = 15,
            prior.dist = "LogNormal",
            prior.par = c(log(15), 0.2)
        ))
    ),
    remnantErrorModel(
        fname = "Config_RemnantSigma5.txt",
        funk = "Constant",
        par = list(parameter(
            name = "intercept",
            init = 15,
            prior.dist = "LogNormal",
            prior.par = c(log(15), 0.2)
        ))
    )
)



############################################
# End Calibration settings
############################################
