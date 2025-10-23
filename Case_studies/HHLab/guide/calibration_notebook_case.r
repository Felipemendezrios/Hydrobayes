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

last <- function(data) {
    tail(data, n = 1)
}

# Interpolation function passing by some specific nodes:
interpolation_specific_points <- function(total_points = 100,
                                          specific_nodes) {
    if (total_points <= length(specific_nodes)) stop("Total points defined is lower or equal to the specific points required. Please either increase the number of total_points or avoid to used interpolation function")

    # Order specific_nodes
    specific_nodes <- sort(specific_nodes)

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
    specific_seq <-
        c(
            specific_seq_temp,
            last(specific_nodes)
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
        # Include left and exclude right
        for (i in 1:nrow(reach_KP_boundaries)) {
            if (i != nrow(reach_KP_boundaries)) {
                idx_position_reaches <- which(
                    grid >= reach_KP_boundaries$KP_start[i] &
                        grid < reach_KP_boundaries$KP_end[i]
                )
            } else {
                idx_position_reaches <- which(
                    grid >= reach_KP_boundaries$KP_start[i] &
                        grid <= reach_KP_boundaries$KP_end[i]
                )
            }

            reaches[idx_position_reaches] <- reach_KP_boundaries$reach[i]
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
    id_shiftSReaches <- which(duplicated(KP_grid) | duplicated(KP_grid, fromLast = TRUE))

    if (length(id_shiftSReaches) != 0) {
        # Identify which of these are the second (or later) occurrences
        is_second_occurrence <- duplicated(KP_grid[id_shiftSReaches])
        # If TRUE, matrix should be corrected to consider shift point coincide with switch of reach
        # Correct the matrix for the second (or later) occurrences
        for (i in which(is_second_occurrence)) {
            id <- id_shiftSReaches[i]
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

RUGFile_constructor_SR_correction <- function(Key_Info_XR_MR) {
    # Calculate the length of RUGFile using Key_Info_XR_MR
    length_RUGFile <- sum(
        sapply(Key_Info_XR_MR, function(id) {
            length(id$reach)
        })
    )
    # Initialize RUGFile_first_part as a data frame with 5 columns
    RUGFile_first_part <- data.frame(
        matrix(
            data = 0,
            nrow = length_RUGFile,
            ncol = 5 # RUGFILE columns : reach, kp_start, kp_end, kmin,Kflood
        )
    )
    colnames(RUGFile_first_part) <- c("reach", "KP_start", "KP_end", "Kmin", "Kflood")
    i1 <- 1
    for (id_Key in Key_Info_XR_MR) {
        i2 <- i1 + length(id_Key$KP_grid)
        RUGFile_first_part[i1:i2, c(1, 2, 3)] <- data.frame(
            id_reach = id_Key$reach[-length(id_Key$reach)],
            KP_start = id_Key$KP_grid[-length(id_Key$KP_grid)],
            KP_end = id_Key$KP_grid[-1]
        )
        i1 <- i2 + 1
    }
    return(RUGFile = RUGFile_first_part)
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


SR_constructor <- function(SR_key_HM,
                           Key_Info_XR_MR) {
    # Initialize structure by fixing layer 1 : XR
    # Layer 2 : Kmin or Kmoy is handle outside this function
    SR <- vector(mode = "list", length = length(SR_key_HM))
    names(SR) <- names(SR_key_HM)

    # Extract only KP boundaries points keeping the same structure of XR
    SR_KP_boundaries_structure <- lapply(SR_key_HM, function(element) {
        lapply(element, function(SR) {
            SR$KP_boundaries_points
        })
    })

    # Loop through XR
    for (id_XR in seq_along(SR_key_HM)) {
        # Get KP boundaries of all SR at a XR
        SR_KP_boundaries_list <- SR_KP_boundaries_structure[[id_XR]]
        RUG_min_boundary_KP <- min(Key_Info_XR_MR[[id_XR]]$RUGFile$KP_start)
        RUG_max_boundary_KP <- max(Key_Info_XR_MR[[id_XR]]$RUGFile$KP_end)

        # Check if KP boundaries of each SR respect the KP boundaries of XR
        if (min(unlist(SR_KP_boundaries_list)) != RUG_min_boundary_KP ||
            max(unlist(SR_KP_boundaries_list)) != RUG_max_boundary_KP) {
            return(stop(sprintf(
                "KP boundary mismatch: data min/max (%s, %s) vs HM min/max (%s, %s)",
                RUG_min_boundary_KP,
                RUG_max_boundary_KP,
                min(unlist(SR_KP_boundaries_list)),
                max(unlist(SR_KP_boundaries_list))
            )))
        }
        data_KP <- Key_Info_XR_MR[[id_XR]]$KP_grid
        data_reaches <- Key_Info_XR_MR[[id_XR]]$reach

        # Adding layer 3 : declare all SRs into the structure
        SR[[id_XR]] <- vector("list", length(SR_KP_boundaries_list))
        names(SR[[id_XR]]) <- names(SR_KP_boundaries_list)
        # Loop through each SR
        for (id_SR in seq_along(SR_KP_boundaries_list)) {
            # Adding layer 4: assign properties of each SR
            # Starting with KP and reach
            boundaries <- SR_KP_boundaries_list[[id_SR]]

            # Include left and exclude right
            if (id_SR != length(SR_KP_boundaries_list)) {
                position_match <- which(data_KP >= boundaries[1] &
                    data_KP < boundaries[2])
            } else {
                position_match <- which(data_KP >= boundaries[1] &
                    data_KP <= boundaries[2])
            }

            SR[[id_XR]][[id_SR]]$KP <- data_KP[position_match]
            SR[[id_XR]][[id_SR]]$reach <- data_reaches[position_match]

            # Add Z and priors
            if (identical(
                SR_key_HM[[id_XR]][[id_SR]]$function_SR,
                getCovariate_piecewise
            )) {
                if (!("shiftPoints" %in% names(SR_key_HM[[id_XR]][[id_SR]]))) {
                    return(stop(paste0(
                        "To apply getCovariate_piecewise function, shiftPoints must be passed as argument. Please check the XR = ", names(SR_key_HM)[[id_XR]], ", SR = ", names(SR_key_HM[[id_XR]])[[id_SR]]
                    )))
                }
                SR[[id_XR]][[id_SR]]$Z <- getCovariate_piecewise(
                    shiftPoints = SR_key_HM[[id_XR]][[id_SR]]$shiftPoints,
                    KP_grid = SR[[id_XR]][[id_SR]]$KP
                )
            } else if (identical(
                SR_key_HM[[id_XR]][[id_SR]]$function_SR,
                getCovariate_Legendre
            )) {
                if (!("max_polynomial_degree" %in% names(SR_key_HM[[id_XR]][[id_SR]]))) {
                    return(stop(paste0(
                        "To apply getCovariate_Legendre function, max_polynomial_degree must be passed as argument. Please check the XR = ", names(SR_key_HM)[[id_XR]], ", SR = ", names(SR_key_HM[[id_XR]])[[id_SR]]
                    )))
                }
                SR[[id_XR]][[id_SR]]$Z <- getCovariate_Legendre(
                    max_polynomial_degree = SR_key_HM[[id_XR]][[id_SR]]$max_polynomial_degree,
                    covariate_discretization = SR[[id_XR]][[id_SR]]$KP
                )
            } else {
                return(stop("function_SR given in input is not supported. Please select either getCovariate_Legendre or getCovariate_piecewise"))
            }

            SR[[id_XR]][[id_SR]]$prior <- SR_key_HM[[id_XR]][[id_SR]]$prior

            if (dim(SR[[id_XR]][[id_SR]]$Z)[2] != length(SR[[id_XR]][[id_SR]]$prior)) {
                stop(paste0(
                    "Error identified in XR = ", names(SR_key_HM)[[id_XR]], ", SR = ", names(SR_key_HM[[id_XR]])[[id_SR]], ". Number of columns of Z (", dim(SR[[id_XR]][[id_SR]]$Z)[2], ") must be equal to the length(prior) = ", length(SR[[id_XR]][[id_SR]]$prior)
                ))
            }
        }
    }
    return(SR)
}

############################################
# End source functions
############################################


############################################
# Modules
############################################

############################################
# Module 1 : hydraulic model (HM) environment
############################################
# Input information about reaches and boundaries KP of the HM.
# MR : Model reach interpreted by the HM.
# This input information must be coherent with HM specifications.
Input_MR <- data.frame(
    reach = c(1, 2, 3, 4),
    KP_start = c(0, 5, 10, 25),
    KP_end = c(10, 15, 20, 35)
)


############################################
# Module 2 : BaM environment
############################################

# First layer: XR
############################################
# All reaches in the model (MR) defined within the same coordinate reference

Input_XR <- list(
    main_channel = c(1, 3),
    tributary_1 = 2,
    tributary_2 = 4
)

# At this level, all information is already available to link BaM and HM
# Objective: common information to share in BaM and HM environment
# This step is so important to ensure the transmission of information

# Define the grid of interpolation by XR
Key_Info_XR_MR <- vector(mode = "list", length = length(Input_XR))
names(Key_Info_XR_MR) <- names(Input_XR)
for (i in seq_along(Key_Info_XR_MR)) {
    # Get information from MH by XR
    mask_MR_by_XR <- Input_MR[which(Input_MR$reach %in% Input_XR[[i]]), ]

    # Get nodes of MR
    MR_nodes <- unique(c(
        mask_MR_by_XR$KP_start,
        mask_MR_by_XR$KP_end
    ))

    Key_Info_XR_MR[[i]]$MR_nodes <- MR_nodes

    # Grid from interpolation passing by MR_nodes
    # KP_MR_nodes is defined by XR
    KP_MR_nodes <- interpolation_specific_points(
        total_points = 100,
        specific_nodes = MR_nodes
    )
    # Assign reach at the KP_MR_nodes
    reach_MR_nodes <- assign_reach_from_a_grid(
        reach_KP_boundaries = mask_MR_by_XR,
        grid = KP_MR_nodes
    )
    # These information allow build RUGFile
    RUGFile <- data.frame(
        id_reach = reach_MR_nodes[-length(reach_MR_nodes)],
        KP_start = KP_MR_nodes[-length(KP_MR_nodes)],
        KP_end = KP_MR_nodes[-1],
        Kmin = 0,
        Kflood = 0
    )
    Key_Info_XR_MR[[i]]$RUGFile <- RUGFile
    # Now, the middle value of each interval is taken to get the KP_grid, which will be used for all the calculations

    # Calculate the middle value for each row
    Key_Info_XR_MR[[i]]$KP_grid <- (RUGFile$KP_start + RUGFile$KP_end) / 2
    Key_Info_XR_MR[[i]]$reach <- RUGFile$id_reach
}


########################
# RUGFile
########################
# Extract all RUGFile data frames from the main list
all_RUGFiles <- lapply(Key_Info_XR_MR, function(channel) {
    channel$RUGFile
})

# Row-bind all RUGFile data frames into a single data frame
RUGFile <- do.call(rbind, all_RUGFiles)
# Remove row names
rownames(RUGFile) <- NULL
#######################################
# End RUGFile
#######################################


# Second layer: Parameters for vertical estimation of friction
####################################################################

# In MAGE context:
# Constant piecewise function : Kmin and Kflood

# In Dashflow context:
# K = alpha * (H) ^ beta
# Exponential function : so, alpha and beta are the parameters


############################################
# Kmin environment (encapsulated in Kmin_SR)
############################################

# Third layer: Set all SRs to each XR by Kmin or Kflood
# Fourth layer: Set properties of each SR in terms of KP and spatialisation function
####################################################################
# Key to relate all SR by XR
Input_Kmin_Key_SR_MR <- list(
    # XR information
    main_channel =
        list(
            SR1 = list(
                # KP boundary points
                KP_boundaries_points = c(0, 20),
                # Function to apply at this SR
                function_SR = getCovariate_piecewise,
                # Arguments of this SR
                shiftPoints = c(5, 10, 14, 18),
                prior = list(
                    RBaM::parameter(name = "MC_Kmin_SR1_a0_1", init = 35),
                    RBaM::parameter(name = "MC_Kmin_SR1_a0_2", init = 55),
                    RBaM::parameter(name = "MC_Kmin_SR1_a0_3", init = 30),
                    RBaM::parameter(name = "MC_Kmin_SR1_a0_4", init = 20),
                    RBaM::parameter(name = "MC_Kmin_SR1_a0_5", init = 25)
                )
            )
        ),
    tributary_1 = list(
        SR1 = list(
            KP_boundaries_points = c(5, 15),
            function_SR = getCovariate_Legendre,
            max_polynomial_degree = 2,
            prior = list(
                RBaM::parameter(name = "TB1_Kmin_SR1_a0", init = 30),
                RBaM::parameter(name = "TB1_Kmin_SR1_a1", init = 0),
                RBaM::parameter(name = "TB1_Kmin_SR1_a2", init = 0)
            )
        )
    ),
    tributary_2 = list(
        SR1 = list(
            KP_boundaries_points = c(25, 35),
            function_SR = getCovariate_Legendre,
            max_polynomial_degree = 1,
            prior = list(
                RBaM::parameter(name = "TB2_Kmin_SR1_a0", init = 30),
                RBaM::parameter(name = "TB2_Kmin_SR1_a1", init = 0)
            )
        )
    )
)

# Check if size is respected between ID and Kmin
if (length(Input_XR) != length(Input_Kmin_Key_SR_MR)) stop("Size must be equal between Input_Kmin_Key_SR_MR and Input_XR")

# Assign properties of each SR in XR structure
Kmin_SR <- SR_constructor(
    SR_key_HM = Input_Kmin_Key_SR_MR,
    Key_Info_XR_MR = Key_Info_XR_MR
)

# Get Z file for all XR in Kmin
Z_all_Kmin_SR <- lapply(Kmin_SR, function(channel) {
    lapply(channel, function(SR) {
        SR$Z
    })
})
# Flatten the nested list of matrices
flat_Kmin_SR <- unlist(Z_all_Kmin_SR, recursive = FALSE)

# Pass the flattened list to block_diagonal_matrix
zFileKmin <- do.call(block_diagonal_matrix, flat_Kmin_SR)

# Check size between RUGFile and Z_file
if (nrow(RUGFile) != nrow(zFileKmin)) stop("RugFile must have the same size as Z file (spatialisation)")

############################################
# End Kmin environment
############################################

############################################
# Kflood environment (encapsulated in Kflood_SR)
############################################

# Third layer: Set all SRs to each XR by Kflood
# Fourth layer: Set properties of each SR in terms of KP and spatialisation function
####################################################################
# Key to relate all SR by XR
Input_Kflood_Key_SR_MR <- list(
    # XR information
    main_channel =
        list(
            SR1 = list(
                # KP boundary points
                KP_boundaries_points = c(0, 10),
                # Function to apply at this SR
                function_SR = getCovariate_piecewise,
                # Arguments of this SR
                shiftPoints = c(2.5, 7),
                prior = list(
                    RBaM::parameter(name = "MC_Kflood_SR1_a0_1", init = 35),
                    RBaM::parameter(name = "MC_Kflood_SR1_a0_2", init = 55),
                    RBaM::parameter(name = "MC_Kflood_SR1_a0_3", init = 30)
                )
            ),
            SR2 = list(
                KP_boundaries_points = c(10, 20),
                function_SR = getCovariate_Legendre,
                max_polynomial_degree = 2,
                prior = list(
                    RBaM::parameter(name = "MC_Kflood_SR2_a0", init = 35),
                    RBaM::parameter(name = "MC_Kflood_SR2_a1", init = 0),
                    RBaM::parameter(name = "MC_Kflood_SR2_a2", init = 0)
                )
            )
        ),
    tributary_1 = list(
        SR1 = list(
            KP_boundaries_points = c(5, 15),
            function_SR = getCovariate_piecewise,
            shiftPoints = c(8),
            prior = list(
                RBaM::parameter(name = "TB1_Kflood_SR1_a0_1", init = 30),
                RBaM::parameter(name = "TB1_Kflood_SR1_a0_2", init = 20)
            )
        )
    ),
    tributary_2 = list(
        SR1 = list(
            KP_boundaries_points = c(25, 35),
            function_SR = getCovariate_Legendre,
            max_polynomial_degree = 0,
            prior = list(
                RBaM::parameter(name = "TB2_Kflood_SR1_a0", init = 30)
            )
        )
    )
)

# Check if size is respected between ID and Kmin
if (length(Input_XR) != length(Input_Kflood_Key_SR_MR)) stop("Size must be equal between Input_Kflood_Key_SR_MR and Input_XR")

# Assign properties of each SR in XR structure
Kflood_SR <- SR_constructor(
    SR_key_HM = Input_Kflood_Key_SR_MR,
    Key_Info_XR_MR = Key_Info_XR_MR
)

# Get Z file for all XR in Kflood
Z_all_Kflood_SR <- lapply(Kflood_SR, function(channel) {
    lapply(channel, function(SR) {
        SR$Z
    })
})
# Flatten the nested list of matrices
flat_Kflood_SR <- unlist(Z_all_Kflood_SR, recursive = FALSE)

# Pass the flattened list to block_diagonal_matrix
zFileKflood <- do.call(block_diagonal_matrix, flat_Kflood_SR)

# Check size between RUGFile and Z_file
if (nrow(RUGFile) != nrow(zFileKflood)) stop("RugFile must have the same size as Z file (spatialisation)")

############################################
# End Kflood environment
############################################

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

# Extract all third-level sublists from the nested list
third_level_sublists <- lapply(Kmin_SR, function(channel) {
    lapply(channel, function(SR) {
        # Extract the third-level sublists (priors)
        SR$prior
    })
})

# Flatten the list to remove the second and first levels
Kmin_prior <- unlist(third_level_sublists, recursive = FALSE)
# Extract all third-level sublists from the nested list
third_level_sublists <- lapply(Kflood_SR, function(channel) {
    lapply(channel, function(SR) {
        # Extract the third-level sublists (priors)
        SR$prior
    })
})

# Flatten the list to remove the second and first levels
Kflood_prior <- unlist(third_level_sublists, recursive = FALSE)

theta_param <- c(Kmin_prior, Kflood_prior)
MAGE_executable <- "/home/famendezrios/Documents/Softwares/pamhyr2/mage8/mage"

xtra <- xtraModelInfo(
    fname = "Config_setup.txt",
    object = list(
        exeFile = MAGE_executable,
        version = "8",
        mageDir = mageDir,
        repFile = paste0(mage_projet_name, ".REP"),
        zKmin = matrix_zFileKmin,
        zFileKmin = zFileKmin,
        doExpKmin = FALSE,
        zKmoy = matrix_zFileKmoy,
        zFileKmoy = zFileKmoy,
        doExpKmoy = FALSE
    )
)

mod_polynomials <- model(
    ID = "MAGE_ZQV",
    nX = 4,
    nY = 5, ,
    par = theta_param,
    xtra = xtra
)

runModel(
    workspace = tempdir(),
    mod = mod_polynomials,
    X = X,
    stout = NULL
)
