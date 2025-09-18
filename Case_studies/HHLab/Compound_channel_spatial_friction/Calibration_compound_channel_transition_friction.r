rm(list = ls())
graphics.off()

# function_list <- list.files("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/Functions/", full.names = TRUE)
# for (i in function_list) {
#     source(i)
# }

#############################################
# Load libraries
#############################################

library(RBaM)
library(dplyr)
library(patchwork)
library(tidyr)
library(ggplot2)
library(stringr)
#############################################
# End load libraries
#############################################
# Load functions
# function_list <- list.files("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/Functions/", full.names = TRUE)
# for (i in function_list) {
#     source(i)
# }



#############################################
# Load Functions
#############################################
# Assuming Kflood_SR is your nested list
extract_priors <- function(nested_list) {
    priors <- list()

    # Recursive function to extract 'prior' elements
    extract <- function(x) {
        if (is.list(x)) {
            if ("prior" %in% names(x)) {
                priors <<- c(priors, list(x$prior))
            } else {
                lapply(x, extract)
            }
        }
    }

    # Apply the recursive function
    extract(nested_list)

    # Flatten the list of priors
    flat_priors <- unlist(priors, recursive = FALSE)

    return(flat_priors)
}

get_init_prior <- function(parameter, FIX_dist = FALSE) {
    # Identify if parameter is remnantErrorModel class
    logical_test <- class(parameter[[1]]) == "remnantErrorModel"
    if (logical_test) { # if remnantErrorModel, a list is needed
        init_priors <- list()
    } else { # if not a vector is needed
        init_priors <- numeric(0)
    }

    counter_gamma <- 1
    for (i in parameter) {
        # Handle if parameters is remnantErrorModel
        if (logical_test) {
            param <- i$par
            number_var_error_model <- seq_along(param)

            for (local_counter in number_var_error_model) {
                if (FIX_dist) {
                    init_priors[[counter_gamma]] <- param[[local_counter]]$init
                    counter_gamma <- counter_gamma + 1
                } else {
                    if (param[[local_counter]]$prior$dist != "FIX") {
                        init_priors[[counter_gamma]] <- param[[local_counter]]$init
                        counter_gamma <- counter_gamma + 1
                    }
                }
            }
        } else {
            # Handle if parameters comes from theta
            param <- i
            if (FIX_dist) {
                init_priors <- c(init_priors, param$init)
            } else {
                if (param$prior$dist != "FIX") {
                    init_priors <- c(init_priors, param$init)
                }
            }
        }
    }
    return(init_priors)
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
# Module 1 : set directory paths
############################################

# Set directory
dir_workspace <- here::here()

# Common for observation and calibration folders
all_experiments <- c(
    "piecewise_function"
)
# Experiments input data to be used during calibration setting
all_cal_case <- c(
    "Kmin_n_deg_0_Kflood_PW_function_nK_2.r"
)

# Folder related to the observations (careful with the order!)
all_events <- c(
    "CWMQ18",
    "CWMQ12"
)

dir_exe_BaM <- "/home/famendezrios/Documents/Git/BaM/makefile/"
MAGE_executable <- "/home/famendezrios/Documents/Softwares/pamhyr2/mage8/mage"
file_main_path <- "/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/Case_studies/HHLab/Compound_channel_spatial_friction/Calibration_experiments"
mage_projet_name <- "CWMQ"
############################################
# End module 1: set directories paths
############################################

############################################
# Module 2: calibration data
############################################
# Processed data
load("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/processed_data/HHLab/compound_transition_friction/data_HHLab_uniform_case.RData")

# Observations data input:
# Measurements for calibration by event!
# Ensure that order should be in coherence with mage model set up!
Cal_measures <-
    list(
        CWMQ18 =
            data.frame(
                event = 1,
                reach = 1,
                x = c(2.225, 4.310, 5.270, 6.240, 8.315, 9.115, 9.250, 10.250, 11.250, 12.250, 13.250),
                t = rep(3420, 11)
            ),
        CWMQ12 =
            data.frame(
                event = 2,
                reach = 1,
                x = c(2.225, 2.975, 3.345, 4.095, 5.060, 6.240, 6.990, 8.030, 8.795, 9.065, 9.115, 9.545, 10.030, 11.370, 12.120, 12.325, 13.290),
                t = rep(3420, 17)
            )
    )


# All available data
WSE_data_temp <- all_data_calibration$WSE

WSE_data_temp <- data.frame(WSE_data_temp %>%
    group_by(ID_experiment, x) %>%
    summarise(
        z_mean = mean(z_mean),
        Yu = mean(Yu),
        z_riverbed = mean(z_riverbed),
        .groups = "drop"
    ))

WSE_data_temp <- WSE_data_temp[str_detect(WSE_data_temp$ID_experiment, pattern = mage_projet_name), ]

# This part help to handle when observed data is presented as a list and user need to specify the the order of the events with the same name for matching with all_all_events, where the order is important

# Only giving the name keeping the data order and creating a link with the names of the event and the Cal_measures defined previously
match_case_event <- data.frame(
    idx_data = unique(levels(factor(WSE_data_temp$ID_experiment))),
    idx_event = unique(levels(factor(WSE_data_temp$ID_experiment)))
)

Link_x_t_ind_event <- lapply(names(Cal_measures), function(event) {
    # Find the row in match_case_event where idx_event matches the current event
    row_idx <- which(match_case_event$idx_event == event)

    list(
        idx_data = match_case_event$idx_data[row_idx],
        idx_event = event,
        Cal_measure = Cal_measures[[event]]
    )
})

# # Initialize an empty list to store the results
observed_data <- data.frame()
X <- c()

# Loop over each element in Link_x_t_ind_event
for (i in seq_along(Link_x_t_ind_event)) {
    idx_data <- Link_x_t_ind_event[[i]]$idx_data
    cal_x <- Link_x_t_ind_event[[i]]$Cal_measure$x

    # Get the corresponding data frame from WSE_data_temp
    # wse_df <- WSE_data_temp[[idx_data]]
    wse_df <- WSE_data_temp[WSE_data_temp$ID_experiment == idx_data, ]

    # Filter rows where x is in cal_x
    observed_data <- rbind(
        observed_data,
        data.frame(
            event = Link_x_t_ind_event[[i]]$Cal_measure$event,
            wse_df[wse_df$x %in% cal_x, c("z_riverbed", "z_mean", "Yu")]
        )
    )
    # Extract the Cal_measure data frame
    X <- rbind(
        X,
        Link_x_t_ind_event[[i]]$Cal_measure
    )
}

if (!all(observed_data$event == X$event)) stop("The order of the event are not the same between the X data frame and observed_data data frame")

# Assign manually uncertainty values in meters
do_manual_uncertainty <- TRUE
if (do_manual_uncertainty) {
    uncertainty_in_observation <- rep(0.005, length(observed_data$Yu))
} else {
    uncertainty_in_observation <- observed_data$Yu
}

Y <- data.frame(WSE = observed_data$z_mean)
Y$Discharge <- -9999
Y$Velocity <- -9999
Y$Y_Kmin <- -9999
Y$Y_Kmoy <- -9999

Yu <- data.frame(Yu_z = uncertainty_in_observation)
Yu$Yu_Q <- -9999
Yu$Yu_V <- -9999
Yu$Yu_Kmin <- -9999
Yu$Yu_Kmoy <- -9999


## Plot
CalData_plot_df <- data.frame(
    x = X$x,
    z_mean = Y$WSE * 1000,
    Yu = Yu$Yu_z * 1000,
    z_riverbed = observed_data$z_riverbed * 1000,
    ID = factor(X$event)
)

CalData_plot_export <- CalData_plot(
    data = CalData_plot_df,
    scales_free = ,
    y_label = "Water surface elevation (mm)",
    title_label = "WSE observations",
    col_label = "Events",
    plot_water_depth = FALSE,
    wrap = FALSE
)
############################################
# End module 2: calibration data
############################################

############################################
# Module 3 : hydraulic model (HM) environment
############################################
# Input information about reaches and boundaries KP of the HM.
# MR : Model reach interpreted by the HM.
# This input information must be coherent with HM specifications.
Input_MR <- data.frame(
    reach = c(1),
    KP_start = c(0),
    KP_end = c(18)
)
############################################
# End module 3 : hydraulic model (HM) environment
############################################


############################################
# Module 4: calibration setting
############################################
# Structural error information
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

# Setting jump standard deviation for MCMC sampling
jump_MCMC_theta_param_user <- 8
jump_MCMC_error_model_user <- 0.001
threshold_jump_MCMC_error_model <- 0.5
prior_error_model <- get_init_prior(remant_error_list)


# Logical to trigger some functionalities
do_calibration <- FALSE

############################################
# End module 4: calibration setting
############################################


############################################
# Module 5: Calculations
############################################
mod_polynomials <- list()

# for (Experiment_id in all_experiments) {

# }

Experiment_id <- all_experiments

path_experiment <- file.path(file_main_path, Experiment_id)
counter_model <- 1
if (!dir.exists(path_experiment)) {
    dir.create(path_experiment)
}

ggsave(file.path(path_experiment, "CalData_plot.png"),
    plot = CalData_plot_export,
    width = 20,
    height = 20,
    units = "cm"
)

wood_first <- mage_projet_name == "CWMQ"

for (id_cal_case in 1:length(all_cal_case)) {
    # Source the input data for the experiments
    path_Experiment_Input_Data <- file.path(file_main_path, "Experiments_Input_Data", all_cal_case[id_cal_case])
    if (!file.exists(path_Experiment_Input_Data)) {
        stop(paste0(
            "The experiment input data named : '",
            all_cal_case[id_cal_case],
            "' does not exist"
        ))
    }
    source(path_Experiment_Input_Data)
    path_polynomial <- file.path(
        path_experiment,
        sub("\\.r$", "", all_cal_case[id_cal_case])
    )
    # Path to save results
    workspace_user <- file.path(path_polynomial, "BaM")
    path_post_traitement <- file.path(workspace_user, "post_traitement")
    path_post_traitement_data <- file.path(path_post_traitement, "RData")
    if (!dir.exists(path_polynomial)) {
        dir.create(path_polynomial)
    } else {
        if (do_calibration) {
            file.remove(list.files(path_polynomial, full.names = TRUE, recursive = TRUE))
        }
    }
    if (!dir.exists(workspace_user)) {
        dir.create(workspace_user)
    } else {
        if (do_calibration) {
            file.remove(list.files(workspace_user, full.names = TRUE, recursive = TRUE))
        }
    }
    if (!dir.exists(path_post_traitement)) {
        dir.create(path_post_traitement)
    }

    if (!dir.exists(path_post_traitement_data)) {
        dir.create(path_post_traitement_data)
    }


    ################## To be continued







    if (length(Kmin_prior) != (Nb_reaches_estimation_Kmin * max((n_degree_Kmin + 1)))) stop(paste0("More prior information (", length(Kmin_prior), ") than the number of reaches for estimation (", Nb_reaches_estimation_Kmin, ")"))
    if (length(Kflood_prior) != (Nb_reaches_estimation_Kflood * (n_degree_Kflood + 1))) stop(paste0("More prior information (", length(Kflood_prior), ") than the number of reaches for estimation (", Nb_reaches_estimation_Kflood, ")"))

    theta_param <- c(Kmin_prior, Kflood_prior)

    # Write Legendre vector:
    write.table(
        grid_covariant_discretized_Kmin,
        file = file.path(workspace_user, "legendre_covariate_non_normalized_Kmin.txt"), row.names = F
    )
    # Write Legendre vector:
    write.table(
        grid_covariant_discretized_Kflood,
        file = file.path(workspace_user, "legendre_covariate_non_normalized_Kflood.txt"), row.names = F
    )


    mageDir <- c(paste0(path_polynomial, "/model_mage/", all_events, "/"))

    RUGFiles <- paste0(mageDir, mage_projet_name, ".RUG")

    # j is a index for multi-events
    for (j in seq_along(RUGFiles)) {
        # Calculate original covariant considering that grid_covariant_discretized is in the middle of the real covariant value except on the boundaries
        grid_covariant_meshed_Model <- RUG_KP_end_by_reach <- RUG_KP_start_by_reach <- RUG_id_reach_by_reach <- RUG_Kmin_by_reach <- RUG_Kflood_by_reach <- list()

        for (reach in 1:nrow(specific_points_Model)) {
            covariant_discretized_by_reach <- grid_covariant_discretized_real_bief$interpol_values[which(grid_covariant_discretized_real_bief$real_reach == specific_points_Model$reach[reach])]

            grid_covariant_meshed_Model[[reach]] <-
                data.frame(
                    reach = specific_points_Model$reach[reach],
                    real_values = get_covariant_meshed(
                        covariant_discretized = covariant_discretized_by_reach,
                        specific_points_Model = specific_points_Model[reach, ]
                    )
                )
            RUG_id_reach_by_reach[[reach]] <- rep(specific_points_Model$reach[reach], length(covariant_discretized_by_reach))

            RUG_KP_start_by_reach[[reach]] <- grid_covariant_meshed_Model[[reach]]$real_values[-length(grid_covariant_meshed_Model[[reach]]$real_values)]

            RUG_KP_end_by_reach[[reach]] <- grid_covariant_meshed_Model[[reach]]$real_values[-1]

            RUG_Kmin_by_reach[[reach]] <- rep(RUG_Kmin[reach], length(covariant_discretized_by_reach))
            RUG_Kflood_by_reach[[reach]] <- rep(RUG_Kflood[reach], length(covariant_discretized_by_reach))
        }

        RUG_KP_start <- unlist(RUG_KP_start_by_reach)
        RUG_KP_end <- unlist(RUG_KP_end_by_reach)
        RUG_id_reach <- unlist(RUG_id_reach_by_reach)
        RUG_Kmin_args <- unlist(RUG_Kmin_by_reach)
        RUG_Kflood_args <- unlist(RUG_Kflood_by_reach)


        write_RUGFile(
            RUG_path = RUGFiles[j],
            RUG_id_reach = RUG_id_reach,
            RUG_KP_start = RUG_KP_start,
            RUG_KP_end = RUG_KP_end,
            RUG_Kmin = RUG_Kmin_args,
            RUG_Kmoy = RUG_Kflood_args,
            RUG_format = "%1s%3d      %10.0f%10.0f%10.2f%10.2f"
        )
    }

    # Spatial distributed friction in the main channel
    matrix_zFileKmin <- KFile_spatial(
        reaches = reaches_user_Kmin,
        max_degree = n_degree_Kmin,
        grid_covariant_discretized = grid_covariant_discretized_Kmin
    )

    # Spatial distributed friction in the floodplain
    matrix_zFileKmoy <- KFile_spatial(
        reaches = reaches_user_Kflood,
        max_degree = n_degree_Kflood, grid_covariant_discretized = grid_covariant_discretized_Kflood
    )

    zFileKmin <- file.path(workspace_user, "Zfile_Kmin.txt")
    zFileKmoy <- file.path(workspace_user, "Zfile_Kflood.txt")

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

    mod_polynomials[[counter_model]] <- model(
        ID = "MAGE_ZQV",
        nX = 4,
        nY = 5, ,
        par = theta_param,
        xtra = xtra
    )
    mod <- mod_polynomials[[counter_model]]
    # Re-write RUG file : keep the position of the cross sections as original

    prior_theta_param <- get_init_prior(theta_param)
    prior_all_dist_theta_param <- get_init_prior(theta_param, FIX_dist = TRUE)

    data <- dataset(X = X, Y = Y, Yu = Yu, data.dir = file.path(workspace_user))

    jump_MCMC_theta_param <- ifelse(prior_theta_param != 0,
        prior_theta_param[(prior_theta_param != 0)] * 0.1,
        jump_MCMC_theta_param_user
    )

    jump_MCMC_error_model <- list()
    for (ind in 1:length(prior_error_model)) {
        jump_MCMC_error_model[[ind]] <- ifelse(prior_error_model[[ind]] > threshold_jump_MCMC_error_model,
            prior_error_model[[ind]][(prior_error_model[[ind]] != 0)] * 0.1,
            jump_MCMC_error_model_user
        )
    }

    mcmcOptions_user <- mcmcOptions(
        nCycles = 50,
        nAdapt = 50,
        manualMode = TRUE,
        thetaStd = jump_MCMC_theta_param,
        gammaStd = jump_MCMC_error_model
    )
    mcmcCooking_user <- mcmcCooking(
        burn = 0.2,
        nSlim = 1
    )
    mcmcSummary_user <- mcmcSummary(xtendedMCMC.fname = "Results_xtendedMCMC.txt")

    # Save all data used during calibration for prediction
    save(mod, data, remant_error_list,
        mcmcOptions_user, mcmcCooking_user,
        mcmcSummary_user, workspace_user,
        file = file.path(path_post_traitement_data, "BaM_objects.RData")
    )
    BaM(
        mod = mod,
        data = data,
        remnant = remant_error_list,
        mcmc = mcmcOptions_user,
        cook = mcmcCooking_user,
        summary = mcmcSummary_user,
        residuals = residualOptions(),
        pred = NULL,
        doCalib = TRUE,
        doPred = FALSE,
        na.value = -9999,
        run = FALSE,
        preClean = FALSE,
        workspace = workspace_user,
        # dir.exe = file.path(find.package("RBaM"), "bin"),
        # name.exe = "BaM",
        predMaster_fname = "Config_Pred_Master.txt"
    )
    counter_model <- counter_model + 1

    dir_cf <- file.path(workspace_user, "Config_BaM.txt")

    if (do_calibration) {
        system2(
            command = file.path(dir_exe_BaM, "BaM"),
            args = c("-cf", file.path(workspace_user, "Config_BaM.txt")),
            wait = FALSE
        )
    }
}
