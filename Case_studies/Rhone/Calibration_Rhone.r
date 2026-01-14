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
library(lubridate)
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

last <- function(data) {
    tail(data, n = 1)
}
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
assign_reach_from_a_grid <- function(reach_KP_boundaries, grid, Logical_decreasing) {
    grid_with_reaches <- data.frame(
        reaches = rep(NA, length(grid)),
        grid = rep(NA, length(grid))
    )

    if (nrow(reach_KP_boundaries) == 1) {
        all_set <- dplyr::between(grid, reach_KP_boundaries$KP_start, reach_KP_boundaries$KP_end)
        if (!all(all_set)) stop("reach_KP_boundaries must contain all grid")
        grid_with_reaches$reaches[which(all_set)] <- reach_KP_boundaries$reach[1]
        grid_with_reaches$grid[which(all_set)] <- grid
    } else {
        if (Logical_decreasing) {
            # Include end (right) and exclude start (left)
            for (i in 1:nrow(reach_KP_boundaries)) {
                if (i != nrow(reach_KP_boundaries)) {
                    idx_position_reaches <- which(
                        grid > reach_KP_boundaries$KP_start[i] &
                            grid <= reach_KP_boundaries$KP_end[i]
                    )
                } else {
                    idx_position_reaches <- which(
                        grid >= reach_KP_boundaries$KP_start[i] &
                            grid <= reach_KP_boundaries$KP_end[i]
                    )
                }

                grid_with_reaches$reaches[idx_position_reaches] <- reach_KP_boundaries$reach[i]
                grid_with_reaches$grid[idx_position_reaches] <- grid[idx_position_reaches]
            }
        } else {
            # Include start (left) and exclude end (right)
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


                grid_with_reaches$reaches[idx_position_reaches] <- reach_KP_boundaries$reach[i]
                grid_with_reaches$grid[idx_position_reaches] <- grid[idx_position_reaches]
            }
        }
    }

    if (any(is.na(grid_with_reaches))) {
        return(stop("Some values in the vector reaches have not been assigned"))
    }
    return(grid_with_reaches)
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
    plot_CalData <-
        ggplot(data = data, aes(
            x = x,
            y = WSE,
            col = reach_name
        )) +
        geom_line(
            aes(
                y = Z_thalweg,
                col = "Thalweg"
            )
        )

    if (any(data$Yu != 0)) {
        plot_CalData <- plot_CalData +
            geom_errorbar(aes(
                ymin = WSE - qnorm(0.975) * Yu,
                ymax = WSE + qnorm(0.975) * Yu
            ))
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
            facet_wrap(~name_event, scales = scales_free, ncol = 1)
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

            if (Key_Info_XR_MR[[id_XR]]$Logical_decreasing) {
                # If logical_decreasing is TRUE, boundaries must be reordered to the calculation
                boundaries <- sort(boundaries)
                # Include end (right) and exclude start (left)
                if (id_SR != length(SR_KP_boundaries_list)) {
                    position_match <- which(
                        data_KP > boundaries[1] &
                            data_KP <= boundaries[2]
                    )
                } else {
                    position_match <- which(
                        data_KP >= boundaries[1] &
                            data_KP <= boundaries[2]
                    )
                }
            } else {
                # Include start (left) and exclude end (right)
                if (id_SR != length(SR_KP_boundaries_list)) {
                    position_match <- which(
                        data_KP >= boundaries[1] &
                            data_KP < boundaries[2]
                    )
                } else {
                    position_match <- which(
                        data_KP >= boundaries[1] &
                            data_KP <= boundaries[2]
                    )
                }
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


write_RUGFile <- function(RUG_path,
                          RUGFile_data,
                          RUG_format) {
    # Open a .RUG file for writing
    fileConn <- file(RUG_path, "w")
    # Write the first line as a comment
    writeLines("* This file is generated by PAMHYR, please don't modify", fileConn)
    formatted_lines <- sapply(1:nrow(RUGFile_data), function(i) {
        sprintf(
            RUG_format,
            "K",
            RUGFile_data$id_reach[i],
            RUGFile_data$KP_start[i],
            RUGFile_data$KP_end[i],
            RUGFile_data$Kmin[i],
            RUGFile_data$Kflood[i]
        )
    })

    # Write the formatted lines to the file
    writeLines(formatted_lines, fileConn)


    # Close the file
    close(fileConn)
}

copy_folder <- function(source, destination) {
    # Create the destination directory if it doesn't exist
    if (!dir.exists(destination)) {
        dir.create(destination, recursive = TRUE)
    }

    # List all files and subdirectories in the source
    items <- list.files(source, full.names = TRUE, recursive = FALSE)

    for (item in items) {
        item_relative <- basename(item)
        dest_item <- file.path(destination, item_relative)

        if (dir.exists(item)) {
            # If the item is a subdirectory, recursively copy it
            copy_folder(item, dest_item)
        } else {
            # If the item is a file, copy it
            file.copy(item, dest_item)
        }
    }
}

############################################
# End source functions
############################################


############################################
# Module 1 : set directory paths
############################################

# Set directory
dir_workspace <- here::here()

# Common for observation and calibration folders: names of the experiment
Experiment_id <- c(
    "test_1"
)
# Experiments input data to be used during calibration setting
all_cal_case <- c(
    "Kmin_n_4_Kflood_fix.r"
)

# Folder related to the observations (careful with the order!)
all_events <- c(
    "AIN_90",
    "RHONE_525",
    "RHONE_750"
)

dir_exe_BaM <- "/home/famendezrios/Documents/Git/BaM/makefile/"
MAGE_executable <- "/home/famendezrios/Documents/Softwares/pamhyr2/mage8/mage"
file_main_path <- "/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/Case_studies/Rhone/Calibration_experiments"
MAGE_main_folder <- "/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/Case_studies/Rhone/model_mage/"
mage_projet_name <- "Rhone_PCH_Ain"

if (!dir.exists(MAGE_main_folder)) stop("MAGE_main_folder does not exist")
############################################
# End module 1: set directories paths
############################################

############################################
# Module 2: calibration data
############################################

sd_WSE_fixed <- 0.05 # 5 cm

# Processed data
load("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/processed_data/Rhone/all_observations_AIN.RData")
load("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/processed_data/Rhone/all_observations_Rhone.RData")
# Observations data input:
# Measurements for calibration by event!


# Traitement to each event used during calibration
# Read key information of event
key_info_event_extraction <- read.table("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/processed_data/Rhone/Boundary_conditions/Specific_boundary_conditions/Dates_start_end_extraction.csv", header = TRUE, sep = ",")

#####################################################################
# Case AIN 90
#####################################################################
init_model_date <- key_info_event_extraction %>%
    filter(id_case == "AIN_90") %>%
    pull(model_start)

init_model_date <- as.POSIXct(
    init_model_date,
    tz = "UTC",
    format = "%Y-%m-%d %H:%M:%S"
)
WSE_Ain_90_raw <- all_WSE_Ain %>%
    filter(id_case == 90) %>%
    mutate(
        id_reach_CAL = case_when(
            id_reach == "PCH_PGA" ~ 8,
            id_reach == "PGA_CAIN" ~ 4,
            id_reach == "CAIN_Confluence" ~ 5
        ),
        time_mage = difftime(time, init_model_date, units = c("secs")),
        Yu = sd_WSE_fixed,
        event = 1,
        name_event = "AIN_90"
    )
#################################
# Check if time are within the simulation period
PARFile <- readLines(file.path(MAGE_main_folder, "AIN_90", "Rhone_PCH_Ain.PAR"), n = 7)
end_par_file <- PARFile[3]
timestep_bin_par_file <- PARFile[7]

# Get the mage time in seconds
duration <- sub("^final_time\\s+", "", end_par_file)
parts <- as.numeric(strsplit(duration, ":")[[1]])

end_time_sim <- parts[1] * 86400 + # days
    parts[2] * 3600 + # hours
    parts[3] * 60 + # minutes
    parts[4] # seconds

if (nrow(WSE_Ain_90_raw %>% filter(WSE_Ain_90_raw$time_mage > end_time_sim)) != 0) stop("There are some observed data beyond the simulation time period")

# Round each value in time_mage to the nearest multiple of timestep_bin

timestep_bin <- as.numeric(sub("^timestep_bin\\s+", "", timestep_bin_par_file))
WSE_Ain_90_raw$time_mage_nearest <- round(WSE_Ain_90_raw$time_mage / timestep_bin) * timestep_bin

#################################
set.seed(2026) #  # for reproducibility
# 1. Calibration set (random 70%)
WSE_Ain_90_Cal <- WSE_Ain_90_raw %>%
    group_by(id_reach) %>%
    slice_sample(prop = 0.7) %>%
    ungroup() %>%
    mutate(set = "calibration") %>%
    arrange(KP) # Because order is upstream to downstream and KP is increasing

# If KP duplicated, take it off
if (any(duplicated(WSE_Ain_90_Cal$KP))) {
    WSE_Ain_90_Cal <- WSE_Ain_90_Cal[-which(duplicated(WSE_Ain_90_Cal$KP)), ]
}

# 2. Validation set = rest (30%)
WSE_Ain_90_Val <- WSE_Ain_90_raw %>%
    anti_join(WSE_Ain_90_Cal, by = colnames(WSE_Ain_90_raw)) %>%
    mutate(set = "validation")

# 3. Combine into a single data frame
WSE_Ain_90 <- bind_rows(WSE_Ain_90_Cal, WSE_Ain_90_Val) %>% arrange(KP)

ggplot(WSE_Ain_90, aes(x = time_mage, y = WSE, col = set)) +
    geom_point() +
    geom_errorbar(aes(
        ymin = WSE - qnorm(0.975) * Yu,
        ymax = WSE + qnorm(0.975) * Yu
    ))

ggplot(WSE_Ain_90_Cal, aes(x = time_mage, y = WSE, col = id_reach)) +
    geom_point() +
    geom_errorbar(aes(
        ymin = WSE - qnorm(0.975) * Yu,
        ymax = WSE + qnorm(0.975) * Yu
    ))


#####################################################################
# Case Rhone 525
#####################################################################
init_model_date <- key_info_event_extraction %>%
    filter(id_case == "RHONE_525") %>%
    pull(model_start)

init_model_date <- as.POSIXct(
    init_model_date,
    tz = "UTC",
    format = "%Y-%m-%d %H:%M:%S"
)
WSE_Rhone_525_raw <- all_WSE_Rhone %>%
    filter(id_case == 525) %>%
    mutate(
        id_reach_CAL = case_when(
            id_reach == "ANT_JNS" ~ 3,
            id_reach == "BOU_ANT" ~ 2,
            id_reach == "LGN_BOU" ~ 1
        ),
        time_mage = difftime(time, init_model_date, units = c("secs")),
        Yu = sd_WSE_fixed,
        event = 2,
        name_event = "RHONE_525"
    )
#################################
# Check if time are within the simulation period
PARFile <- readLines(file.path(MAGE_main_folder, "RHONE_525", "Rhone_PCH_Ain.PAR"), n = 7)
end_par_file <- PARFile[3]
timestep_bin_par_file <- PARFile[7]

# Get the mage time in seconds
duration <- sub("^final_time\\s+", "", end_par_file)
parts <- as.numeric(strsplit(duration, ":")[[1]])

end_time_sim <- parts[1] * 86400 + # days
    parts[2] * 3600 + # hours
    parts[3] * 60 + # minutes
    parts[4] # seconds

if (nrow(WSE_Rhone_525_raw %>% filter(WSE_Rhone_525_raw$time_mage > end_time_sim)) != 0) stop("There are some observed data beyond the simulation time period")

# Round each value in time_mage to the nearest multiple of timestep_bin

timestep_bin <- as.numeric(sub("^timestep_bin\\s+", "", timestep_bin_par_file))
WSE_Rhone_525_raw$time_mage_nearest <- round(WSE_Rhone_525_raw$time_mage / timestep_bin) * timestep_bin

#################################
set.seed(2026) #  # for reproducibility
# 1. Calibration set (random 70%)
WSE_Rhone_525_Cal <- WSE_Rhone_525_raw %>%
    group_by(id_reach) %>%
    slice_sample(prop = 0.7) %>%
    ungroup() %>%
    mutate(set = "calibration") %>%
    arrange(desc(KP)) # Because the order is from upstream to downstream

# If KP duplicated, take it off
if (any(duplicated(WSE_Rhone_525_Cal$KP))) {
    WSE_Rhone_525_Cal <- WSE_Rhone_525_Cal[-which(duplicated(WSE_Rhone_525_Cal$KP)), ]
}

# 2. Validation set = rest (30%)
WSE_Rhone_525_Val <- WSE_Rhone_525_raw %>%
    anti_join(WSE_Rhone_525_Cal, by = colnames(WSE_Rhone_525_raw)) %>%
    mutate(set = "validation")

# 3. Combine into a single data frame
WSE_Rhone_525 <- bind_rows(WSE_Rhone_525_Cal, WSE_Rhone_525_Val) %>% arrange(desc(KP))

ggplot(WSE_Rhone_525, aes(x = KP, y = WSE, col = set)) +
    geom_point() +
    geom_errorbar(aes(
        ymin = WSE - qnorm(0.975) * Yu,
        ymax = WSE + qnorm(0.975) * Yu
    ))


ggplot(WSE_Rhone_525_Cal, aes(x = KP, y = WSE, col = id_reach)) +
    geom_point() +
    geom_errorbar(aes(
        ymin = WSE - qnorm(0.975) * Yu,
        ymax = WSE + qnorm(0.975) * Yu
    ))


#####################################################################
# Case Rhone 750
#####################################################################
init_model_date <- key_info_event_extraction %>%
    filter(id_case == "RHONE_750") %>%
    pull(model_start)

init_model_date <- as.POSIXct(
    init_model_date,
    tz = "UTC",
    format = "%Y-%m-%d %H:%M:%S"
)
WSE_Rhone_750_raw <- all_WSE_Rhone %>%
    filter(id_case == 750) %>%
    mutate(
        id_reach_CAL = case_when(
            id_reach == "ANT_JNS" ~ 3,
            id_reach == "BOU_ANT" ~ 2,
            id_reach == "LGN_BOU" ~ 1
        ),
        time_mage = difftime(time, init_model_date, units = c("secs")),
        Yu = sd_WSE_fixed,
        event = 3,
        name_event = "RHONE_750"
    )
#################################
# Check if time are within the simulation period
PARFile <- readLines(file.path(MAGE_main_folder, "RHONE_750", "Rhone_PCH_Ain.PAR"), n = 7)
end_par_file <- PARFile[3]
timestep_bin_par_file <- PARFile[7]

# Get the mage time in seconds
duration <- sub("^final_time\\s+", "", end_par_file)
parts <- as.numeric(strsplit(duration, ":")[[1]])

end_time_sim <- parts[1] * 86400 + # days
    parts[2] * 3600 + # hours
    parts[3] * 60 + # minutes
    parts[4] # seconds

if (nrow(WSE_Rhone_750_raw %>% filter(WSE_Rhone_750_raw$time_mage > end_time_sim)) != 0) stop("There are some observed data beyond the simulation time period")

# Round each value in time_mage to the nearest multiple of timestep_bin

timestep_bin <- as.numeric(sub("^timestep_bin\\s+", "", timestep_bin_par_file))
WSE_Rhone_750_raw$time_mage_nearest <- round(WSE_Rhone_750_raw$time_mage / timestep_bin) * timestep_bin

#################################
set.seed(2026) #  # for reproducibility
# 1. Calibration set (random 70%)
WSE_Rhone_750_Cal <- WSE_Rhone_750_raw %>%
    group_by(id_reach) %>%
    slice_sample(prop = 0.7) %>%
    ungroup() %>%
    mutate(set = "calibration") %>%
    arrange(desc(KP)) # Because the order is from upstream to downstream

# If KP duplicated, take it off
if (any(duplicated(WSE_Rhone_750_Cal$KP))) {
    WSE_Rhone_750_Cal <- WSE_Rhone_750_Cal[-which(duplicated(WSE_Rhone_750_Cal$KP)), ]
}

# 2. Validation set = rest (30%)
WSE_Rhone_750_Val <- WSE_Rhone_750_raw %>%
    anti_join(WSE_Rhone_750_Cal, by = colnames(WSE_Rhone_750_raw)) %>%
    mutate(set = "validation")

# 3. Combine into a single data frame
WSE_Rhone_750 <- bind_rows(WSE_Rhone_750_Cal, WSE_Rhone_750_Val) %>% arrange(desc(KP))

ggplot(WSE_Rhone_750, aes(x = KP, y = WSE, col = set)) +
    geom_point() +
    geom_errorbar(aes(
        ymin = WSE - qnorm(0.975) * Yu,
        ymax = WSE + qnorm(0.975) * Yu
    ))


ggplot(WSE_Rhone_750_Cal, aes(x = KP, y = WSE, col = id_reach)) +
    geom_point() +
    geom_errorbar(aes(
        ymin = WSE - qnorm(0.975) * Yu,
        ymax = WSE + qnorm(0.975) * Yu
    ))


# !!!!!!!!!!!!
# Ensure that order should be in coherence with mage model set up!!!!!!!!!!!!
# !!!!!!!!!!!!

## just to know the order of calibration data
Cal_measures <-
    list(
        AIN_90 =
            data.frame(
                event = WSE_Ain_90_Cal$event,
                reach = WSE_Ain_90_Cal$id_reach_CAL,
                x = WSE_Ain_90_Cal$KP,
                t = WSE_Ain_90_Cal$time_mage
            ),
        RHONE_525 =
            data.frame(
                event = WSE_Rhone_525_Cal$event,
                reach = WSE_Rhone_525_Cal$id_reach_CAL,
                x = WSE_Rhone_525_Cal$KP,
                t = WSE_Rhone_525_Cal$time_mage
            ),
        RHONE_750 =
            data.frame(
                event = WSE_Rhone_750_Cal$event,
                reach = WSE_Rhone_750_Cal$id_reach_CAL,
                x = WSE_Rhone_750_Cal$KP,
                t = WSE_Rhone_750_Cal$time_mage
            )
    )

# Scaled values of WSE, discharge and uncertainties
factor_scaled <- 1

# All available data
WSE_data <- rbind(
    WSE_Ain_90,
    WSE_Rhone_525,
    WSE_Rhone_750
)

# All calibration data
observed_data <- rbind(
    WSE_Ain_90_Cal,
    WSE_Rhone_525_Cal,
    WSE_Rhone_750_Cal
)

X <- observed_data[, c(
    "event",
    "id_reach_CAL",
    "KP",
    "time_mage"
)] %>%
    rename(
        "reach" = "id_reach_CAL",
        "x" = "KP",
        "t" = "time_mage"
    )

# Assign manually uncertainty values in meters
do_manual_uncertainty <- TRUE
if (do_manual_uncertainty) {
    uncertainty_in_observation <- rep(sd_WSE_fixed, length(observed_data$Yu))
} else {
    # uncertainty_in_observation <- observed_data$Yu
    stop("not supported")
}

Y <- data.frame(WSE = observed_data$WSE)
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
    WSE = Y$WSE * factor_scaled,
    Yu = Yu$Yu_z * factor_scaled,
    Z_thalweg = observed_data$Z_thalweg * factor_scaled,
    ID = factor(X$event),
    reach_name = observed_data$id_reach,
    reach_number = observed_data$id_reach_CAL,
    name_event = observed_data$name_event
)

CalData_plot_export <- CalData_plot(
    data = CalData_plot_df,
    scales_free = "free",
    y_label = "Water surface elevation (m)",
    title_label = "WSE observations",
    col_label = "Events",
    plot_water_depth = FALSE,
    wrap = TRUE
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

# KP start and end must be reliable to the model. To define KP range, it must be from upstream (start) from upstream to downstream (end)

Input_MR <- data.frame(
    reach = c(1, 2, 3, 4, 5, 6, 7, 8),
    KP_start = c(55900, 36250, 34500, 37491, 41211, 0, 20, 22333),
    KP_end = c(36250, 34500, 26750, 41211, 41461, 20, 40, 37491)
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
# do_calibration <- TRUE

############################################
# End module 4: calibration setting
############################################


############################################
# Module 5: Calculations
############################################
mod_polynomials <- list()

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
    MAGE_polynomial_subfolder <- file.path(path_polynomial, "model_mage")

    if (!dir.exists(MAGE_polynomial_subfolder)) {
        dir.create(MAGE_polynomial_subfolder)
    }
    # Copy Main MAGE folder to local folder for calibration
    copy_folder(MAGE_main_folder, MAGE_polynomial_subfolder)

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


    # Write covariate discretization, available for Kmin and Kflood:
    write.table(
        covariate_grid,
        file = file.path(workspace_user, "grid_covariate_non_normalized.txt"), row.names = F
    )
    mageDir <- c(paste0(path_polynomial, "/model_mage/", all_events, "/"))

    RUGFile_paths <- paste0(mageDir, mage_projet_name, ".RUG")

    # j is a index for multi-events
    for (id_multi_event in seq_along(RUGFile_paths)) {
        write_RUGFile(
            RUG_path = RUGFile_paths[id_multi_event],
            RUGFile_data = RUGFile,
            RUG_format = "%1s%3d      %10.3f%10.3f%10.2f%10.2f"
        )
    }

    zFileKmin <- file.path(workspace_user, "Zfile_Kmin.txt")
    zFileKmoy <- file.path(workspace_user, "Zfile_Kflood.txt")

    xtra <- xtraModelInfo(
        fname = "Config_setup.txt",
        object = list(
            exeFile = MAGE_executable,
            version = "8",
            mageDir = mageDir,
            repFile = paste0(mage_projet_name, ".REP"),
            zKmin = Z_MatrixKmin,
            zFileKmin = zFileKmin,
            doExpKmin = FALSE,
            zKmoy = Z_MatrixKflood,
            zFileKmoy = zFileKmoy,
            doExpKmoy = FALSE
        )
    )

    theta_param <- c(Kmin_prior, Kflood_prior)

    mod_polynomials[[counter_model]] <- model(
        ID = "MAGE_ZQV",
        nX = 4,
        nY = 5, ,
        par = theta_param,
        xtra = xtra
    )
    mod <- mod_polynomials[[counter_model]]

    prior_theta_param <- get_init_prior(theta_param)

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



#################################
# PLOTS:
#################################
## Functions:
# Plot spatialisation of k in main channel or floodplain
K_plot <- function(
    matrix_spatialisation,
    mcmc,
    n_param_Kmin,
    n_param_Kflood,
    covariate_discretization,
    MAP_param_vector,
    main_channel = TRUE) {
    if (main_channel) {
        indx <- 1:n_param_Kmin
        label_title <- "Friction coefficient estimation in the main channel \nwith parametric uncertainty"
    } else {
        indx <- (n_param_Kmin + 1):(n_param_Kmin + n_param_Kflood)
        label_title <- "Friction coefficient estimation in the floodplain \nwith parametric uncertainty"
    }
    k_estimated_MCMC <- as.data.frame(as.matrix(matrix_spatialisation) %*% as.matrix(t(mcmc[, indx])))

    k_estimated_MCMC$KP <- covariate_discretization[[1]]

    # Convert to long format
    df_MCMC_sampling <- tidyr::pivot_longer(
        k_estimated_MCMC,
        cols = -KP,
        values_to = "Value"
    ) %>%
        select(KP, Value) %>%
        mutate(ID = "MCMC Sampling")

    k_estimated_MAP <- as.matrix(matrix_spatialisation) %*% MAP_param_vector[indx]

    df_MAP <- data.frame(
        KP = covariate_discretization[[1]],
        Value = k_estimated_MAP,
        ID = "MAP"
    )

    # Get 95% uncertainty for envelope curve : create ribbon data from MCMC
    df_envelope <- df_MCMC_sampling %>%
        filter(ID == "MCMC Sampling") %>%
        group_by(KP) %>%
        summarise(
            ymin = quantile(Value, probs = 0.025, na.rm = TRUE),
            ymax = quantile(Value, probs = 0.975, na.rm = TRUE),
            ID = "Parametric\nuncertainty", # so we can map to fill
            .groups = "drop"
        )

    K_plot <- ggplot() +
        geom_ribbon(
            data = df_envelope,
            aes(x = KP, ymin = ymin, ymax = ymax, fill = ID)
        ) +
        geom_line(
            data = df_MAP,
            aes(x = KP, y = Value, color = ID)
        ) +
        scale_fill_manual(values = c("Parametric\nuncertainty" = "pink")) +
        scale_color_manual(values = c("MAP" = "black")) +
        labs(
            title = label_title,
            x = "Streamwise position (meters)",
            y = expression("Friction coefficient (m"^
                {
                    1 / 3
                } * "/s)"),
            fill = "95% credibility\ninterval",
            color = NULL
        ) +
        theme_bw() +
        theme(
            plot.title = element_text(hjust = 0.5),
            legend.title = element_text(hjust = 0.5)
        )


    return(list(df_MAP, K_plot, df_envelope))
}

read_fortran_data <- function(file_path, col_widths_RUGFile, skip = 0) {
    # Read the file with the fixed-width format
    data <- read.fwf(file_path, widths = col_widths_RUGFile, header = FALSE, skip = skip)
    data <- data[, -3]
    colnames(data) <-
        c(
            "",
            "id_reach",
            "KP_start",
            "KP_end",
            "Kmin",
            "Kflood"
        )
    return(data)
}
# Plot the simulation of ZQdX
ZQdX_sim <- function(
    sim_event_mm_m3_s,
    Q_input = NULL,
    Qplot = TRUE,
    title_label,
    ylabel) {
    if (Qplot == TRUE && is.null(Q_input)) {
        stop("Q_input must be provided if Qplot is TRUE")
    }

    # Replace -9999 and -1e9 with NA
    sim_event_mm_m3_s[sim_event_mm_m3_s == -9999] <- NA
    sim_event_mm_m3_s[sim_event_mm_m3_s == -1e9] <- NA

    ZQplot <- ggplot(sim_event_mm_m3_s, aes(x = X3_obs)) +
        labs(title = title_label, x = "Streamwise position (meters)", y = ylabel) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))
    if (Qplot == TRUE) {
        ZQplot <- ZQplot +
            geom_point(aes(y = Y2_sim, col = "sim")) +
            scale_color_manual(values = c("sim" = "blue", "boundary\nconditions" = "red")) +
            facet_wrap(~X1_obs, ncol = 2)

        # Create a data frame for the horizontal lines
        hline_data <- data.frame(
            X1_obs = 1:length(Q_input),
            yintercept = Q_input,
            linetype = 1:length(Q_input), # or use named linetypes like c("solid", "dashed", "dotted")
            col = "boundary\nconditions" # same color for all lines
        )

        ZQplot <- ZQplot +
            geom_hline(
                data = hline_data,
                aes(yintercept = yintercept, linetype = factor(linetype), color = col),
                show.legend = TRUE
            )

        ZQplot <- ZQplot +
            scale_linetype_manual(values = c("1" = "dashed", "2" = "solid", "3" = "dotted", "4" = "twodash")) +
            labs(colour = NULL, linetype = "boundary\ncondition\nevents")
    } else {
        ZQplot <- ZQplot +
            geom_point(aes(y = Y1_sim, col = "sim")) +
            geom_point(aes(y = Y1_obs, col = "obs"))

        if (any(sim_event_mm_m3_s$Yu_z != 0)) {
            ZQplot <- ZQplot +
                geom_errorbar(aes(
                    ymin = Y1_obs - qnorm(0.975) * Yu_z,
                    ymax = Y1_obs + qnorm(0.975) * Yu_z, col = "obs"
                ))
        }

        ZQplot <- ZQplot +
            scale_color_manual(values = c("sim" = "blue", "obs" = "red")) +
            labs(colour = NULL) +
            facet_wrap(~X1_obs, ncol = 2)
    }
    return(ZQplot)
}

plot_DIC <- function(
    dir_polynomial,
    DIC_criterion = "DIC3") {
    all_files <- list.files(
        file.path(
            dir_polynomial
        ),
        recursive = TRUE
    )

    DIC_path_Files <- all_files[grep(all_files, pattern = "Results_DIC.txt")]
    if (length(DIC_path_Files) == 0) stop("Any Results_DIC.txt file was found")
    DIC_results <- c()


    for (i in 1:length(DIC_path_Files)) {
        DIC_by_degree <- read.table(file.path(dir_polynomial, DIC_path_Files[i]), col.names = c("Criteria", "Value"))
        # Match the criterion chosen
        DIC_by_degree <- DIC_by_degree[which(DIC_by_degree[, 1] == DIC_criterion), ]
        # Assign the polynomial degree
        extraction <- strsplit(DIC_path_Files[i], "/")[[1]][1]
        # # Extraire le chiffre après le underscore
        # chiffre <- sub(".*_(\\d+)", "\\1", extraction)

        # # Convertir en numérique (optionnel)
        DIC_by_degree$case <- extraction

        DIC_results <- rbind(DIC_results, DIC_by_degree)
    }
    min_local <- DIC_results[which.min(DIC_results$Value), ]

    DIC_plot <- ggplot(DIC_results, aes(x = factor(case), y = Value, col = factor(Criteria))) +
        geom_point(size = 3) +
        geom_point(data = min_local, aes(x = , factor(case), y = Value), col = "blue", size = 3) +
        annotate("segment",
            x = factor(min_local$Degree), y = min_local$Value * 1.005, xend = factor(min_local$Degree), yend = min_local$Value * 1.002,
            linewidth = 2, linejoin = "mitre",
            arrow = arrow(type = "closed", length = unit(0.01, "npc"))
        ) +
        theme_bw() +
        labs(
            x = "Calibration case",
            y = "Value",
            col = "Criterion",
            title = "DIC criterion :",
            subtitle = "Comparison of several polynomial degrees"
        ) +
        theme(
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5)
        )

    return(DIC_plot)
}

###############
# Theorical values
####################

Kmin_literature <- data.frame(
    x_start = c(0),
    x_end = c(18),
    min_value = c(1 / 0.013),
    max_value = c(1 / 0.009),
    mean_value = c(1 / 0.010)
)


Kflood_literature <- data.frame(
    x_start = c(0, 9.05),
    x_end = c(9.05, 18),
    min_value = c(mean(c(1 / 0.010, 1 / 0.012)), mean(c(1 / 0.025, 1 / 0.030))),
    mean_value = c(mean(c(1 / 0.015, 1 / 0.012)), mean(c(1 / 0.025, 1 / 0.035))),
    max_value = c(mean(c(1 / 0.014, 1 / 0.018)), mean(c(1 / 0.035, 1 / 0.05)))
)


segment_layer_reference <- function(K_literature) {
    K_segment_long <- data.frame()

    for (i in 1:nrow(K_literature)) {
        row <- K_literature[i, ]

        # Add min segment
        K_segment_long <- rbind(
            K_segment_long,
            data.frame(
                x_start = row$x_start,
                x_end = row$x_end,
                y_value = row$min_value,
                line_type = "Chow (1959)"
            )
        )

        # Add mean segment
        K_segment_long <- rbind(
            K_segment_long,
            data.frame(
                x_start = row$x_start,
                x_end = row$x_end,
                y_value = row$mean_value,
                line_type = "Chow (1959)"
            )
        )

        # Add max segment
        K_segment_long <- rbind(
            K_segment_long,
            data.frame(
                x_start = row$x_start,
                x_end = row$x_end,
                y_value = row$max_value,
                line_type = "Chow (1959)"
            )
        )
    }
    K_segment_layer <- geom_segment(
        data = K_segment_long,
        aes(x = x_start, y = y_value, xend = x_end, yend = y_value, linetype = line_type),
        color = "gray"
    )

    return(K_segment_layer)
}

Kmin_segment_layer <- segment_layer_reference(K_literature = Kmin_literature)
Kflood_segment_layer <- segment_layer_reference(K_literature = Kflood_literature)


#######################
plotDIC <- plot_DIC(path_experiment)
ggsave(
    file.path(
        path_experiment,
        "DIC.png"
    ),
    plotDIC,
    width = 20,
    height = 20,
    units = "cm"
)

# At each eid_cal_case
Q_observed <- c(162, 162)
final_calibraion <- FALSE
counter <- 1
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
    path_temp_plots <- file.path(path_polynomial, "BaM")
    path_post_traitement <- file.path(path_temp_plots, "post_traitement")
    path_post_traitement_data <- file.path(path_post_traitement, "RData")

    if (!dir.exists(path_post_traitement)) {
        dir.create(path_post_traitement)
    }

    if (!dir.exists(path_post_traitement_data)) {
        dir.create(path_post_traitement_data)
    }

    path_model_mage <- c(paste0(path_polynomial, "/model_mage/", all_events, "/"))

    mcmc_not_cooked <- readMCMC(file.path(path_temp_plots, "Results_MCMC.txt"))
    plots <- tracePlot(mcmc_not_cooked)

    mcmcplot <- wrap_plots(plots, ncol = 3)
    # Save plot and data
    ggsave(
        file.path(
            path_post_traitement,
            "MCMC_not_cooked.png"
        ),
        mcmcplot,
        width = 20,
        height = 20,
        units = "cm"
    )

    # Density plot for each parameter
    plots <- densityPlot(mcmc_not_cooked)
    pdf_plot <- wrap_plots(plots, ncol = 3)

    ggsave(
        file.path(
            path_post_traitement,
            "densityplot_not_cooked.png"
        ),
        pdf_plot,
        width = 20,
        height = 20,
        units = "cm"
    )
    if (final_calibraion) {
        if (!file.exists(file.path(path_temp_plots, "Results_Cooking.txt"))) stop("MCMC is still running or calculation is not going to the end. Verify if calibration is already finished or verify that calibration has not error messages. Please put final_results = FALSE as input")

        mcmc <- readMCMC(file.path(path_temp_plots, "Results_Cooking.txt"))
        plots <- tracePlot(mcmc)

        mcmcplot <- wrap_plots(plots, ncol = 3)


        ggsave(
            file.path(
                path_post_traitement,
                "MCMC_Cooked.png"
            ),
            mcmcplot,
            width = 20,
            height = 20,
            units = "cm"
        )

        # Density plot for each parameter
        plots <- densityPlot(mcmc)
        pdf_plot <- wrap_plots(plots, ncol = 3)

        ggsave(
            file.path(
                path_post_traitement,
                "densityplot_Cooked.png"
            ),
            pdf_plot,
            width = 20,
            height = 20,
            units = "cm"
        )


        png(
            file.path(path_post_traitement, "corelation_cooked.png"),
            width = 800,
            height = 800,
            res = 120
        )
        pairs(mcmc)
        dev.off()

        getSummary <- read.table(
            file = file.path(path_temp_plots, "Results_Summary.txt"),
            header = TRUE,
            stringsAsFactors = FALSE
        )

        # Values of error model in meter for WSE, in m3/s for discharge and m/s for velocity.
        knitr::kable(getSummary,
            align = "c"
        )


        # Zoom into the MAP and standard deviation of the error model
        getSummary_zoom <- getSummary[c(11, 16), ]
        getSummary_zoom[, c("Y1_intercept")] <- getSummary_zoom[, c("Y1_intercept")] * factor_scaled # WSE in mm

        # Get MAP simulation
        MAP_param_matrix <- as.numeric(getSummary_zoom[2, c(1:(length(Kmin_prior) + length(Kflood_prior)))])
    } else {
        results_MCMC_sampling <- read.table(file.path(path_temp_plots, "Results_MCMC.txt"), header = TRUE)
        mcmc <- data.frame(results_MCMC_sampling[, 1:(length(Kmin_prior) + length(Kflood_prior))])

        # Get MAP simulation: from current analysis without burning and slim
        MAP_param_matrix <- as.numeric(results_MCMC_sampling[which.max(results_MCMC_sampling$LogPost), 1:(length(Kmin_prior) + length(Kflood_prior))])
    }
    Kmin_calculations <- K_plot(
        matrix_spatialisation = Z_MatrixKmin,
        mcmc = mcmc,
        n_param_Kmin = length(Kmin_prior),
        n_param_Kflood = length(Kflood_prior),
        covariate_discretization = covariate_grid,
        MAP_param_vector = MAP_param_matrix,
        main_channel = TRUE
    )
    df_MAP_Kmin <- Kmin_calculations[[1]]
    Kmin_plot <- Kmin_calculations[[2]]
    df_envelope_Kmin <- Kmin_calculations[[3]]

    ls_spatial_friction_Kmin <- list(
        df_envelope = df_envelope_Kmin,
        df_MAP = df_MAP_Kmin
    )
    save(ls_spatial_friction_Kmin,
        file = file.path(path_post_traitement_data, "Data_friction_estimation_ls_spatial_friction_Kmin.RData")
    )

    final_Kmin_plot <- Kmin_plot + Kmin_segment_layer
    scale_linetype_manual(
        name = "Reference\nvalues",
        values = c("Chow (1959)" = "dashed")
    )

    ggsave(
        file.path(
            path_post_traitement,
            "Kmin.png"
        ),
        final_Kmin_plot,
        width = 20,
        height = 20,
        units = "cm"
    )
    save(final_Kmin_plot,
        file = file.path(path_post_traitement_data, "Plot_friction_estimation_plot_Kmin_plot.RData")
    )

    Kflood_calculations <- K_plot(
        matrix_spatialisation = Z_MatrixKflood,
        mcmc = mcmc,
        n_param_Kmin = length(Kmin_prior),
        n_param_Kflood = length(Kflood_prior),
        covariate_discretization = covariate_grid,
        MAP_param_vector = MAP_param_matrix,
        main_channel = FALSE
    )
    df_MAP_Kflood <- Kflood_calculations[[1]]
    Kflood_plot <- Kflood_calculations[[2]]
    df_envelope_Kflood <- Kflood_calculations[[3]]

    ls_spatial_friction_Kflood <- list(
        df_envelope = df_envelope_Kflood,
        df_MAP = df_MAP_Kflood
    )
    save(ls_spatial_friction_Kflood,
        file = file.path(path_post_traitement_data, "Data_friction_estimation_ls_spatial_friction_Kflood.RData")
    )

    final_Kflood_plot <- Kflood_plot + Kflood_segment_layer +
        scale_linetype_manual(
            name = "Reference\nvalues",
            values = c("Chow (1959)" = "dashed")
        )

    ggsave(
        file.path(
            path_post_traitement,
            "Kflood.png"
        ),
        final_Kflood_plot,
        width = 20,
        height = 20,
        units = "cm"
    )
    save(final_Kflood_plot,
        file = file.path(path_post_traitement_data, "Plot_friction_estimation_plot_Kflood_plot.RData")
    )
    ####################################
    # Residuals
    ###################################
    CalData <- read.table(file.path(path_temp_plots, "CalibrationData.txt"), header = TRUE)
    # Read residuals
    if (final_calibraion) {
        residuals <- read.table(
            file = file.path(path_temp_plots, "Results_Residuals.txt"),
            header = TRUE,
            stringsAsFactors = FALSE
        )

        # Convert residuals to mm and m3/s
        residuals_mm_m3_s <- data.frame(
            residuals[, 1:4],
            residuals[, c(9, 19, 20, 24, 29)] * factor_scaled
        )
        sim_event_mm_m3_s <- data.frame(residuals_mm_m3_s, Yu_z = CalData$Yu_z * factor_scaled)
    } else {
        # Residuals file is not ready, so I need to create by myself with MAP estimation from sampled data while MCMC is turning
        # Get data format from mage results
        temporal_dir <- tempdir()
        files <- list.files(temporal_dir, full.names = TRUE)

        exclude_file <- file.path(temporal_dir, "vscode-R")

        # List all files in the temp directory
        files <- list.files(temporal_dir, full.names = TRUE)

        # Exclude the specific file
        files_to_delete <- setdiff(files, exclude_file)

        # Remove remaining files if any
        if (length(files_to_delete) > 0) {
            unlink(files_to_delete, recursive = TRUE)
        }
        file.copy(
            from = path_model_mage,
            to = temporal_dir, recursive = TRUE
        )
        temp_path <- file.path(temporal_dir, basename(path_model_mage))

        path_RUGFile <- file.path(
            temp_path,
            list.files(temp_path)[grep(list.files(temp_path), pattern = ".RUG")]
        )
        MAP_RUGFile <- lapply(path_RUGFile, function(file) {
            read_fortran_data(
                file_path = file,
                col_widths_RUGFile = c(1, 3, 6, 10, 10, 10, 10),
                skip = 1
            )
        })
        df_MAP <- list(
            df_MAP_Kmin,
            df_MAP_Kflood
        )
        for (i in seq_along(MAP_RUGFile)) {
            if (nrow(df_MAP[[i]]) != nrow(MAP_RUGFile[[i]])) stop("Mage project has not the same size as the value estimated based on the MAP estimator")
            MAP_RUGFile[[i]]$V6 <- df_MAP[[1]]$Value

            write_RUGFile(
                RUG_path = path_RUGFile[i],
                RUGFile_data = MAP_RUGFile[[i]],
                RUG_format = "%1s%3d      %10.3f%10.3f%10.2f%10.2f"
            )
        }

        REPFile <- list.files(temp_path)[grep(list.files(temp_path), pattern = ".REP")]
        REPFile <- str_remove(REPFile, pattern = ".REP")

        for (i in 1:length(temp_path)) {
            setwd(temp_path[i])
            system2(MAGE_executable, args = REPFile[i], wait = TRUE)
        }

        runModel(
            workspace = temporal_dir,
            mod = mod_polynomials[[counter]],
            X = X,
            stout = NULL
        )
        counter <- counter + 1

        setwd(temporal_dir)

        sim <- read.table("Y.txt", header = TRUE)
        obs <- CalData[, 5:9]
        obs[obs == -9999] <- NA
        u_obs <- CalData[, 10:14]
        u_obs[u_obs == -9999] <- NA
        res <- obs - sim

        sim_event_mm_m3_s <- data.frame(
            X1_obs = X[, 1],
            X2_obs = X[, 2],
            X3_obs = X[, 3],
            X4_obs = X[, 4],
            Y1_obs = obs[, 1],
            Y1_sim = sim[, 1],
            Y2_sim = sim[, 2],
            Y1_res = res[, 1],
            Yu_z = u_obs[, 1]
        )
        setwd(dir_workspace)
    }

    Q_sim_vs_obs <- ZQdX_sim(
        sim_event_mm_m3_s = sim_event_mm_m3_s,
        Q_input = Q_observed,
        Qplot = TRUE,
        title_label = "MAP simulations vs boundary condition \n(Discharge upstream)",
        ylabel = "Q (L/s)"
    )
    ggsave(
        file.path(
            path_post_traitement,
            "Q_residuals.png"
        ),
        Q_sim_vs_obs,
        width = 20,
        height = 20,
        units = "cm"
    )
    save(Q_sim_vs_obs,
        file = file.path(path_post_traitement_data, "Plot_Q_sim_vs_obs.RData")
    )
    Z_sim_vs_obs <- ZQdX_sim(
        sim_event_mm_m3_s = sim_event_mm_m3_s,
        Q_input = NULL,
        Qplot = FALSE,
        title_label = "MAP simulations vs Observations (WSE)",
        ylabel = "Water surface elevation (mm)"
    )
    ggsave(
        file.path(
            path_post_traitement,
            "Plot_Z_residuals.png"
        ),
        Z_sim_vs_obs,
        width = 20,
        height = 20,
        units = "cm"
    )
    save(Z_sim_vs_obs,
        file = file.path(path_post_traitement_data, "Plot_Z_sim_vs_obs.RData")
    )
    save(sim_event_mm_m3_s,
        file = file.path(path_post_traitement_data, "Data_Z_sim_vs_obs.RData")
    )
}
