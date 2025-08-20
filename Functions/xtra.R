interpolation_specific_points <- function(total_points = 100,
                                          all_specific_points) {
    if (total_points <= length(all_specific_points)) stop("Total points defined is lower or equal to the specific points required. Please either increase the number of total_points or not used interpolation function")
    # Calculate the intervals between the specific points
    intervals <- diff(all_specific_points)
    # Distribute the number of points proportionally to each interval
    total_local_points <- total_points - 1 # last value during diff function
    points_per_interval_float <- (intervals / sum(intervals)) * total_local_points

    points_per_interval <- round(points_per_interval_float)
    # Due to round function, we need to check the size
    # Check if more points than expected, need to remove some of them
    if (sum(points_per_interval) > total_local_points) {
        id_order <- order(points_per_interval_float %% 1)
        i <- 1
        while (sum(points_per_interval) > total_local_points) {
            id_position <- which(i == id_order)
            points_per_interval[id_position] <- points_per_interval[id_position] - 1
            i <- i + 1
        }
        # Check if less points than expected, need to add some of them
    } else if (sum(points_per_interval) < total_local_points) {
        id_order <- order(points_per_interval_float %% 1)
        i <- 1
        while (sum(points_per_interval) < total_local_points) {
            id_position <- which(i == id_order)
            points_per_interval[id_position] <- points_per_interval[id_position] + 1
            i <- i + 1
        }
    }

    # Generate sub-sequences for each interval, excluding the duplicate endpoint
    specific_seq_temp <- unlist(mapply(
        function(start, end, n) seq(start, end, length.out = n + 1)[-(n + 1)], # Exclude the endpoint
        all_specific_points[-length(all_specific_points)], # Start of each interval
        all_specific_points[-1], # End of each interval
        points_per_interval # Number of points per interval
    ))

    #####  reviser le fait que j'ai à la fin 101 valeurs dans la séquance, alors que je demande 100 points, mais les boundary font le bordel
    # Add the final endpoint manually
    specific_seq <- c(specific_seq_temp, all_specific_points[length(all_specific_points)])

    if (!all(all_specific_points %in% specific_seq)) stop("Distance so close, increase number of interpolate KP to get all points required")

    if (length(specific_seq) != total_points) stop("Something is wrong with the interpolation")
    return(specific_seq)
}



getlegendre <- function(
    max_degree,
    covariate_discretization) {
    normalized_values <- 2 * (covariate_discretization - min(covariate_discretization)) / (max(covariate_discretization) - min(covariate_discretization)) - 1

    # Write co variant matrix for main channel
    # Create an empty data frame and fill it with polynomial values
    legendre_df_covariate <- data.frame(x = normalized_values) # Start with x values
    legendre_df <- legendre_df_covariate

    if (max_degree < 0) stop("The maximum degree must be positive")
    if (any(!dplyr::between(normalized_values, -1, 1))) stop("Range of normalized_values should be between [-1,1]")

    if (max_degree == 0) {
        return(rep(1, length(normalized_values)))
    }
    if (max_degree == 1) {
        return(normalized_values)
    }

    # P_0(normalized_values) = 1
    Pn_1 <- rep(1, length(normalized_values))
    # P_1(normalized_values) = normalized_values
    Pn <- normalized_values

    for (k in 1:(max_degree - 1)) {
        P_next <- ((2 * k + 1) * normalized_values * Pn - k * Pn_1) / (k + 1)
        Pn_1 <- Pn
        Pn <- P_next
    }
    return(Pn)
}


KFile_spatial <- function(
    reaches,
    max_degree,
    grid_covariant_discretized) {
    # Create matrix for spatialisation
    matrix_spatialisation <- as.data.frame(matrix(
        data = 0, nrow = nrow(grid_covariant_discretized),
        ncol = (max_degree + 1) * (length(reaches))
    ))
    # Create the sequences
    first_seq <- paste0("P", seq(0, max_degree))
    second_seq <- seq(length(reaches))

    # Use outer to combine all elements of the sequences
    names_columns <- outer(first_seq, second_seq, FUN = paste, sep = "_Reach_")
    colnames(matrix_spatialisation) <- as.vector(names_columns)

    counter_position <- 1
    # Loop for the reach
    for (reach in 1:length(reaches)) {
        ZFile_id_row <- which(grid_covariant_discretized$Reach == reach)

        cov_values_by_reach <- grid_covariant_discretized$Covariate[ZFile_id_row]

        # Loop for the polynomial degree
        for (degree in 0:max_degree) {
            pn_values <- getlegendre(
                max_degree = degree,
                covariate_discretization = cov_values_by_reach
            )
            matrix_spatialisation[ZFile_id_row, counter_position] <- pn_values

            counter_position <- counter_position + 1
        }
    }

    return(matrix_spatialisation)
}



get_covariant_meshed <- function(covariant_discretized, specific_points_Model) {
    # Reconstruction of covariant_meshed_temps :
    covariant_meshed_temps <- numeric()

    # 1. Find the limit points for the first two data, then for the rest of the data
    covariant_meshed_temps[1] <- sum(covariant_discretized[c(1, 2)]) / 2

    for (i in 2:(length(covariant_discretized) - 1)) {
        # Loop starts at 2 because first data is already given
        # Loop ends at n-2 number of covariant_discretized because
        # to find limits points, at least two points are required and the limit values will be added manually
        covariant_meshed_temps[i] <- sum(covariant_discretized[c(i + 1, i)]) / 2
    }
    # Add model's specific points to the mesh. Necessary condition for running Mage model
    covariant_meshed <- unique(sort(c(
        specific_points_Model$KP_start[1],
        specific_points_Model$KP_end,
        covariant_meshed_temps
    )))
    return(covariant_meshed)
}

# Function to write RUG file
write_RUGFile <- function(RUG_path,
                          RUG_id_reach,
                          RUG_KP_start,
                          RUG_KP_end,
                          RUG_Kmin,
                          RUG_Kmoy,
                          RUG_format) {
    # Open a .RUG file for writing
    fileConn <- file(RUG_path, "w")
    # Write the first line as a comment
    writeLines("* This file is generated by PAMHYR, please don't modify", fileConn)
    # Loop to write each line in Fortran format
    for (i in 1:length(RUG_KP_start)) {
        # Format the line according to Fortran style
        line <- sprintf(
            "%1s%3d      %10.3f%10.3f%10.2f%10.2f",
            "K",
            RUG_id_reach[i],
            RUG_KP_start[i],
            RUG_KP_end[i],
            RUG_Kmin,
            RUG_Kmoy
        )
        # Write to file
        writeLines(line, fileConn)
    }
    # Close the file
    close(fileConn)
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

read_fortran_data <- function(file_path, col_widths_RUGFile, skip = 0) {
    # Read the file with the fixed-width format
    data <- read.fwf(file_path, widths = col_widths_RUGFile, header = FALSE, skip = skip)
    return(data)
}
