rm(list = ls())
graphics.off()

load("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/processed_data/HHLab/compound_transition_friction/data_HHLab_uniform_case.RData")



library(RBaM)
library(dplyr)
library(stringr)
# setPathToBaM("/home/famendezrios/Documents/Git/BaM/makefile/")

path_temp_results <- "/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/Case_studies/HHLab/Compound_channel_spatial_friction_test/Calibration"

mage_projet_name <- "CWMQ"

WSE_data_temp <- all_data_calibration$WSE

WSE_data_temp <- data.frame(WSE_data_temp %>%
    group_by(ID_experiment, x) %>%
    summarise(
        z_mean_all = mean(z_mean),
        Yu_mean_all = mean(Yu), .groups = "drop"
    ))

case <- "CWMQ"
WSE_data <- WSE_data_temp[str_detect(WSE_data_temp$ID_experiment, pattern = case), ]


Cal_measures <- list(
    CWMQ18 = c(2.225, 4.310, 5.270, 6.240, 8.315, 9.115, 9.250, 10.250, 11.250, 12.250, 13.250),
    CWMQ12 = c(2.225, 2.975, 3.345, 4.095, 5.060, 6.240, 6.990, 8.030, 8.795, 9.065, 9.115, 9.545, 10.030, 11.370, 12.120, 12.325, 13.290)
)

CalData_sample <- c()
for (i in 1:length(Cal_measures)) {
    if (i == 1) {
        CalData_sample <- WSE_data[WSE_data$ID_experiment %in% names(Cal_measures)[i] &
            WSE_data$x %in% Cal_measures[[i]], ]
    } else {
        CalData_sample <- rbind(
            CalData_sample,
            WSE_data[WSE_data$ID_experiment %in% names(Cal_measures)[i] &
                WSE_data$x %in% Cal_measures[[i]], ]
        )
    }
}


X <- data.frame(
    event = c(rep(1, lengths(Cal_measures)[1]), rep(2, lengths(Cal_measures)[2])),
    reach = rep(1, sum(lengths(Cal_measures))),
    x = CalData_sample$x,
    t = rep(3420, nrow(CalData_sample))
)

# position_direction
# increasing: positions are increasing in flow direction
# decreasing: positions are decreasing in flow direction
position_direction <- "increasing"

# Give the reaches connected
# Reaches for spatially distributed the friction : the order here will be important (be careful depending on the mage model set up). That will impact the way that positions are understanded.

case <- "independant"
# case <- "follow" independant monoreach

if (case == "independant") {
    reaches_user <- list(
        1,
        # c(1, 2)
        2
    )
} else if (case == "follow") {
    reaches_user <- list(
        # 1,
        c(1, 2)
    )
} else if (case == "monoreach") {
    reaches_user <- list(
        1
    )
} else {
    stop("ERROR")
}


# Number of different reaches
Nb_reaches_estimation <- length(reaches_user)

# For example: position_direction = increasing and reaches_user is c(1,2)

# Only a unit will be considered, composed of two biefs, one and two. A check will be performed to check if all positions are increasing as defined in input data

# Second example: position_direction = decreasing and reaches_user is c(2,1)

# In this case, reach 2 is downstream and reach 1 is upstream. A check will be performed to be sure that position will decrease by reach

n_degree_Kmin <- 0
n_degree_Kflood <- 0

# Input parameters in main channel:
# Reach 1 : main channel (wood)
a0_min_reach_1 <- RBaM::parameter(
    name = "a0_min_glass",
    init = 1 / 0.010, # Initial guess
    prior.dist = "FlatPrior+", # Prior distribution
    prior.par = NULL
) # Parameters of the prior distribution


# reach 1 : floodplain
a0_flood_reach_1 <- RBaM::parameter(
    name = "a0_flood_wood",
    init = 33,
    prior.dist = "FlatPrior+",
    prior.par = NULL
)

# reach 2 : floodplain
a0_flood_reach_2 <- RBaM::parameter(
    name = "a0_flood_grass",
    init = 1 / 0.013,
    prior.dist = "FlatPrior+",
    prior.par = NULL
)

if (n_degree_Kmin == 0 && n_degree_Kflood == 0) {
    Kmin_prior <- list(
        a0_min_reach_1
    )
    Kflood_prior <- list(
        a0_flood_reach_1,
        a0_flood_reach_2
    )
} else {
    stop("case not supported")
}


# if (length(Kmin_prior) != (Nb_reaches_estimation * (n_degree_Kmin + 1))) stop(paste0("More prior information (", length(Kmin_prior), ") than the number of reaches for estimation (", Nb_reaches_estimation, ")"))

theta_param <- c(Kmin_prior, Kflood_prior)

mageDir <- c(
    "/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/Case_studies/HHLab/Compound_channel_spatial_friction_test/CWMQ18/",
    "/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/Case_studies/HHLab/Compound_channel_spatial_friction_test/CWMQ12/"
)
### interpolation

specific_points_Model <- data.frame(
    reach = c(1, 2),
    KP_start = c(0, 9.05),
    KP_end = c(9.05, 18)
)



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


# Real reaches from Mage projet (interpolation by reach)
grid_covariant_discretized_by_reach <- list()
for (i in 1:nrow(specific_points_Model)) {
    boundaries_points <- c(specific_points_Model$KP_start[i], specific_points_Model$KP_end[i])

    grid_covariant_discretized_by_reach[[i]] <- data.frame(
        real_reach = specific_points_Model$reach[i],
        interpol_values = interpolation_specific_points(
            total_points = 100,
            all_specific_points = boundaries_points
        )
    )
}
grid_covariant_discretized_real_bief <- bind_rows(grid_covariant_discretized_by_reach)

grid_covariant_discretized <- c()
save_original_data <- list()
# Change to reach for spatialisation
for (i in 1:length(reaches_user)) {
    mask <- which(grid_covariant_discretized_real_bief$real_reach %in% reaches_user[[i]])

    if (length(mask) == 0) stop(paste0("Real reaches are ", specific_points_Model$reach, " and the reaches defined by the user are not matching: ", reaches_user))

    df_covariant_by_reach <- grid_covariant_discretized_real_bief[mask, ]

    save_original_data[[i]] <- df_covariant_by_reach

    if (length(reaches_user[[i]]) != 1) {
        # Case union of real bief for spatialisation
        if (position_direction == "increasing") {
            # Check if all interpol_values are positive for each real_reach
            check <- df_covariant_by_reach %>%
                group_by(real_reach) %>%
                summarize(check = all(diff(interpol_values) > 0, na.rm = TRUE))


            if (!any(check$check)) stop("All position values should increase in flow direction")
        } else if (position_direction == "decreasing") {
            # Check if all interpol_values are positive for each real_reach
            check <- df_covariant_by_reach %>%
                group_by(real_reach) %>%
                summarize(check = all(diff(interpol_values) < 0, na.rm = TRUE))

            if (!all(check$check)) stop(paste0("All position values should decrease in flow direction. Please check reach(es) ", check$real_reach[check$check == FALSE]))
        } else {
            stop(paste0("The position_direction variable should be either increasing or decreasing in flow direction all the time"))
        }
    }
    #### Save data independently of the case
    if (i == 1) {
        grid_covariant_discretized <- data.frame(
            Reach = i,
            Covariate = df_covariant_by_reach$interpol_values
        )
    } else {
        grid_covariant_discretized_temp <- data.frame(
            Reach = i,
            Covariate = df_covariant_by_reach$interpol_values
        )
        grid_covariant_discretized <- rbind(
            grid_covariant_discretized, grid_covariant_discretized_temp
        )
    }
}

min_max_values <- aggregate(grid_covariant_discretized[, 2] ~ grid_covariant_discretized[, 1],
    data = grid_covariant_discretized,
    FUN = function(x) c(min = min(x), max = max(x))
)

min_max_values <- do.call(data.frame, min_max_values)
colnames(min_max_values) <- c("reach", "min_value", "max_value")

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

# Write Legendre vector:
write.table(
    grid_covariant_discretized,
    file = file.path(path_temp_results, "vector_legendre.txt"), row.names = F
)



# Spatial distributed friction in the main channel

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

grid_covariant_discretized_kmin <- grid_covariant_discretized
grid_covariant_discretized_kmin$Reach <- 1
# Spatial distributed friction in the main channel
matrix_zFileKmin <- KFile_spatial(
    reaches = reaches_user[[1]],
    max_degree = n_degree_Kmin,
    grid_covariant_discretized = grid_covariant_discretized_kmin
)

# Spatial distributed friction in the floodplain
matrix_zFileKmoy <- KFile_spatial(
    reaches = reaches_user,
    max_degree = n_degree_Kflood, grid_covariant_discretized = grid_covariant_discretized
)

zFileKmin <- file.path(path_temp_results, "Zfile_Kmin.txt")
zFileKmoy <- file.path(path_temp_results, "Zfile_Kflood.txt")


xtra <- xtraModelInfo(
    fname = "Config_setup.txt",
    object = list(
        exeFile = "/home/famendezrios/Documents/Softwares/pamhyr2/mage8/mage",
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

mod <- model(
    ID = "MAGE_ZQV",
    nX = 4,
    nY = 5, ,
    par = theta_param,
    xtra = xtra
)

# Re-write RUG file : keep the position of the cross sections as original

col_widths_RUGFile <- c(1, 3, 6, 10, 10, 10, 10)
read_fortran_data <- function(file_path, col_widths_RUGFile, skip = 0) {
    # Read the file with the fixed-width format
    data <- read.fwf(file_path, widths = col_widths_RUGFile, header = FALSE, skip = skip)
    return(data)
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

prior_theta_param <- get_init_prior(theta_param)
prior_all_dist_theta_param <- get_init_prior(theta_param, FIX_dist = TRUE)

RUGFiles <- paste0(mageDir, mage_projet_name, ".RUG")

# Harcode values, it does not matter, it just useful for give the size of the .RUG File
RUG_Kmin <- 20
RUG_Kmoy <- 10

# j is a index for multi-events
for (j in seq_along(RUGFiles)) {
    # Calculate original covariant considering that grid_covariant_discretized is in the middle of the real covariant value except on the boundaries
    grid_covariant_meshed_Model <- RUG_KP_end_by_reach <- RUG_KP_start_by_reach <- RUG_id_reach_by_reach <- RUG_Kmin_by_reach <- RUG_Kmoy_by_reach <- list()

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
        RUG_id_reach_by_reach[[reach]] <- rep(1, length(covariant_discretized_by_reach))
        # RUG_id_reach_by_reach[[reach]] <- rep(specific_points_Model$reach[reach], length(covariant_discretized_by_reach))

        RUG_KP_start_by_reach[[reach]] <- grid_covariant_meshed_Model[[reach]]$real_values[-length(grid_covariant_meshed_Model[[reach]]$real_values)]

        RUG_KP_end_by_reach[[reach]] <- grid_covariant_meshed_Model[[reach]]$real_values[-1]

        RUG_Kmin_by_reach[[reach]] <- rep(RUG_Kmin[reach], length(covariant_discretized_by_reach))
        RUG_Kmoy_by_reach[[reach]] <- rep(RUG_Kmoy[reach], length(covariant_discretized_by_reach))
    }

    RUG_KP_start <- unlist(RUG_KP_start_by_reach)
    RUG_KP_end <- unlist(RUG_KP_end_by_reach)
    RUG_id_reach <- unlist(RUG_id_reach_by_reach)
    RUG_Kmin_args <- unlist(RUG_Kmin_by_reach)
    RUG_Kmoy_args <- unlist(RUG_Kmoy_by_reach)


    write_RUGFile(
        RUG_path = RUGFiles[j],
        RUG_id_reach = RUG_id_reach,
        RUG_KP_start = RUG_KP_start,
        RUG_KP_end = RUG_KP_end,
        RUG_Kmin = RUG_Kmin,
        RUG_Kmoy = RUG_Kmoy,
        RUG_format = "%1s%3d      %10.0f%10.0f%10.2f%10.2f"
    )
}

Y <- data.frame(WSE = CalData_sample$z_mean_all)
Y$Discharge <- -9999
Y$Velocity <- -9999
Y$Y_Kmin <- -9999
Y$Y_Kmoy <- -9999
Yu <- data.frame(Yu_z = rep(0, length(CalData_sample$Yu_mean_all)))
# Yu <- data.frame(Yu_z = CalData_sample$Yu_mean_all)
Yu$Yu_Q <- -9999
Yu$Yu_V <- -9999
Yu$Yu_Kmin <- -9999
Yu$Yu_Kmoy <- -9999

# Y$Y_Kmin[4] <- 100
# Yu$Yu_Kmin[4] <- 0.1

# Y$Y_Kmoy[which(X$x %in% c(6.5))] <- 33
# Yu$Yu_Kmoy[which(X$x %in% c(6.5))] <- 0.5
# Y$Y_Kmoy[which(X$x %in% c(3.5, 5, 6.5, 7, 7.5, 9.25, 11.75, 12.5, 13))] <- 33
# Yu$Yu_Kmoy[which(X$x %in% c(3.5, 5, 6.5, 7, 7.5, 9.25, 11.75, 12.5, 13))] <- 5

data <- dataset(X = X, Y = Y, Yu = Yu, data.dir = file.path(path_temp_results))

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
            init = log(10),
            prior.dist = "LogNormal",
            prior.par = c(log(10), 1)
        ))
    ),
    remnantErrorModel(
        fname = "Config_RemnantSigma4.txt",
        funk = "Constant",
        par = list(parameter(
            name = "intercept",
            init = log(10),
            prior.dist = "LogNormal",
            prior.par = c(log(10), 0.2)
        ))
    )
)

prior_error_model <- get_init_prior(remant_error_list)

jump_MCMC_theta_param_user <- c(8, 8, 5)
jump_MCMC_error_model_user <- 0.003
threshold_jump_MCMC_error_model <- 0.5

jump_MCMC_theta_param <- ifelse(prior_theta_param != 0,
    prior_theta_param[(prior_theta_param != 0)] * 0.1,
    jump_MCMC_theta_param_user
)

# jump_MCMC_theta_param <- c(15, 10)

jump_MCMC_error_model <- list()
for (ind in 1:length(prior_error_model)) {
    jump_MCMC_error_model[[ind]] <- ifelse(prior_error_model[[ind]] > threshold_jump_MCMC_error_model,
        prior_error_model[[ind]][(prior_error_model[[ind]] != 0)] * 0.1,
        jump_MCMC_error_model_user
    )
}

mcmcOptions_user <- mcmcOptions(
    nCycles = 25,
    nAdapt = 25,
    manualMode = TRUE,
    thetaStd = jump_MCMC_theta_param,
    gammaStd = jump_MCMC_error_model
)
mcmcCooking_user <- mcmcCooking(
    burn = 0,
    nSlim = 1
)
mcmcSummary_user <- mcmcSummary(xtendedMCMC.fname = "Results_xtendedMCMC.txt")

# runModel(workspace = tempdir(), mod = mod, X = X)

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
    run = TRUE,
    preClean = FALSE,
    workspace = path_temp_results,
    # dir.exe = file.path(find.package("RBaM"), "bin"),
    # name.exe = "BaM",
    predMaster_fname = "Config_Pred_Master.txt"
)


library(patchwork)

mcmc <- readMCMC(file.path(path_temp_results, "Results_MCMC.txt"))
plots <- tracePlot(mcmc)

mcmcplot <- wrap_plots(plots, ncol = 3)

library(ggplot2)

ggsave(
    file.path(
        path_temp_results,
        "MCMC.png"
    ),
    mcmcplot
)


png(
    file.path(path_temp_results, "corelation.png"),
    width = 800,
    height = 800,
    res = 120
)
pairs(mcmc[, c("a0_min_glass", "a0_flood_wood", "a0_flood_gras", "Y1_intercept", "Y4_intercept")])
dev.off()


residuals <- read.table(file.path(path_temp_results, "Results_Residuals.txt"), header = TRUE)

ggplot(residuals, aes(x = X3_obs)) +
    geom_point(aes(y = Y1_obs, col = "z_obs")) +
    geom_point(aes(y = Y1_sim, col = "z_sim")) +
    facet_wrap(~X1_obs, scales = "free_y")


getSummary <- read.table(
    file = file.path(path_temp_results, "Results_Summary.txt"),
    header = TRUE,
    stringsAsFactors = FALSE
)

# Values of error model in meter for WSE, in m3/s for discharge and m/s for velocity.
knitr::kable(getSummary,
    align = "c"
)

library(dplyr)
library(tidyr)
# Zoom into the MAP and standard deviation of the error model
getSummary_zoom <- getSummary[c(11, 16), ]
getSummary_zoom[, c(5, 6)] <- getSummary_zoom[, c(5, 6)] * 1000 # Convert to mm for WSE and discharge
getSummary_zoom[, 7] <- getSummary_zoom[, 7] * 100 # Convert to cm/s for velocity

n_degree <- 0
n_degree_flood <- 1

ks_literature <- data.frame(min = 1 / 0.013, max = 1 / 0.009, mean = 1 / 0.010)

# Get MAP simulation
MAP_param_matrix <- as.numeric(getSummary_zoom[2, c(1:(n_degree + 1 + n_degree_flood + 1))])

matrix_zFileKmin <- read.table(file.path(path_temp_results, "Zfile_Kmin.txt"), header = TRUE)

position <- read.table(file.path(path_temp_results, "vector_legendre.txt"), header = TRUE)

k_estimated_all <- as.data.frame(as.matrix(matrix_zFileKmin) %*% as.matrix(t(mcmc[, 1:(n_degree + 1)])))

k_estimated_all$KP <- position[, 2]

# Convert to long format
df_MCMC_sampling <- pivot_longer(
    k_estimated_all,
    cols = -KP,
    values_to = "Value"
) %>%
    select(KP, Value) %>%
    mutate(ID = "MCMC Sampling")

k_estimated_MAP <- as.matrix(matrix_zFileKmin) %*% MAP_param_matrix[1:(n_degree + 1)]

df_MAP <- data.frame(
    KP = position[, 2],
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

kmin_plot <- ggplot() +
    geom_ribbon(
        data = df_envelope,
        aes(x = KP, ymin = ymin, ymax = ymax, fill = ID)
    ) +
    geom_line(
        data = df_MAP,
        aes(x = KP, y = Value, color = ID)
    ) +
    geom_hline(aes(yintercept = ks_literature$min, linetype = "ASCI (1980)"), color = "gray") +
    geom_hline(aes(yintercept = ks_literature$max, linetype = "ASCI (1980)"), color = "gray") +
    geom_hline(aes(yintercept = ks_literature$mean, linetype = "ASCI (1980)"), color = "gray") +
    scale_fill_manual(values = c("Parametric\nuncertainty" = "pink")) +
    scale_color_manual(values = c("MAP" = "black")) +
    # Échelle de type de ligne
    scale_linetype_manual(
        name = "Reference\nvalues",
        values = c("ASCI (1980)" = "dashed")
    ) +
    labs(
        title = "Friction coefficient estimation in the main channel \nwith parametric uncertainty",
        x = "Lengthwise position (meters)",
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

ggsave(
    file.path(
        path_temp_results,
        "kmin.png"
    ),
    kmin_plot
)


ks_literature_floodplain <- data.frame(min = 1 / 0.05, max = 1 / 0.025, mean = 1 / 0.033)

matrix_zFileKflood <- read.table(file.path(path_temp_results, "Zfile_Kflood.txt"), header = TRUE)

k_estimated_all_floodplain <- as.data.frame(as.matrix(matrix_zFileKflood) %*% as.matrix(t(mcmc[, (n_degree + 1 + 1):(n_degree + 1 + n_degree_flood + 1)])))

k_estimated_all_floodplain$KP <- position[, 2]

# Convert to long format
df_MCMC_sampling <- pivot_longer(
    k_estimated_all_floodplain,
    cols = -KP,
    values_to = "Value"
) %>%
    select(KP, Value) %>%
    mutate(ID = "MCMC Sampling")

k_estimated_MAP <- as.matrix(matrix_zFileKflood) %*% MAP_param_matrix[(n_degree + 1 + 1):(n_degree + 1 + 1 + n_degree_flood)]

df_MAP <- data.frame(
    KP = position[, 2],
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

kmoy_plot <- ggplot() +
    geom_ribbon(
        data = df_envelope,
        aes(x = KP, ymin = ymin, ymax = ymax, fill = ID)
    ) +
    geom_line(
        data = df_MAP,
        aes(x = KP, y = Value, color = ID)
    ) +
    geom_hline(aes(yintercept = ks_literature_floodplain$min, linetype = "Chow (1959)"), color = "gray") +
    geom_hline(aes(yintercept = ks_literature_floodplain$max, linetype = "Chow (1959)"), color = "gray") +
    geom_hline(aes(yintercept = ks_literature_floodplain$mean, linetype = "Chow (1959)"), color = "gray") +
    scale_fill_manual(values = c("Parametric\nuncertainty" = "pink")) +
    scale_color_manual(values = c("MAP" = "black")) +
    # Échelle de type de ligne
    scale_linetype_manual(
        name = "Reference\nvalues",
        values = c("Chow (1959)" = "dashed")
    ) +
    labs(
        title = "Friction coefficient estimation in the floodplain \nwith parametric uncertainty",
        x = "Lengthwise position (meters)",
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

ggsave(
    file.path(
        path_temp_results,
        "kmoy.png"
    ),
    kmoy_plot
)

# Replace -9999 and -1e9 with NA
residuals[residuals == -9999] <- NA
residuals[residuals == -1e9] <- NA

# Convert residuals to mm and m3/s
residuals_mm_m3_s <- data.frame(
    residuals[, 1:4],
    residuals[, c(9, 19, 20, 24, 29)] * 1000
)

residuals_event_mm_m3_s <- data.frame(residuals_mm_m3_s, Yu_z = Yu$Yu_z * 1000)

# Plot WSE residuals
WSE_plot <- ggplot(residuals_event_mm_m3_s, aes(x = X3_obs)) +
    geom_point(aes(y = Y1_sim, col = "sim")) +
    labs(title = "MAP simulations vs Observations (WSE)", x = "Position (m)", y = "Water surface elevation (mm)") +
    geom_point(aes(y = Y1_obs, col = "obs")) +
    geom_errorbar(aes(
        ymin = Y1_obs - 1.96 * Yu_z,
        ymax = Y1_obs + 1.96 * Yu_z, col = "obs"
    )) +
    scale_color_manual(values = c("sim" = "blue", "obs" = "red")) +
    labs(colour = NULL, shape = "Events") +
    theme_bw() +
    facet_wrap(~X1_obs, ncol = 2) +
    theme(plot.title = element_text(hjust = 0.5), )

# Plot Q residuals
Q_plot <- ggplot(residuals_event_mm_m3_s, aes(x = X3_obs)) +
    geom_point(aes(y = Y2_sim, col = "sim")) +
    labs(title = "MAP simulations vs boundary condition in Q", x = "Position (m)", y = "Q (L/s)") +
    geom_hline(aes(yintercept = 162, col = "boundary\nconditions", linetype = "1")) +
    scale_color_manual(values = c("sim" = "blue", "boundary\nconditions" = "red")) +
    scale_linetype_manual(values = c("1" = "dashed", "2" = "solid", "3" = "dotted", "4" = "twodash")) +
    labs(colour = NULL, shape = "Simulated\nevents", linetype = "boundary\ncondition\nevents") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))


ggsave(
    file.path(
        path_temp_results,
        "residual_WSE.png"
    ),
    WSE_plot
)

ggsave(
    file.path(
        path_temp_results,
        "residual_Q.png"
    ),
    Q_plot
)
