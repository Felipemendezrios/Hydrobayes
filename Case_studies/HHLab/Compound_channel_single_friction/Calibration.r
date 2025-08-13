rm(list = ls())
graphics.off()

load("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/processed_data/HHLab/compound_single_friction/data_HHLab_uniform_case.RData")



library(RBaM)
library(dplyr)

# setPathToBaM("/home/famendezrios/Documents/Git/BaM/makefile/")

path_temp_results <- "/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/Case_studies/HHLab/Compound_channel_single_friction/Calibration"

mage_projet_name <- "compound_channel"

WSE_data_temp <- all_data_calibration$WSE

WSE_data <- data.frame(WSE_data_temp %>%
    group_by(x) %>%
    summarise(
        z_mean_all = mean(z_mean),
        Yu_mean_all = mean(Yu)
    ))

Cal_measures <- c(3.5, 5, 6.5, 7, 7.5, 9.25, 11.75, 12.5, 13)

CalData_sample <- WSE_data[WSE_data$x %in% Cal_measures, ]

X <- data.frame(
    event = rep(1, nrow(CalData_sample)),
    reach = rep(1, nrow(CalData_sample)),
    x = CalData_sample$x,
    t = rep(3420, nrow(CalData_sample))
)

# position_direction
# increasing: positions are increasing in flow direction
# decreasing: positions are decreasing in flow direction
position_direction <- "increasing"

# Give the reaches connected
# Reaches for spatially distributed the friction : the order here will be important (be careful depending on the mage model set up). That will impact the way that positions are interpretated.

case <- "monoreach"
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
n_degree_Kflood <- 1

# Input parameters in main channel:
# Reach 1 : main channel
a0_min_reach_1 <- RBaM::parameter(
    name = "a0_min_r_1",
    init = 1 / 0.010, # Initial guess
    prior.dist = "FlatPrior+", # Prior distribution
    prior.par = NULL
) # Parameters of the prior distribution

a1_min_reach_1 <- RBaM::parameter(
    name = "a1_min_r_1",
    init = 0,
    prior.dist = "FlatPrior",
    prior.par = NULL
)

# reach 1 : floodplain
a0_flood_reach_1 <- RBaM::parameter(
    name = "a0_flood_r_1",
    init = 33,
    prior.dist = "FlatPrior+",
    prior.par = NULL
)

a1_flood_reach_1 <- RBaM::parameter(
    name = "a1_flood_r_1",
    init = 0,
    prior.dist = "FlatPrior",
    prior.par = NULL
)


### Give prior about the reach for estimation.
# First all information about a reach and then the next one
if (n_degree_Kmin == 0 && n_degree_Kflood == 0) {
    Kmin_prior <- list(
        a0_min_reach_1
        # a1_min_reach_1
    )
    Kflood_prior <- list(
        a0_flood_reach_1
        # a1_flood_reach_1
    )
} else if (n_degree_Kmin == 0 && n_degree_Kflood == 1) {
    Kmin_prior <- list(
        a0_min_reach_1
        # a1_min_reach_1
    )
    Kflood_prior <- list(
        a0_flood_reach_1,
        a1_flood_reach_1
    )
} else {
    stop("case not supported")
}


if (length(Kmin_prior) != (Nb_reaches_estimation * (n_degree_Kmin + 1))) stop(paste0("More prior information (", length(Kmin_prior), ") than the number of reaches for estimation (", Nb_reaches_estimation, ")"))

theta_param <- c(Kmin_prior, Kflood_prior)

mageDir <- c("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/Case_studies/HHLab/Compound_channel_single_friction/model_mage/")
### interpolation

specific_points_Model <- data.frame(
    reach = c(1),
    KP_start = c(0),
    KP_end = c(18)
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


# Spatial distributed friction in the main channel
matrix_zFileKmin <- KFile_spatial(
    reaches = reaches_user,
    max_degree = n_degree_Kmin,
    grid_covariant_discretized = grid_covariant_discretized
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
        RUG_id_reach_by_reach[[reach]] <- rep(specific_points_Model$reach[reach], length(covariant_discretized_by_reach))

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
Yu <- data.frame(Yu_z = CalData_sample$Yu_mean_all)
Yu$Yu_Q <- -9999
Yu$Yu_V <- -9999
Yu$Yu_Kmin <- -9999
Yu$Yu_Kmoy <- -9999

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
            init = 15,
            prior.dist = "LogNormal",
            prior.par = c(log(15), 0.2)
        ))
    ),
    remnantErrorModel(
        fname = "Config_RemnantSigma3.txt",
        funk = "Constant",
        par = list(parameter(
            name = "intercept",
            init = 15,
            prior.dist = "LogNormal",
            prior.par = c(log(15), 0.2)
        ))
    )
)

prior_error_model <- get_init_prior(remant_error_list)

jump_MCMC_theta_param_user <- 8
jump_MCMC_error_model_user <- 0.003
threshold_jump_MCMC_error_model <- 0.5

jump_MCMC_theta_param <- ifelse(prior_theta_param != 0,
    prior_theta_param[(prior_theta_param != 0)] * 0.1,
    jump_MCMC_theta_param_user
)

jump_MCMC_theta_param <- c(15, 10, 0.5)

jump_MCMC_error_model <- list()
for (ind in 1:length(prior_error_model)) {
    jump_MCMC_error_model[[ind]] <- ifelse(prior_error_model[[ind]] > threshold_jump_MCMC_error_model,
        prior_error_model[[ind]][(prior_error_model[[ind]] != 0)] * 0.1,
        jump_MCMC_error_model_user
    )
}

mcmcOptions_user <- mcmcOptions(
    nCycles = 75,
    nAdapt = 75,
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

wrap_plots(plots, ncol = 3)

residuals <- read.table(file.path(path_temp_results, "Results_Residuals.txt"), header = TRUE)

library(ggplot2)

ggplot(residuals, aes(x = X3_obs)) +
    geom_point(aes(y = Y1_obs, col = "z_obs")) +
    geom_point(aes(y = Y1_sim, col = "z_sim")) +
    facet_wrap(~X1_obs, scales = "free_y")
