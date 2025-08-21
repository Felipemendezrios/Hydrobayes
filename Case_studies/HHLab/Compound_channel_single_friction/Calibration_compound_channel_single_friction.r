rm(list = ls())
graphics.off()

function_list <- list.files("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/Functions/", full.names = TRUE)
for (i in function_list) {
    source(i)
}

# Processed data
load("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/processed_data/HHLab/compound_single_friction/data_HHLab_uniform_case.RData")

# Libraries
library(RBaM)
library(dplyr)
library(patchwork)
library(tidyr)
library(ggplot2)

# setPathToBaM("/home/famendezrios/Documents/Git/BaM/makefile/")

# Experiment cases:

# 1_WSE_floodplain_low_uncertainty: a single event with WSE in the floodplain with low uncertainty values
# 1_WSE_floodplain_realistic_uncertainty: a single event with WSE in the floodplain with realistic uncertainty values
# 1_WSE_floodplain_1_Kmin_low_uncertainty: a single event with WSE in the floodplain and one Kmin tight
# 1_WSE_floodplain_several_Kmin: a single event with WSE in the floodplain and several Kmin with a realistic uncertainties values
# 1_WSE_floodplain_1_Kmoy_low_uncertainty: a single event with WSE in the floodplain and one Kmoy tight

all_experiments <- c(
    "1_WSE_floodplain_low_uncertainty",
    "1_WSE_floodplain_realistic_uncertainty",
    "1_WSE_floodplain_1_Kmin_low_uncertainty",
    "1_WSE_floodplain_several_Kmin",
    "1_WSE_floodplain_1_Kmoy_low_uncertainty"
)

dir_exe_BaM <- "/home/famendezrios/Documents/Git/BaM/makefile/"
do_calibration <- FALSE

for (Experiment_id in all_experiments) {
    # Path to save results
    path_temp_results <- file.path("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/Case_studies/HHLab/Compound_channel_single_friction/Calibration", Experiment_id)

    if (!dir.exists(path_temp_results)) {
        dir.create(path_temp_results)
    } else {
        if (do_calibration) {
            file.remove(list.files(path_temp_results, full.names = TRUE, recursive = TRUE))
        }
    }
    # Name Mage projet
    mage_projet_name <- "compound_channel"

    # Measurements for calibration
    Cal_measures <- c(3.5, 5, 6.5, 7, 7.5, 9.25, 11.75, 12.5, 13)

    WSE_data_temp <- all_data_calibration$WSE

    WSE_data <- data.frame(WSE_data_temp %>%
        group_by(x) %>%
        summarise(
            z_mean_all = mean(z_mean),
            Yu_mean_all = mean(Yu),
            z_riverbed = mean(z_riverbed)
        ))

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
    n_degree_Kflood <- 0

    # Input parameters in main channel:
    # Reach 1 : main channel
    a0_min_reach_1 <- RBaM::parameter(
        name = "a0_min",
        init = 1 / 0.010, # Initial guess
        prior.dist = "FlatPrior+", # Prior distribution
        prior.par = NULL
    ) # Parameters of the prior distribution

    # reach 1 : floodplain
    a0_flood_reach_1 <- RBaM::parameter(
        name = "a0_flood",
        init = 33,
        prior.dist = "FlatPrior+",
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
    } else {
        stop("case not supported")
    }


    if (length(Kmin_prior) != (Nb_reaches_estimation * (n_degree_Kmin + 1))) stop(paste0("More prior information (", length(Kmin_prior), ") than the number of reaches for estimation (", Nb_reaches_estimation, ")"))

    theta_param <- c(Kmin_prior, Kflood_prior)

    mageDir <- c(paste0("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/Case_studies/HHLab/Compound_channel_single_friction/model_mage/", Experiment_id, "/"))
    ### interpolation

    specific_points_Model <- data.frame(
        reach = c(1),
        KP_start = c(0),
        KP_end = c(18)
    )


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

    # Write Legendre vector:
    write.table(
        grid_covariant_discretized,
        file = file.path(path_temp_results, "vector_legendre.txt"), row.names = F
    )





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

    Yu <- data.frame(Yu_z = rep(0, length(CalData_sample$Yu_mean_all)))
    Yu$Yu_Q <- -9999
    Yu$Yu_V <- -9999
    Yu$Yu_Kmin <- -9999
    Yu$Yu_Kmoy <- -9999


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

    if (Experiment_id == "1_WSE_floodplain_low_uncertainty") {

    } else if (Experiment_id == "1_WSE_floodplain_realistic_uncertainty") {
        Yu$Yu_z <- CalData_sample$Yu_mean_all
    } else if (Experiment_id == "1_WSE_floodplain_1_Kmin_low_uncertainty") {
        Y$Y_Kmin[3] <- 100
        Yu$Yu_Kmin[3] <- 0

        remant_error_list[[4]] <- remnantErrorModel(
            fname = "Config_RemnantSigma4.txt",
            funk = "Constant",
            par = list(parameter(
                name = "intercept",
                init = 0.001,
                prior.dist = "LogNormal",
                prior.par = c(log(0.001), 1)
            ))
        )
    } else if (Experiment_id == "1_WSE_floodplain_several_Kmin") {
        Y$Y_Kmin[1:nrow(CalData_sample)] <- 100
        Yu$Yu_Kmin[1:nrow(CalData_sample)] <- 3

        remant_error_list[[4]] <- remnantErrorModel(
            fname = "Config_RemnantSigma4.txt",
            funk = "Constant",
            par = list(parameter(
                name = "intercept",
                init = 0.001,
                prior.dist = "LogNormal",
                prior.par = c(log(0.001), 1)
            ))
        )
    } else if (Experiment_id == "1_WSE_floodplain_1_Kmoy_low_uncertainty") {
        Y$Y_Kmoy[3] <- 33
        Yu$Yu_Kmoy[3] <- 0

        remant_error_list[[5]] <- remnantErrorModel(
            fname = "Config_RemnantSigma5.txt",
            funk = "Constant",
            par = list(parameter(
                name = "intercept",
                init = 0.001,
                prior.dist = "LogNormal",
                prior.par = c(log(0.001), 1)
            ))
        )
    } else {
        stop("Experiment_id not supported")
    }

    CalData_plot_df <- data.frame(
        x = X$x,
        z_mean = Y$WSE * 1000,
        Yu = Yu$Yu_z * 1000,
        z_riverbed = CalData_sample$z_riverbed * 1000,
        ID = "Averaged obs"
    )

    CalData_plot_export <- CalData_plot(
        data = CalData_plot_df,
        scales_free = "free_y",
        y_label = "Water surface elevation (mm)",
        title_label = "WSE observations",
        col_label = NULL,
        plot_water_depth = FALSE
    )

    ggsave(file.path(path_temp_results, "CalData_plot.png"),
        plot = CalData_plot_export
    )

    data <- dataset(X = X, Y = Y, Yu = Yu, data.dir = file.path(path_temp_results))

    prior_error_model <- get_init_prior(remant_error_list)

    jump_MCMC_theta_param_user <- 8
    jump_MCMC_error_model_user <- 0.001
    threshold_jump_MCMC_error_model <- 0.5

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
        run = FALSE,
        preClean = FALSE,
        workspace = path_temp_results,
        # dir.exe = file.path(find.package("RBaM"), "bin"),
        # name.exe = "BaM",
        predMaster_fname = "Config_Pred_Master.txt"
    )

    # Move Config_BaM.txt file to the result folder
    dir_cf <- file.path(path_temp_results, "Config_BaM.txt")

    if (do_calibration) {
        system2(
            command = file.path(dir_exe_BaM, "BaM"),
            args = c("-cf", file.path(path_temp_results, "Config_BaM.txt")),
            wait = FALSE
        )
    }
}

for (Experiment_id in all_experiments) {
    path_temp_results <- file.path("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/Case_studies/HHLab/Compound_channel_single_friction/Calibration", Experiment_id)

    ### PLOTS
    mcmc_not_cooked <- readMCMC(file.path(path_temp_results, "Results_MCMC.txt"))
    plots <- tracePlot(mcmc_not_cooked)

    mcmcplot <- wrap_plots(plots, ncol = 3)

    ggsave(
        file.path(
            path_temp_results,
            "MCMC_not_cooked.png"
        ),
        mcmcplot
    )

    # Density plot for each parameter
    plots <- densityPlot(mcmc_not_cooked)
    pdf_plot <- wrap_plots(plots, ncol = 3)

    ggsave(
        file.path(
            path_temp_results,
            "densityplot_not_cooked.png"
        ),
        pdf_plot
    )

    mcmc <- readMCMC(file.path(path_temp_results, "Results_Cooking.txt"))
    plots <- tracePlot(mcmc)

    mcmcplot <- wrap_plots(plots, ncol = 3)


    ggsave(
        file.path(
            path_temp_results,
            "MCMC_Cooked.png"
        ),
        mcmcplot
    )

    # Density plot for each parameter
    plots <- densityPlot(mcmc)
    pdf_plot <- wrap_plots(plots, ncol = 3)

    ggsave(
        file.path(
            path_temp_results,
            "densityplot_Cooked.png"
        ),
        pdf_plot
    )


    png(
        file.path(path_temp_results, "corelation_cooked.png"),
        width = 800,
        height = 800,
        res = 120
    )
    pairs(mcmc)
    dev.off()

    getSummary <- read.table(
        file = file.path(path_temp_results, "Results_Summary.txt"),
        header = TRUE,
        stringsAsFactors = FALSE
    )

    # Values of error model in meter for WSE, in m3/s for discharge and m/s for velocity.
    knitr::kable(getSummary,
        align = "c"
    )


    # Zoom into the MAP and standard deviation of the error model
    getSummary_zoom <- getSummary[c(11, 16), ]
    getSummary_zoom[, c("Y1_intercept")] <- getSummary_zoom[, c("Y1_intercept")] * 1000 # WSE in mm

    # friction for glass
    ks_literature_glass <- data.frame(min = 1 / 0.013, max = 1 / 0.009, mean = 1 / 0.010)

    # Get MAP simulation
    MAP_param_matrix <- as.numeric(getSummary_zoom[2, c(1:(n_degree_Kmin + 1 + 1 + n_degree_Kflood))])

    matrix_zFileKmin <- read.table(file.path(path_temp_results, "Zfile_Kmin.txt"), header = TRUE)

    position <- read.table(file.path(path_temp_results, "vector_legendre.txt"), header = TRUE)

    kmin_plot <- k_plot(
        matrix_spatialisation = matrix_zFileKmin,
        mcmc = mcmc,
        n_degree_kmin = n_degree_Kmin,
        n_degree_kflood = n_degree_kflood,
        covariate_discretization = position$Covariate,
        MAP_param_vector = MAP_param_matrix,
        main_channel = TRUE,
        ks_literature = ks_literature_glass
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

    kmoy_plot <- k_plot(
        matrix_spatialisation = matrix_zFileKflood,
        mcmc = mcmc,
        n_degree_kmin = n_degree_Kmin,
        n_degree_kflood = n_degree_Kflood,
        covariate_discretization = position$Covariate,
        MAP_param_vector = MAP_param_matrix,
        main_channel = FALSE,
        ks_literature = ks_literature_floodplain
    )

    ggsave(
        file.path(
            path_temp_results,
            "kmoy.png"
        ),
        kmoy_plot
    )

    # Read residuals
    residuals <- read.table(
        file = file.path(path_temp_results, "Results_Residuals.txt"),
        header = TRUE,
        stringsAsFactors = FALSE
    )

    # Convert residuals to mm and m3/s
    residuals_mm_m3_s <- data.frame(
        residuals[, 1:4],
        residuals[, c(9, 19, 20, 24, 29)] * 1000
    )

    CalData <- read.table(file.path(path_temp_results, "CalibrationData.txt"), header = TRUE)


    residuals_event_mm_m3_s <- data.frame(residuals_mm_m3_s, Yu_z = CalData$Yu_z * 1000)

    Q_residuals <- ZQdX_residuals(
        residuals_event_mm_m3_s = residuals_event_mm_m3_s,
        Q_input = 114,
        Qplot = TRUE,
        title_label = "MAP simulations vs boundary condition \n(Discharge upstream)",
        ylabel = "Q (L/s)"
    )

    ggsave(
        file.path(
            path_temp_results,
            "Q_residuals.png"
        ),
        Q_residuals
    )

    Z_residuals <- ZQdX_residuals(
        residuals_event_mm_m3_s = residuals_event_mm_m3_s,
        Q_input = NULL,
        Qplot = FALSE,
        title_label = "MAP simulations vs Observations (WSE)",
        ylabel = "Water surface elevation (mm)"
    )

    ggsave(
        file.path(
            path_temp_results,
            "Z_residuals.png"
        ),
        Z_residuals
    )
}
