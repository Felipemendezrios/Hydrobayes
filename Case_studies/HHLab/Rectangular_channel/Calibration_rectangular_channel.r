rm(list = ls())
graphics.off()

function_list <- list.files("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/Functions/", full.names = TRUE)
for (i in function_list) {
    source(i)
}

# Processed data
load("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/processed_data/HHLab/smooth_bed/data_HHLab_all_cases.RData")

# Libraries
library(RBaM)
library(dplyr)
library(patchwork)
library(tidyr)
library(ggplot2)

# Common for observation and calibration
all_experiments <- c(
    "4_WSE_main_channel_real_uncertainty"
)
# Folder related to the polynomial degree calibration
all_polynomial_degree <- c(
    "n_0",
    "n_1",
    "n_2",
    "n_3"
)

# Folder related to the observations (careful with the order!)
all_events <- c(
    "Q_120",
    "Q_60",
    "Q_30",
    "Q_14"
)


Do_Correction_Cross_section <- TRUE

dir_exe_BaM <- "/home/famendezrios/Documents/Git/BaM/makefile/"
do_calibration <- FALSE

# Observations data input:
# Measurements for calibration by event!
# Ensure that order should be in coherence with mage model set up!
Cal_measures <-
    list(
        Q_120 =
            data.frame(
                event = 1,
                reach = 1,
                x = c(
                    0.060,
                    0.650,
                    1.650,
                    3.650,
                    5.650,
                    7.650,
                    9.650,
                    11.650
                ),
                t = rep(3420, 8)
            ),
        Q_60 = data.frame(
            event = 2,
            reach = 1,
            x = c(
                0.060,
                0.650,
                1.650,
                3.650,
                5.650,
                7.650,
                9.650,
                11.650
            ),
            t = rep(3420, 8)
        ),
        Q_30 = data.frame(
            event = 3,
            reach = 1,
            x = c(
                0.060,
                0.650,
                1.650,
                3.650,
                5.650,
                7.650,
                9.650,
                11.650
            ),
            t = rep(3420, 8)
        ),
        Q_14 = data.frame(
            event = 4,
            reach = 1,
            x = c(
                0.060,
                0.650,
                1.650,
                3.650,
                5.650,
                7.650,
                9.650,
                11.650
            ),
            t = rep(3420, 8)
        )
    )

# All available data
WSE_data_temp <- all_data_calibration$WSE

# Only giving the name keeping the data order and creating a link with the names of the event and the Cal_measures defined previously
match_case_event <- data.frame(
    idx_data = names(WSE_data_temp),
    idx_event = c("Q_14", "Q_30", "Q_60", "Q_120")
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

# Initialize an empty list to store the results
observed_data <- data.frame()
X <- c()

# Loop over each element in Link_x_t_ind_event
for (i in seq_along(Link_x_t_ind_event)) {
    idx_data <- Link_x_t_ind_event[[i]]$idx_data
    cal_x <- Link_x_t_ind_event[[i]]$Cal_measure$x

    # Get the corresponding data frame from WSE_data_temp
    wse_df <- WSE_data_temp[[idx_data]]

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

Y <- data.frame(WSE = observed_data$z_mean)
Y$Discharge <- -9999
Y$Velocity <- -9999
Y$Y_Kmin <- -9999
Y$Y_Kmoy <- -9999

Yu <- data.frame(Yu_z = observed_data$Yu)
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
########

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


## Calibration setting

# Input parameters in main channel:
# Reach 1 : main channel
a0_min_reach_1 <- RBaM::parameter(
    name = "a0_min",
    init = 1 / 0.010, # Initial guess
    prior.dist = "FlatPrior+", # Prior distribution
    prior.par = NULL
) # Parameters of the prior distribution

a1_min_reach_1 <- RBaM::parameter(
    name = "a1_min",
    init = 0, # Initial guess
    prior.dist = "FlatPrior", # Prior distribution
    prior.par = NULL
)

a2_min_reach_1 <- RBaM::parameter(
    name = "a2_min",
    init = 0, # Initial guess
    prior.dist = "FlatPrior", # Prior distribution
    prior.par = NULL
)

a3_min_reach_1 <- RBaM::parameter(
    name = "a3_min",
    init = 0, # Initial guess
    prior.dist = "FlatPrior", # Prior distribution
    prior.par = NULL
)

# reach 1 : floodplain
a0_flood_reach_1 <- RBaM::parameter(
    name = "a0_flood",
    init = 10,
    prior.dist = "FIX",
    prior.par = NULL
)

n_degree_Kflood <- 0

Kflood_prior <- list(
    a0_flood_reach_1
)


## Cross-section treatment:
### interpolation

specific_points_Model <- data.frame(
    reach = c(1),
    KP_start = c(0),
    KP_end = c(16)
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

# Common name for all mage projects
mage_projet_name <- "rectangular_channel"


##################################
## Do not touch !
# Harcode values, it does not matter, it just useful for give the size of the .RUG File
RUG_Kmin <- 20
RUG_Kmoy <- 10

jump_MCMC_theta_param_user <- 8
jump_MCMC_error_model_user <- 0.001
threshold_jump_MCMC_error_model <- 0.5
prior_error_model <- get_init_prior(remant_error_list)
##################################

for (Experiment_id in all_experiments) {
    if (Do_Correction_Cross_section) {
        path_Cross_section_to_remplace <- c(paste0("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/Case_studies/HHLab/Rectangular_channel/Calibration_experiments/", Experiment_id, "/", all_polynomial_degree, "/model_mage/", all_events[1], "/net/Reach_002.ST"))

        for (temp_path in path_Cross_section_to_remplace) {
            Cross_sections_interpolation_and_export_rectangular_channel(
                So = 1.04e-3,
                Pk = c(
                    seq(0, 16.5040, by = 1)
                ),
                positions_banks = c(0, 2),
                borders_heigh = 1,
                path_export = temp_path
            )
        }
    }

    path_temp_experiment <- file.path("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/Case_studies/HHLab/Rectangular_channel/Calibration_experiments", Experiment_id)

    ggsave(file.path(path_temp_experiment, "CalData_plot.png"),
        plot = CalData_plot_export,
        width = 20,
        height = 20,
        units = "cm"
    )

    for (polynomial_id in all_polynomial_degree) {
        # That will work only for exploring several polynomial degrees in the main chanil, fixing the floodplain
        if (polynomial_id == "n_0") {
            n_degree_Kmin <- 0
            Kmin_prior <- list(
                a0_min_reach_1
            )
        } else if (polynomial_id == "n_1") {
            n_degree_Kmin <- 1
            Kmin_prior <- list(
                a0_min_reach_1,
                a1_min_reach_1
            )
        } else if (polynomial_id == "n_2") {
            n_degree_Kmin <- 2
            Kmin_prior <- list(
                a0_min_reach_1,
                a1_min_reach_1,
                a2_min_reach_1
            )
        } else if (polynomial_id == "n_3") {
            n_degree_Kmin <- 3
            Kmin_prior <- list(
                a0_min_reach_1,
                a1_min_reach_1,
                a2_min_reach_1,
                a3_min_reach_1
            )
        } else {
            stop("polynomial_id not supported yet")
        }

        path_polynomial <- file.path(path_temp_experiment, polynomial_id)
        # Path to save results
        path_temp_results <- file.path(path_polynomial, "Calibration")

        if (!dir.exists(path_temp_results)) {
            dir.create(path_temp_results)
        } else {
            if (do_calibration) {
                file.remove(list.files(path_temp_results, full.names = TRUE, recursive = TRUE))
            }
        }


        if (length(Kmin_prior) != (Nb_reaches_estimation * (n_degree_Kmin + 1))) stop(paste0("More prior information (", length(Kmin_prior), ") than the number of reaches for estimation (", Nb_reaches_estimation, ")"))

        theta_param <- c(Kmin_prior, Kflood_prior)

        # Write Legendre vector:
        write.table(
            grid_covariant_discretized,
            file = file.path(path_temp_results, "legendre_covariate_non_normalized.txt"), row.names = F
        )

        mageDir <- c(paste0(path_polynomial, "/model_mage/", all_events, "/"))

        RUGFiles <- paste0(mageDir, mage_projet_name, ".RUG")

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

        prior_theta_param <- get_init_prior(theta_param)
        prior_all_dist_theta_param <- get_init_prior(theta_param, FIX_dist = TRUE)

        data <- dataset(X = X, Y = Y, Yu = Yu, data.dir = file.path(path_temp_results))

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

        dir_cf <- file.path(path_temp_results, "Config_BaM.txt")

        if (do_calibration) {
            system2(
                command = file.path(dir_exe_BaM, "BaM"),
                args = c("-cf", file.path(path_temp_results, "Config_BaM.txt")),
                wait = FALSE
            )
        }
    }
}


for (Experiment_id in all_experiments) {
    for (polynomial_id in all_polynomial_degree) {
        # That will work only for exploring several polynomial degrees in the main chanil, fixing the floodplain
        if (polynomial_id == "n_0") {
            n_degree_Kmin <- 0
            Kmin_prior <- list(
                a0_min_reach_1
            )
        } else if (polynomial_id == "n_1") {
            n_degree_Kmin <- 1
            Kmin_prior <- list(
                a0_min_reach_1,
                a1_min_reach_1
            )
        } else if (polynomial_id == "n_2") {
            n_degree_Kmin <- 2
            Kmin_prior <- list(
                a0_min_reach_1,
                a1_min_reach_1,
                a2_min_reach_1
            )
        } else if (polynomial_id == "n_3") {
            n_degree_Kmin <- 3
            Kmin_prior <- list(
                a0_min_reach_1,
                a1_min_reach_1,
                a2_min_reach_1,
                a3_min_reach_1
            )
        } else {
            stop("polynomial_id not supported yet")
        }

        path_temp_plots <- file.path("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/Case_studies/HHLab/Rectangular_channel/Calibration_experiments", Experiment_id, polynomial_id, "Calibration/")

        ### PLOTS
        mcmc_not_cooked <- readMCMC(file.path(path_temp_plots, "Results_MCMC.txt"))
        plots <- tracePlot(mcmc_not_cooked)

        mcmcplot <- wrap_plots(plots, ncol = 3)

        ggsave(
            file.path(
                path_temp_plots,
                "MCMC_not_cooked.png"
            ),
            mcmcplot
        )

        # Density plot for each parameter
        plots <- densityPlot(mcmc_not_cooked)
        pdf_plot <- wrap_plots(plots, ncol = 3)

        ggsave(
            file.path(
                path_temp_plots,
                "densityplot_not_cooked.png"
            ),
            pdf_plot
        )

        mcmc <- readMCMC(file.path(path_temp_plots, "Results_Cooking.txt"))
        plots <- tracePlot(mcmc)

        mcmcplot <- wrap_plots(plots, ncol = 3)


        ggsave(
            file.path(
                path_temp_plots,
                "MCMC_Cooked.png"
            ),
            mcmcplot
        )

        # Density plot for each parameter
        plots <- densityPlot(mcmc)
        pdf_plot <- wrap_plots(plots, ncol = 3)

        ggsave(
            file.path(
                path_temp_plots,
                "densityplot_Cooked.png"
            ),
            pdf_plot
        )


        png(
            file.path(path_temp_plots, "corelation_cooked.png"),
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
        getSummary_zoom[, c("Y1_intercept")] <- getSummary_zoom[, c("Y1_intercept")] * 1000 # WSE in mm

        # friction for glass
        ks_literature_glass <- data.frame(min = 1 / 0.013, max = 1 / 0.009, mean = 1 / 0.010)

        # Get MAP simulation
        MAP_param_matrix <- as.numeric(getSummary_zoom[2, c(1:(n_degree_Kmin + 1 + 1 + n_degree_Kflood))])

        matrix_zFileKmin <- read.table(file.path(path_temp_plots, "Zfile_Kmin.txt"), header = TRUE)

        position <- read.table(file.path(path_temp_plots, "legendre_covariate_non_normalized.txt"), header = TRUE)

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
                path_temp_plots,
                "kmin.png"
            ),
            kmin_plot
        )

        # ks_literature_floodplain <- data.frame(min = 1 / 0.05, max = 1 / 0.025, mean = 1 / 0.033)

        # matrix_zFileKflood <- read.table(file.path(path_temp_plots, "Zfile_Kflood.txt"), header = TRUE)

        # kmoy_plot <- k_plot(
        #     matrix_spatialisation = matrix_zFileKflood,
        #     mcmc = mcmc,
        #     n_degree_kmin = n_degree_Kmin,
        #     n_degree_kflood = n_degree_Kflood,
        #     covariate_discretization = position$Covariate,
        #     MAP_param_vector = MAP_param_matrix,
        #     main_channel = FALSE,
        #     ks_literature = ks_literature_floodplain
        # )

        # ggsave(
        #     file.path(
        #         path_temp_plots,
        #         "kmoy.png"
        #     ),
        #     kmoy_plot
        # )

        # Read residuals
        residuals <- read.table(
            file = file.path(path_temp_plots, "Results_Residuals.txt"),
            header = TRUE,
            stringsAsFactors = FALSE
        )

        # Convert residuals to mm and m3/s
        residuals_mm_m3_s <- data.frame(
            residuals[, 1:4],
            residuals[, c(9, 19, 20, 24, 29)] * 1000
        )

        CalData <- read.table(file.path(path_temp_plots, "CalibrationData.txt"), header = TRUE)


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
                path_temp_plots,
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
                path_temp_plots,
                "Z_residuals.png"
            ),
            Z_residuals
        )
    }
}
