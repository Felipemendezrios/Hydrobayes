rm(list = ls())
graphics.off()

# Set directory (root of the repository)
dir_workspace <- "/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git"
setwd(dir_workspace)

function_list <- list.files("R", full.names = TRUE)
for (i in function_list) {
    source(i)
}

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
library(utils)
#############################################
# End load libraries
#############################################

############################################
# Module 0 : General settings
############################################

# Logical
do_calibration <- FALSE
do_plot_calibration <- FALSE
do_prediction <- FALSE


# Setting jump standard deviation for MCMC sampling if initial guess = 0
jump_MCMC_theta_param_user_regression <- 5
jump_MCMC_theta_param_user_coeff <- 0.1
jump_MCMC_error_model_user <- 0.001
threshold_jump_MCMC_error_model <- 0.5

############################################
# Module 1 : set directory paths
############################################

# Name of the experiment. All scenarios will be used the same calibration data

# Common for observation and calibration folders: names of the experiment
Experiment_id <- c(
    "4_WSE_3_Q_Seine_2015"
)

# Calibration case: SU distribution
SU_distribution <- "1SU"


# Experiments input data to be used during calibration setting
all_cal_case <- c(
    "Kmin_SU1_n6_4WSE_3Q.r",
    "Kmin_SU1_n7_4WSE_3Q.r",
    "Kmin_SU1_n8_4WSE_3Q.r",
    "Kmin_SU1_n9_4WSE_3Q.r",
    "Kmin_SU1_n10_4WSE_3Q.r",
    "Kmin_SU1_n11_4WSE_3Q.r",
    "Kmin_SU1_n12_4WSE_3Q.r"
)


# Folder related to the observations (careful with the order!)
all_events <- c(
    "Piney_2015"
)
# Date of reference
date_ref <- as.POSIXct(
    "2015-09-28 01:10:00",
    format = "%Y-%m-%d %H:%M:%S",
    tz = "UTC"
)

command_line_MAGE <- ""


file_main_path <- file.path(dir_workspace, "scripts/Case_studies/Lower_Seine", SU_distribution, "Calibration_experiments")
MAGE_main_folder <- file.path(dir_workspace, "scripts/Case_studies/Lower_Seine", SU_distribution, "model_mage")

mage_projet_name <- "Lower_Seine"

check_mage_folder(name_folder = MAGE_main_folder)
############################################
# End module 1: set directories paths
############################################


############################################
# Module 2 : hydraulic model (HM) environment
############################################
# Input information about reaches and boundaries KP of the HM.
# Model_Reach : Model reach interpreted by the HM.
# This input information must be coherent with HM specifications.

# KP start and end must be reliable to the model. To define KP range, it does not matter the direction, it must be coherent to the HM

Input_Model_Reach <- data.frame(
    reach = c(1, 2),
    KP_start = c(0, 5064),
    KP_end = c(5064, 137983)
)

# First layer: Typology
# ###########################################
# All reaches in the model (Model_Reach) defined within the same coordinate reference. Useful to separate main reach to tributary or even a diversion channel

# Must be careful with the order of the reaches, it must be given upstream to downstream

Input_Typology <- list(
    Seine = c(1, 2)
)

############################################
# End module 2 : hydraulic model (HM) environment
############################################


############################################
# Module 3: calibration data
############################################

# Processed data
load("data/processed_data/Lower_Seine/Seine_2015/observed_data.RData")


results_CalData <- constructor_CalData_latest(
    observed_data = observed_data,
    date_ref = date_ref,
    do_manual_uncertainty = FALSE,
    variables = c("WSE", "Q", "V", "Kmin", "Kflood")
)

X_plot <- data.frame(results_CalData$X)

X <- X_plot[, !names(X_plot) %in% c("target_datetime", "var", "id_campaign", "date_time_format")]
Y <- data.frame(results_CalData$Y)
Yu <- data.frame(results_CalData$Yu)

CalData <- cbind(X, Y, Yu)

CalData_plot <- cbind(X_plot, Y, Yu)

path_experiment <- file.path(file_main_path, Experiment_id)

station_KP <- data.frame(
    station = c("Poses", "Elbeuf", "Oissel", "Rouen", "Petit_Couronne", "La_Bouille", "Val_des_leux", "Duclair", "Mesnil_sous_Jumieges", "Heurteauville", "Caudebec", "Vatteville", "Aizier", "Saint_Leonard", "Tancarville"),
    KP = c(0, 16823, 27439, 41568, 51186, 58878, 65055, 78482, 85578, 97091, 111227, 116862, 123562, 130623, 137983),
    id_reach_HM = c(1, rep(2, 14))
)


if (do_plot_calibration) {
    # Customized the order of the reach (HM) for plotting
    CalData_ordered <- CalData_plot %>%
        mutate(reach = factor(reach,
            levels = c(2)
        ))
    plots_CalData <- plot_CalData_lastest(
        CalData = CalData_ordered,
        scales_free = "free",
        wrap = FALSE
    )

    # Add bathymetry
    bathymetry_upstream <- read.table("data/data-raw/Lower_Seine/Bathymetry/ST_Bathy_upstream.txt", header = TRUE, sep = "")
    bathymetry_downstream <- read.table("data/data-raw/Lower_Seine/Bathymetry/ST_Bathy_downstream.txt", header = TRUE, sep = "")

    all_bathy <- rbind(
        bathymetry_upstream,
        bathymetry_downstream
    )

    target_dates <- CalData_ordered %>%
        mutate(WSE = convert_9999_to_NA(WSE)) %>%
        filter(!is.na(WSE)) %>%
        distinct(event, target_datetime)

    profiles <- all_bathy %>%
        mutate(reach = ifelse(X <= 5064, 1, 2)) %>%
        group_by(profile_id, reach) %>%
        slice_min(order_by = Z, n = 1, with_ties = FALSE) %>%
        summarise(
            KP = first(X),
            Z_thalweg = first(Z),
            .groups = "drop" # This ungroups after summarising
        )

    model_profile_id_reach_KP <- profiles %>%
        crossing(target_dates) %>%
        filter(reach %in% unique(CalData_ordered$reach))

    plot_WSE_Thalweg <- plots_CalData$plot_WSE +
        geom_line(data = model_profile_id_reach_KP, aes(x = KP, y = Z_thalweg), color = "black") +
        facet_wrap(~event,
            labeller = labeller(
                event = as_labeller(
                    c(
                        "1" = "Event~1:~2015"
                    ),
                    label_parsed
                )
            )
        ) +
        scale_x_continuous(labels = scales::label_number())

    plots_CalData$plot_WSE_Thalweg <- plot_WSE_Thalweg

    # Customize the plot
    plots_CalData$plot_WSE <- plots_CalData$plot_WSE +
        facet_wrap(
            ~ event + reach + target_datetime,
            labeller = labeller(
                event = as_labeller(
                    c(
                        "1" = "Event~1:~2015"
                    ),
                    label_parsed
                ),
                reach = as_labeller(c(
                    "2" = "Reach~2:~Eure~-~Tancarville"
                ), label_parsed)
            ),
            scales = "free",
            ncol = 2
        ) + scale_x_continuous(labels = scales::label_number())

    # Customize the plot
    station_labeller <- function(labels) {
        labels$reach <- sapply(labels$reach, function(r) {
            kp <- if (r == 1) {
                min(station_KP$KP[station_KP$id_reach_HM == 1])
            } else {
                max(station_KP$KP[station_KP$id_reach_HM == 2])
            }

            station_KP$station[which.min(abs(station_KP$KP - kp))]
        })

        labels
    }

    plots_CalData$plot_Q <- plots_CalData$plot_Q +
        facet_wrap(
            ~ event + reach + x,
            labeller = labeller(
                event = as_labeller(
                    c(
                        "1" = "Event~1:~2015"
                    ),
                    label_parsed
                ),
                reach = as_labeller(c(
                    "2" = "Reach~2:~Eure~-~Tancarville"
                ), label_parsed),
                x = function(x) {
                    sapply(x, function(kp) {
                        station_KP$station[
                            which.min(abs(station_KP$KP - as.numeric(kp)))
                        ]
                    })
                }
            ),
            scales = "free",
            ncol = 2
        )

    if (!dir.exists(path_experiment)) {
        dir.create(path_experiment)
    }

    for (name in names(plots_CalData)) {
        plot <- plots_CalData[[name]]
        if (!is.null(plot)) {
            ggsave(file.path(path_experiment, paste0(name, ".png")),
                plot = plot,
                width = 35,
                height = 20,
                units = "cm"
            )

            save(
                file = file.path(path_experiment, paste0(name, ".RData")),
                plot
            )
        }
    }
}
############################################
# End module 2: calibration data
############################################

############################################
# Module 4: calibration setting
############################################

######################################################
# Set Remnant error model
# Structural error information
# Put remnantErrorModel_default on the output variable without calibration data
remant_error_list <- list(
    # WSE
    RBaM::remnantErrorModel(
        fname = "Config_RemnantSigma.txt",
        funk = "Constant",
        par = list(parameter(
            name = "intercept",
            init = 0.25,
            prior.dist = "Uniform",
            prior.par = c(-10, 10)
        ))
    ),
    # Q
    RBaM::remnantErrorModel(
        fname = "Config_RemnantSigma2.txt",
        funk = "Constant",
        par = list(parameter(
            name = "intercept",
            init = 500,
            prior.dist = "Uniform",
            prior.par = c(-1000, 1000)
        ))
    ),
    # V
    remnantErrorModel_default(name = "Config_RemnantSigma3.txt"),
    # Both Kmin and Kflood, they must always be forced to 0. Information will be passed by pseudo-obs
    # Kmin
    remnantErrorModel_default(name = "Config_RemnantSigma4.txt"),
    # Kflood
    remnantErrorModel_default(name = "Config_RemnantSigma5.txt")
)

# Get initial prior of structural error
prior_error_model <- get_init_prior(remant_error_list)

############################################
# End module 4: calibration setting
############################################



############################################
# Module 5: BaM environment
############################################

Key_Info_Typology_Model_Reach <- get_Key_Info_Typology_Model_Reach(
    Input_Typology = Input_Typology,
    Input_Model_Reach = Input_Model_Reach,
    total_points_discretization = 200
)

############################################
# Module 6: Calibration
############################################
list_mod_polynomials <- list_Z_MatrixKmin <- list_Z_MatrixKflood <- list_Kmin_prior <- list_Kflood_prior <- list_Kmin_SU <- list_Kflood_SU <- list_ref_Matrix_Prior_Correlation <- list_summary_SU_Kflood <- list_summary_SU_Kmin <- list()

for (id_cal_case in 1:length(all_cal_case)) {
    # Load experiment
    paths <- load_experiment(
        file_main_path = file_main_path,
        cal_case = all_cal_case[[id_cal_case]],
        path_experiment = path_experiment,
        all_events = all_events
    )

    results_estimation <- Estimation_Mage(
        paths = paths,
        Key_Info_Typology_Model_Reach = Key_Info_Typology_Model_Reach,
        Input_Typology = Input_Typology,
        MAGE_main_folder = MAGE_main_folder,
        do_calibration = do_calibration,
        command_line_MAGE = command_line_MAGE,
        ID_model_BaM = ID_model_BaM,
        nX_BaM = nX_BaM,
        nY_BaM = nY_BaM,
        mage_projet_name = mage_projet_name,
        mcmcCooking = RBaM::mcmcCooking(burn = 0.5, nSlim = 10),
        mcmcOptions = RBaM::mcmcOptions(nAdapt = 60, nCycles = 60),
        mcmcSummary = RBaM::mcmcSummary(xtendedMCMC.fname = "Results_xtendedMCMC.txt"),
        remant_error_list = remant_error_list
    )
    script_path <- file.path(paths$path_BaM_folder, "run_BaM.sh")

    writeLines(
        c(
            "#!/bin/bash",
            paste(
                shQuote(file.path(RBaM::getPathToBaM(), "BaM")),
                "-cf",
                shQuote(file.path(paths$path_BaM_folder, "Config_BaM.txt"))
            )
        ),
        script_path
    )

    Sys.chmod(script_path, "0755")
    # Find path pstree -ap | grep BaM
    # Run outside of Vscodium. To kill a job : pkill -f BaM

    if (do_calibration) {
        system2(
            "nohup",
            args = c("bash", script_path),
            stdout = FALSE,
            stderr = FALSE,
            wait = FALSE
        )
    }

    Kmin_SU_RData <- results_estimation$Kmin_SU
    Kflood_SU_RData <- results_estimation$Kflood_SU
    save(Kmin_SU_RData,
        file = file.path(paths$path_RData, "Kmin_SU_RData.RData")
    )
    save(Kflood_SU_RData,
        file = file.path(paths$path_RData, "Kflood_SU_RData.RData")
    )

    list_Z_MatrixKmin[[id_cal_case]] <- results_estimation$Z_MatrixKmin
    list_Z_MatrixKflood[[id_cal_case]] <- results_estimation$Z_MatrixKflood
    list_Kmin_prior[[id_cal_case]] <- results_estimation$Kmin_prior
    list_Kflood_prior[[id_cal_case]] <- results_estimation$Kflood_prior
    list_Kmin_SU[[id_cal_case]] <- results_estimation$Kmin_SU
    list_Kflood_SU[[id_cal_case]] <- results_estimation$Kflood_SU
    list_mod_polynomials[[id_cal_case]] <- results_estimation$mod
    list_summary_SU_Kmin[[id_cal_case]] <- results_estimation$summary_SU_Kmin
    list_summary_SU_Kflood[[id_cal_case]] <- results_estimation$summary_SU_Kflood
    list_ref_Matrix_Prior_Correlation[[id_cal_case]] <- results_estimation$ref_Matrix_Prior_Correlation
}

# Plot DIC
if (do_plot_calibration) {
    plotDIC <- plot_DIC(dir_polynomial = c(file.path(
        path_experiment, sub("\\.r$", "", all_cal_case)
    )))
    ggsave(
        file.path(
            path_experiment,
            "DIC.png"
        ),
        plotDIC,
        width = 30,
        height = 20,
        units = "cm"
    )

    save(
        plotDIC,
        file = file.path(
            path_experiment,
            "DIC.RData"
        )
    )
}
synthetic_case <- FALSE

################################
# POSTPROCESS CALIBRATION WORKFLOW
################################
final_calibration <- FALSE

for (id_cal_case in 1:length(all_cal_case)) {
    # Load experiment
    paths <- load_experiment(
        file_main_path = file_main_path,
        cal_case = all_cal_case[[id_cal_case]],
        path_experiment = path_experiment,
        all_events = all_events
    )

    results_postprocess <- postprocess_calibration(
        paths = paths,
        X_input = X_plot,
        Y_observations = Y,
        Yu_observations = Yu,
        final_calibration = final_calibration,
        Key_Info_Typology_Model_Reach = Key_Info_Typology_Model_Reach,
        summary_SU_Kmin = list_summary_SU_Kmin[[id_cal_case]],
        summary_SU_Kflood = list_summary_SU_Kflood[[id_cal_case]],
        Kmin_prior = list_Kmin_prior[[id_cal_case]],
        Kflood_prior = list_Kflood_prior[[id_cal_case]],
        Kmin_SU = list_Kmin_SU[[id_cal_case]],
        Kflood_SU = list_Kflood_SU[[id_cal_case]],
        Z_MatrixKmin = list_Z_MatrixKmin[[id_cal_case]],
        Z_MatrixKflood = list_Z_MatrixKflood[[id_cal_case]],
        mod_polynomials = list_mod_polynomials[[id_cal_case]],
        Kmin_segment_layer = NULL,
        Kflood_segment_layer = NULL,
        command_line_MAGE = command_line_MAGE,
        dir_workspace = dir_workspace
    )

    CalData_updated <- results_postprocess$CalData_updated
    residuals <- results_postprocess$residuals
    plot_Kmin_without_obs <- results_postprocess$plots_param$Kmin$plot_without_obs

    plot_Kmin_with_obs <- results_postprocess$plots_param$Kmin$plot_with_obs

    if (any(!is.na(CalData_updated[, c("Kmin")]))) {
        plot_Kmin_with_obs <-
            plot_Kmin_with_obs +
            geom_point(
                data = CalData_updated,
                aes(x = x, y = Kmin, col = "Pseudo \nobs", group = id_reach_SU_Kmin),
                alpha = 0.2
                # inherit.aes = FALSE
            ) +
            geom_errorbar(
                data = CalData_updated,
                aes(x = x, ymin = Kmin - 1.96 * Yu_Kmin, ymax = Kmin + 1.96 * Yu_Kmin, col = "Pseudo \nobs", group = id_reach_SU_Kmin),
                alpha = 0.2
            ) +
            scale_color_manual(values = c(
                "MAP" = "black",
                "Pseudo \nobs" = "blue"
            ))
    }

    plot_Kflood_without_obs <- results_postprocess$plots_param$Kflood$plot_without_obs
    plot_Kflood_with_obs <- results_postprocess$plots_param$Kflood$plot_with_obs
    if (any(!is.na(CalData_updated[, c("Kflood")]))) {
        plot_Kflood_with_obs <- plot_Kflood_with_obs +
            geom_point(
                data = CalData_updated, aes(x = x, y = Kflood, col = "Pseudo \nobs", group = "id_reach_SU_Kflood"), alpha = 0.2,
            ) +
            geom_errorbar(data = CalData_updated, aes(x = x, ymin = Kflood - 1.96 * Yu_Kflood, ymax = Kflood + 1.96 * Yu_Kflood, col = "Pseudo \nobs", group = "id_reach_SU_Kflood"), alpha = 0.2) +
            scale_color_manual(values = c(
                "MAP" = "black",
                "Pseudo \nobs" = "blue"
            ))
    }

    plots_MAP_output_variables <- results_postprocess$plots_MAP_output_variables

    if (!is.null(plot_Kmin_with_obs)) {
        ggsave(
            filename = file.path(
                paths$path_plot_folder,
                paste0("plot_Kmin_with_true_values_generated_obs.png")
            ),
            plot = plot_Kmin_with_obs,
            width = 20,
            height = 20,
            units = "cm"
        )
    }
    if (!is.null(plot_Kflood_with_obs)) {
        ggsave(
            filename = file.path(
                paths$path_plot_folder,
                paste0("plot_Kflood_with_true_values_generated_obs.png")
            ),
            plot = plot_Kflood_with_obs,
            width = 20,
            height = 20,
            units = "cm"
        )
    }


    for (i in seq_along(plots_MAP_output_variables)) {
        ggsave(
            filename = file.path(
                paths$path_plot_folder,
                paste0("plot_obs_sim_MAP_Y", i, ".png")
            ),
            plot = plots_MAP_output_variables[[i]],
            width = 20,
            height = 20,
            units = "cm"
        )
    }
    save(plots_MAP_output_variables,
        file = file.path(paths$path_RData, "plots_MAP_output_variables.RData")
    )

    # Plot prior vs posterior
    # Parameter to parameter Kmin
    if (!is.null(results_postprocess$plots_prior_vs_posterior$param$Kmin)) {
        ggsave(
            filename = file.path(
                paths$path_plot_folder,
                paste0("plot_prior_vs_posterior_Kmin.png")
            ),
            plot = results_postprocess$plots_prior_vs_posterior$param$Kmin,
            width = 20,
            height = 20,
            units = "cm"
        )
    }
    # Parameter to parameter Kflood
    if (!is.null(results_postprocess$plots_prior_vs_posterior$param$Kflood)) {
        ggsave(
            filename = file.path(
                paths$path_plot_folder,
                paste0("plot_prior_vs_posterior_Kflood.png")
            ),
            plot = results_postprocess$plots_prior_vs_posterior$param$Kflood,
            width = 20,
            height = 20,
            units = "cm"
        )
    }
    # Parameter to parameter: K(x)
    if (!is.null(results_postprocess$plots_prior_vs_posterior$KdX$Kmin)) {
        ggsave(
            filename = file.path(
                paths$path_plot_folder,
                paste0("plot_prior_vs_posterior_Kmin_KdX.png")
            ),
            plot = results_postprocess$plots_prior_vs_posterior$KdX$Kmin,
            width = 20,
            height = 20,
            units = "cm"
        )
    }
    # Parameter to parameter: K(x)
    if (!is.null(results_postprocess$plots_prior_vs_posterior$KdX$Kflood)) {
        ggsave(
            filename = file.path(
                paths$path_plot_folder,
                paste0("plot_prior_vs_posterior_Kflood_KdX.png")
            ),
            plot = results_postprocess$plots_prior_vs_posterior$KdX$Kflood,
            width = 20,
            height = 20,
            units = "cm"
        )
    }

    # Specific case of synthetic case
    if (synthetic_case) {
        plot_output_with_synthetic_data <-
            plots_MAP_output_variables[[1]] +
            geom_point(
                data = real_synt_data,
                aes(x = KP, y = WSE_real_obs, col = "synthetic data", group = id_reach_CAL), shape = 2
            ) +
            scale_color_manual(
                values =
                    c(
                        "sim" = "black",
                        "obs" = "blue",
                        "synthetic data" = "purple"
                    )
            ) +
            facet_wrap(
                ~X1_obs,
                # labeller = labeller(
                #     X1_obs = c(
                #         "1" = "Main reach (MR)",
                #         "2" = "Tributary (TR)"
                #     )
                # ),
                scales = "free",
                ncol = 1
            )
        save(plot_output_with_synthetic_data,
            file = file.path(paths$path_RData, "plot_output_with_synthetic_data.RData")
        )
        ggsave(
            filename = file.path(
                paths$path_plot_folder,
                paste0("plot_output_with_synthetic_data_Y1.png")
            ),
            plot = plot_output_with_synthetic_data,
            width = 20,
            height = 20,
            units = "cm"
        )
    }
}
################################
# Prediction
################################

# Define the calculation grid:
# It must respect the information in the ST File for the space and Time simulation for the time

# ZdX: t fixed and x variable (SWOT, measure wse, drone, etc)
# QdXT: t fixed and x fixed (gaugings)
# VdXT: t fixed and x fixed (radar)
# ZdT: t variable, x fixed (hydrometric stations)

# Grid is possible to be spatial and temporal, in prediction step, this is not a problem.
# But i need to fixed time or space to vary the other one ! (could be improved in the future)

# In this case, prediction will be performed only at the time and space of calibration data

min_max_Q_event_1 <- CalData %>%
    filter(
        Q != -9999
    ) %>%
    group_by(x) %>%
    summarise(
        min = min(t),
        max = max(t)
    )

# A grid by Typology !
info_events_reaches <-
    list(
        # 1st event: Seine 2015
        Seine = list(
            # Pred WSE
            pred_WSE = list(
                type = "ZdX",
                INFO = list(
                    # First prediction: WSE (2015-09-29 07:35:00)
                    pred_1 = data.frame(
                        event = c(1),
                        reach = c(1, 2),
                        xmin = c(0, 5064),
                        xmax = c(5064, 137983),
                        tmin = c(109500),
                        tmax = c(109500),
                        nb_discretization = c(10, 200),
                        target_datetime = as.POSIXct(
                            "2015-09-29 07:35:00",
                            format = "%Y-%m-%d %H:%M:%S",
                            tz = "UTC"
                        ),
                        id_campaign = "campaign_WSE_1"
                    ),
                    # Second prediction: WSE (2015-09-30 00:15:00)
                    pred_2 = data.frame(
                        event = c(1),
                        reach = c(1, 2),
                        xmin = c(0, 5064),
                        xmax = c(5064, 137983),
                        tmin = c(169500),
                        tmax = c(169500),
                        nb_discretization = c(10, 200),
                        target_datetime = as.POSIXct(
                            "2015-09-30 00:15:00",
                            format = "%Y-%m-%d %H:%M:%S",
                            tz = "UTC"
                        ),
                        id_campaign = "campaign_WSE_2"
                    ),
                    # Third prediction: WSE (2015-09-30 09:55:00)
                    pred_3 = data.frame(
                        event = c(1),
                        reach = c(1, 2),
                        xmin = c(0, 5064),
                        xmax = c(5064, 137983),
                        tmin = c(204300),
                        tmax = c(204300),
                        nb_discretization = c(10, 200),
                        target_datetime = as.POSIXct(
                            "2015-09-30 09:55:00",
                            format = "%Y-%m-%d %H:%M:%S",
                            tz = "UTC"
                        ),
                        id_campaign = "campaign_WSE_3"
                    ),
                    # Fourth prediction: WSE (2015-09-30 14:10:00)
                    pred_4 = data.frame(
                        event = c(1),
                        reach = c(1, 2),
                        xmin = c(0, 5064),
                        xmax = c(5064, 137983),
                        tmin = c(219600),
                        tmax = c(219600),
                        nb_discretization = c(10, 200),
                        target_datetime = as.POSIXct(
                            "2015-09-30 14:10:00",
                            format = "%Y-%m-%d %H:%M:%S",
                            tz = "UTC"
                        ),
                        id_campaign = "campaign_WSE_4"
                    )
                )
            ),
            # Pred Q
            pred_Q = list(
                type = "QdT",
                INFO = list(
                    # Prediction: Q (Rouen)
                    pred_1 = data.frame(
                        event = c(1),
                        reach = c(2),
                        xmin = c(41568),
                        xmax = c(41568),
                        tmin = floor(
                            min_max_Q_event_1 %>% filter(x == 41568) %>% pull(min) - (60 * 5)
                        ),
                        tmax = ceiling(
                            min_max_Q_event_1 %>% filter(x == 41568) %>% pull(max) + (60 * 5)
                        ),
                        nb_discretization = c(200),
                        target_datetime = NA,
                        id_campaign = "campaign_Q_3"
                    ),
                    # Prediction: Q (Heurteauville)
                    pred_2 = data.frame(
                        event = c(1),
                        reach = c(2),
                        xmin = c(97091),
                        xmax = c(97091),
                        tmin = floor(
                            min_max_Q_event_1 %>% filter(x == 97091) %>% pull(min) - (60 * 5)
                        ),
                        tmax = ceiling(
                            min_max_Q_event_1 %>% filter(x == 97091) %>% pull(max) + (60 * 5)
                        ),
                        nb_discretization = c(200),
                        target_datetime = NA,
                        id_campaign = "campaign_Q_2"
                    ),
                    # Prediction: Q (Aizier)
                    pred_3 = data.frame(
                        event = c(1),
                        reach = c(2),
                        xmin = c(123562),
                        xmax = c(123562),
                        tmin = floor(
                            min_max_Q_event_1 %>% filter(x == 123562) %>% pull(min) - (60 * 5)
                        ),
                        tmax = ceiling(
                            min_max_Q_event_1 %>% filter(x == 123562) %>% pull(max) + (60 * 5)
                        ),
                        nb_discretization = c(200),
                        target_datetime = NA,
                        id_campaign = "campaign_Q_1"
                    )
                )
            )
        )
    )

X_pred <- grid_user_lastest(
    info_events_reaches = info_events_reaches,
    date_ref = date_ref,
    X_Cal = X_plot
)

for (id_cal_case in 1:length(all_cal_case)) {
    # Load experiment
    paths <- load_experiment(
        file_main_path = file_main_path,
        cal_case = all_cal_case[[id_cal_case]],
        path_experiment = path_experiment,
        all_events = all_events
    )

    # Load data and model used during calibration
    load(file.path(paths$path_RData, "BaM_objects.RData"))

    # Run prediction
    return_prediction <- prediction_MAGE(
        cal_case = all_cal_case[[id_cal_case]],
        paths = paths,
        BaM_data = data,
        prediction_file = c("ParamU", "Maxpost", "TotalU"),
        names_file_prediction = c("WSE", "Q", "V", "Kmin", "Kflood"),
        do_prediction = do_prediction,
        X_pred = X_pred,
        mod = list_mod_polynomials[[id_cal_case]],
        remant_error_list = remant_error_list,
        mcmcOptions = mcmcOptions,
        mcmcCooking = mcmcCooking,
        mcmcSummary = mcmcSummary,
        nsim_prior = 500
    )

    cf_file <- return_prediction$cf_file

    script_path_pred <- file.path(paths$path_BaM_folder, "run_pred_BaM.sh")
    # write the bash script
    lines <- c(
        "#!/bin/bash", # start with bash header
        "exec > /dev/null 2>&1", # From now on, send all output (stdout and stderr) to nowhere.
        ""
    )

    for (cf in cf_file) {
        # add a line for each BaM run
        lines <- c(
            lines,
            paste(
                "nohup", # run in background
                shQuote(file.path(RBaM::getPathToBaM(), "BaM")),
                "-cf",
                shQuote(cf),
                ">",
                shQuote(paste0(cf, ".log")),
                "2>&1",
                "&"
            )
        )
    }
    writeLines(lines, script_path_pred)

    Sys.chmod(script_path_pred, "0755")
    # Run outside of Vscodium. To kill a job: pkill -f '/Git/BaM/makefile/BaM'
    # See if the runs are in parallel: pgrep -af BaM

    if (do_prediction) {
        system2(
            "bash",
            args = script_path_pred,
            wait = FALSE
        )
    }
}


##########################################
### Plots
##########################################

for (id_cal_case in 1:length(all_cal_case)) {
    # Load experiment
    paths <- load_experiment(
        file_main_path = file_main_path,
        cal_case = all_cal_case[[id_cal_case]],
        path_experiment = path_experiment,
        all_events = all_events
    )

    results_postprocess <- postprocess_prediction(
        paths = paths,
        X_input = X_plot,
        date_ref = date_ref,
        Y_observations = Y,
        Yu_observations = Yu,
        conf_level = 0.95,
        summary_SU_Kmin = list_summary_SU_Kmin[[id_cal_case]],
        summary_SU_Kflood = list_summary_SU_Kflood[[id_cal_case]],
        grid = X_pred,
        Input_Typology = Input_Typology,
        suffix_patterns = c("_WSE", "_Q", "_V", "_Kmin", "_Kflood"),
        desired_order = c("Total", "Parametric", "Maxpost", "Observations")
    )
}
