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
    "1_WSE_AIN_90_2_WSE_RHONE_525_750"
)


# Experiments input data to be used during calibration setting
all_cal_case <- c(
    "Kmin_Rh_SU1_SU2_n0_Ain_SU1_n0.r",
    "Kmin_Rh_SU1_n4_SU2_n2_Ain_SU1_n4.r",
    "Kmin_Rh_SU1_n4_SU2_n0_Ain_SU1_n5.r"
)

# Folder related to the observations (careful with the order!)

all_events <- c(
    "AIN_90",
    "RHONE_525",
    "RHONE_750"
)

command_line_MAGE <- "-fp=2 -LC=0 -eps=5"
file_main_path <- file.path(dir_workspace, "scripts/Case_studies/Rhone/Real_Condition_model/Calibration_experiments")
MAGE_main_folder <- "/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/scripts/Case_studies/Rhone/Real_Condition_model/model_mage"
mage_projet_name <- "Rhone_PCH_Ain"

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
    reach = c(1, 2, 3, 4, 5, 6, 7, 8),
    KP_start = c(55900, 36250, 34500, 37491, 41211, 0, 20, 22333),
    KP_end = c(36250, 34500, 26750, 41211, 41461, 20, 40, 37491)
)

# First layer: Typology
# ###########################################
# All reaches in the model (Model_Reach) defined within the same coordinate reference. Useful to separate main reach to tributary or even a diversion channel

# Must be careful with the order of the reaches, it must be given upstream to downstream
Input_Typology <- list(
    Rhone = c(1, 2, 3),
    Ain = c(8, 4, 5),
    Cassier = c(6, 7)
)

############################################
# End module 2 : hydraulic model (HM) environment
############################################


############################################
# Module 3: calibration data
############################################

# Processed data
if (Experiment_id %in% c("1_WSE_AIN_90_2_WSE_RHONE_525_750")) {
    load("data/processed_data/Rhone/Ain_90_RH_525_750/observed_data.RData")
} else {
    stop("Experiment id in not correct")
}


######################################## UNTIL HERE

X <- observed_data[, c(
    "event",
    "id_reach_CAL",
    "KP",
    "time"
)] %>%
    rename(
        "reach" = "id_reach_CAL",
        "x" = "KP",
        "t" = "time"
    )


results_CalData <- constructor_CalData(
    observed_data = observed_data,
    do_manual_uncertainty = FALSE,
    sd_WSE_fixed = 0.05, # in meters
    sd_Q_fixed = 8, # in %
    sd_V_fixed = 3, # in %
    sd_Kmin_fixed = 5, # in m1/3/s
    sd_Kflood_fixed = 5 # in m1/3/s
)

Y <- results_CalData$Y
Yu <- results_CalData$Yu

CalData <- cbind(X, Y, Yu)


##########################################
# Specificities : add pseudo obs of friction
# obs always has a gaussian distribution
Kmin_prior_u_Yu <- c(30, 10)
Kflood_prior_u_Yu <- c(25, 8)


#############################
# Add x,t coordinates to spatially distributed friction observations
#############################

# Specific of each case: here, i do not want assign pseudo obs at the reservoir
nodes_to_add_pseudo_obs <- unlist(c(Input_Typology$Rhone, Input_Typology$Ain))

# Distance between two consecutive pseudo obs
step_size <- 150

K_pseudo_obs <- data.frame(Input_Model_Reach %>%
    filter(reach %in% nodes_to_add_pseudo_obs) %>%
    rowwise() %>%
    mutate(KP = list(
        seq(KP_start,
            KP_end,
            length.out = ceiling(abs(KP_end - KP_start) / step_size) + 1
        )
    )) %>%
    unnest(KP) %>%
    select(reach, KP) %>%
    rename(x = KP) %>%
    mutate(
        event = 1,
        t = CalData$t[1]
    ) %>%
    select(event, reach, x, t) %>%
    mutate(
        WSE = -9999,
        Q = -9999,
        V = -9999,
        Kmin = Kmin_prior_u_Yu[1],
        Kflood = Kflood_prior_u_Yu[1],
        Yu_WSE = -9999,
        Yu_Q = -9999,
        Yu_V = -9999,
        Yu_Kmin = Kmin_prior_u_Yu[2],
        Yu_Kflood = Kflood_prior_u_Yu[2]
    ))


CalData <- rbind(CalData, K_pseudo_obs)

path_experiment <- file.path(file_main_path, Experiment_id)

if (do_plot_calibration) {
    CalData_ordered <- CalData %>%
        mutate(reach = factor(reach,
            levels = c(1, 2, 3, 8, 4, 5)
        ))

    plots_CalData <- plot_CalData(
        CalData = CalData_ordered,
        scales_free = "free",
        wrap = TRUE
    )
    # Customize the plot
    plots_CalData$plot_WSE <- plots_CalData$plot_WSE +
        facet_wrap(
            ~ event + reach,
            labeller = labeller(
                event = as_labeller(
                    c(
                        "1" = "Event~1:~Q(Ain):~90~m^3/s",
                        "2" = "Event~2:~Q(Rhone):~525~m^3/s",
                        "3" = "Event~3:~Q(Rhone):~750~m^3/s"
                    ),
                    label_parsed
                ),
                reach = as_labeller(c(
                    "1" = "Reach~1:~Lagnieu~-~Bourbe",
                    "2" = "Reach~2:~Bourbe~-~Anthon",
                    "3" = "Reach~3:~Anthon~-~Jons",
                    "4" = "Reach~4:~Port~Galland~-~reservoir",
                    "5" = "Reach~5:~reservoir~-~Anthon",
                    "8" = "Reach~8:~Pont~Chazey~-~Port~Galland"
                ), label_parsed)
            ),
            scales = "free",
            ncol = 3
        )

    plots_CalData$plot_Kmin <- plots_CalData$plot_Kmin +
        facet_wrap(
            ~reach,
            labeller = labeller(
                reach = as_labeller(c(
                    "1" = "Reach~1:~Lagnieu~-~Bourbe",
                    "2" = "Reach~2:~Bourbe~-~Anthon",
                    "3" = "Reach~3:~Anthon~-~Jons",
                    "4" = "Reach~4:~Port~Galland~-~reservoir",
                    "5" = "Reach~5:~reservoir~-~Anthon",
                    "8" = "Reach~8:~Pont~Chazey~-~Port~Galland"
                ), label_parsed)
            ),
            scales = "free",
            ncol = 3
        )

    plots_CalData$plot_Kflood <- plots_CalData$plot_Kflood +
        facet_wrap(
            ~reach,
            labeller = labeller(
                reach = as_labeller(c(
                    "1" = "Reach~1:~Lagnieu~-~Bourbe",
                    "2" = "Reach~2:~Bourbe~-~Anthon",
                    "3" = "Reach~3:~Anthon~-~Jons",
                    "4" = "Reach~4:~Port~Galland~-~reservoir",
                    "5" = "Reach~5:~reservoir~-~Anthon",
                    "8" = "Reach~8:~Pont~Chazey~-~Port~Galland"
                ), label_parsed)
            ),
            scales = "free",
            ncol = 3
        )

    obs_adapted <- observed_data %>%
        rename("reach" = "id_reach_CAL")

    plot_WSE_Thalweg <- plots_CalData$plot_WSE +
        geom_line(data = obs_adapted, aes(x = KP, y = Z_thalweg), color = "black")

    plots_CalData$plot_WSE_Thalweg <- plot_WSE_Thalweg

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
            init = 0.0005,
            prior.dist = "FlatPrior"
        ))
    ),
    # Q
    remnantErrorModel_default(name = "Config_RemnantSigma2.txt"),
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
    total_points_discretization = 100
)

############################################
# Module 6: Calibration
############################################
list_mod_polynomials <- list_Z_MatrixKmin <- list_Z_MatrixKflood <- list_Kmin_prior <- list_Kflood_prior <- list_Kmin_SU <- list_Kflood_SU <- list_summary_SU_Kflood <- list_summary_SU_Kmin <- list()

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
        mcmcCooking = RBaM::mcmcCooking(burn = 0, nSlim = 1),
        mcmcOptions = RBaM::mcmcOptions(nAdapt = 20, nCycles = 10),
        mcmcSummary = RBaM::mcmcSummary(xtendedMCMC.fname = "Results_xtendedMCMC.txt"),
        remant_error_list = remant_error_list
    )
    if (do_calibration) {
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
        # Run outside of Vscodium. To kill a job : pkill -f BaM
        system2(
            "nohup",
            args = c("bash", script_path),
            stdout = FALSE,
            stderr = FALSE,
            wait = FALSE
        )
    }
    list_Z_MatrixKmin[[id_cal_case]] <- results_estimation$Z_MatrixKmin
    list_Z_MatrixKflood[[id_cal_case]] <- results_estimation$Z_MatrixKflood
    list_Kmin_prior[[id_cal_case]] <- results_estimation$Kmin_prior
    list_Kflood_prior[[id_cal_case]] <- results_estimation$Kflood_prior
    list_Kmin_SU[[id_cal_case]] <- results_estimation$Kmin_SU
    list_Kflood_SU[[id_cal_case]] <- results_estimation$Kflood_SU
    list_mod_polynomials[[id_cal_case]] <- results_estimation$mod
    list_summary_SU_Kmin[[id_cal_case]] <- results_estimation$summary_SU_Kmin
    list_summary_SU_Kflood[[id_cal_case]] <- results_estimation$summary_SU_Kflood
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


################################
# POSTPROCESS CALIBRATION WORKFLOW
################################
final_calibration <- TRUE

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
        X_input = X,
        Y_observations = Y,
        Yu_observations = Yu,
        type = "dx",
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

    if (do_plot_calibration) {
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
}
