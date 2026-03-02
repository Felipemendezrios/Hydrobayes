rm(list = ls())
graphics.off()

# Set directory (root of the repository)
dir_workspace <- here::here()
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
do_plot_calibration <- TRUE
do_prediction <- TRUE


# Setting jump standard deviation for MCMC sampling
jump_MCMC_theta_param_user <- 8
jump_MCMC_error_model_user <- 0.001
threshold_jump_MCMC_error_model <- 0.5


############################################
# Module 1 : set directory paths
############################################


# Name of the experiment. All scenarios will be used the same calibration data

# 2_WSE_KMIN_high: all CalData have high unc, KMIN is overall distributed
# 2_WSE_KMIN_low: all CalData have low unc, KMIN is in some positions
# 2_WSE_high_KMIN_low: WSE has high unc, KMIN is in some positions
# 2_WSE_low_KMIN_high: WSE has low unc, KMin is overall distributed

Experiment_id <- c(
    # "2_WSE_KMIN_high" # ok
    # "2_WSE_low_KMIN_high" #ok
    # "2_WSE_KMIN_low" # ok
    # "2_WSE_high_KMIN_low" # ok
    "2_WSE_high_KMIN_low_distributed" # only tested in piecewise function
    # "2_WSE_KMIN_low_distributed" # only tested in piecewise function
)

# Experiments input data to be used during calibration setting
all_cal_case <- c(
    # "Kmin_n_0_Kflood_n_0.r",
    # "Kmin_n_1_Kflood_n_0.r",
    # "Kmin_n_2_Kflood_n_0.r",
    # "Kmin_n_3_Kflood_n_0.r",
    # "Kmin_n_4_Kflood_n_0.r",
    # "Kmin_n_5_Kflood_n_0.r",
    # "Kmin_n_6_Kflood_n_0.r",
    "Piecewise_Kmin_n_0_Kflood_n_0.r"
)



# Folder related to the observations (careful with the order!)
all_events <- c(
    "5_1",
    "3_2"
)

command_line_MAGE <- ""
file_main_path <- file.path(dir_workspace, "scripts/Case_studies/Synthetic_case/Rectangular_MC/Calibration_experiments")
MAGE_main_folder <- file.path(dir_workspace, "scripts/Case_studies/Synthetic_case/Rectangular_MC/model_mage")
mage_projet_name <- "synt_rect"


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
    reach = c(1, 2, 3),
    KP_start = c(0, 0, 18),
    KP_end = c(18, 20, 25)
)

# First layer: Typology
# ###########################################
# All reaches in the model (Model_Reach) defined within the same coordinate reference. Useful to separate main channel to tributary or even a diversion channel

# Must be careful with the order of the reaches, it must be given upstream to downstream
Input_Typology <- list(
    Main_Reach = c(1, 3),
    Tributary = c(2)
)

############################################
# End module 2 : hydraulic model (HM) environment
############################################


############################################
# Module 3: calibration data
############################################
# Processed data
if (Experiment_id %in% c("2_WSE_KMIN_high", "2_WSE_high_KMIN_low")) {
    load("data/processed_data/Synthetic_case/Rectangular_MC/High_uncertainty/WSE_synthetic_rectangular_MC.RData")
} else if (Experiment_id %in% c("2_WSE_KMIN_low", "2_WSE_low_KMIN_high", "2_WSE_high_KMIN_low_distributed", "2_WSE_KMIN_low_distributed")) {
    load("data/processed_data/Synthetic_case/Rectangular_MC/Low_uncertainty/WSE_synthetic_rectangular_MC.RData")
} else {
    stop("Experiment id in not correct")
}

# Observations data input:
# Measurements for calibration by event!

# Read key information of event
#####################################################################
# Event 1: WSE in the main reach 5 m3/s
#####################################################################


WSE_synthetic_simplified_raw_event_1 <- WSE_synthetic_simplified %>%
    filter(id_case == "5_1", id_reach_CAL %in% Input_Typology[[1]]) %>%
    mutate(
        id_reach = case_when(
            id_reach_CAL == 1 ~ "Main reach upstream",
            id_reach_CAL == 2 ~ "Tributary",
            id_reach_CAL == 3 ~ "Main reach downstream"
        ),
        event = 1,
        name_event = "5_MR_1_TR"
    )

check_simulation_time(
    MAGE_main_folder = MAGE_main_folder,
    mage_projet_name = mage_projet_name,
    Observations = WSE_synthetic_simplified_raw_event_1,
    event = all_events[1]
)
#################################
set.seed(2026) #  # for reproducibility
# 1. Calibration set
WSE_synthetic_simplified_Cal_event_1 <- WSE_synthetic_simplified_raw_event_1 %>%
    group_by(id_reach_CAL) %>%
    # slice_sample(prop = 0.8) %>%
    ungroup() %>%
    mutate(set = "calibration") %>%
    arrange(KP + id_reach_CAL) %>% # Because order is upstream to downstream and KP is increasing
    mutate(Reach_groupped_Cal = NA_character_)


WSE_Event_1 <- assign_calibration_and_validation_data(
    Input_Typology = Input_Typology,
    Input_Model_Reach = Input_Model_Reach,
    CalData = WSE_synthetic_simplified_Cal_event_1,
    All_observations = WSE_synthetic_simplified_raw_event_1
)

CalData_event_1 <- WSE_Event_1 %>% filter(set == "calibration")

ggplot(WSE_Event_1, aes(x = KP, y = WSE, col = set)) +
    geom_point() +
    geom_errorbar(aes(
        ymin = WSE - qnorm(0.975) * Yu_WSE,
        ymax = WSE + qnorm(0.975) * Yu_WSE
    )) +
    facet_wrap(~id_reach_CAL)

ggplot(WSE_Event_1, aes(x = KP, y = WSE, col = factor(id_reach_CAL))) +
    geom_point() +
    geom_errorbar(aes(
        ymin = WSE - qnorm(0.975) * Yu_WSE,
        ymax = WSE + qnorm(0.975) * Yu_WSE
    ))



#####################################################################
# Event 2: WSE in the tributary
#####################################################################

WSE_synthetic_simplified_raw_event_2 <- WSE_synthetic_simplified %>%
    filter(id_case == "3_2", id_reach_CAL %in% Input_Typology[[2]]) %>%
    mutate(
        id_reach = case_when(
            id_reach_CAL == 1 ~ "Main reach upstream",
            id_reach_CAL == 2 ~ "Tributary",
            id_reach_CAL == 3 ~ "Main reach downstream"
        ),
        event = 2,
        name_event = "3_MR_2_TR"
    )
#################################

check_simulation_time(
    MAGE_main_folder = MAGE_main_folder,
    mage_projet_name = mage_projet_name,
    Observations = WSE_synthetic_simplified_raw_event_2,
    event = all_events[2]
)

#################################
set.seed(2026) #  # for reproducibility
# 1. Calibration set
WSE_synthetic_simplified_Cal_event_2 <- WSE_synthetic_simplified_raw_event_2 %>%
    group_by(id_reach_CAL) %>%
    # slice_sample(prop = 0.7) %>%
    ungroup() %>%
    mutate(set = "calibration") %>%
    arrange(KP + id_reach_CAL) %>% # Because order is upstream to downstream and KP is increasing
    mutate(Reach_groupped_Cal = NA_character_)


WSE_Event_2 <- assign_calibration_and_validation_data(
    Input_Typology = Input_Typology,
    Input_Model_Reach = Input_Model_Reach,
    CalData = WSE_synthetic_simplified_Cal_event_2,
    All_observations = WSE_synthetic_simplified_raw_event_2
)

CalData_event_2 <- WSE_Event_2 %>% filter(set == "calibration")

ggplot(WSE_Event_2, aes(x = KP, y = WSE, col = set)) +
    geom_point() +
    geom_errorbar(aes(
        ymin = WSE - qnorm(0.975) * Yu_WSE,
        ymax = WSE + qnorm(0.975) * Yu_WSE
    )) +
    facet_wrap(~id_reach_CAL)

ggplot(CalData_event_2, aes(x = KP, y = WSE, col = factor(set))) +
    geom_point() +
    geom_errorbar(aes(
        ymin = WSE - qnorm(0.975) * Yu_WSE,
        ymax = WSE + qnorm(0.975) * Yu_WSE
    ))

# All available data
WSE_data <- rbind(
    WSE_Event_1,
    WSE_Event_2
)
# All calibration data
observed_data <- rbind(
    CalData_event_1,
    CalData_event_2
)

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

##########################################
# Specificities : add pseudo obs of friction

if (Experiment_id %in% c("2_WSE_KMIN_high", "2_WSE_low_KMIN_high")) {
    Y$Kmin <- ifelse(X$reach == 1 | X$reach == 3, 34, 28)

    Yu_Kmin <- ifelse(X$reach == 1 | X$reach == 3, 8, 8)

    toto <- cbind(X, Yu)
    toto$Yu_Kmin <- Yu_Kmin

    # Step 1: assign group based on Input_XR
    toto_grouped <- toto %>%
        rowwise() %>%
        mutate(group = case_when(
            reach %in% Input_Typology$Main_Reach ~ "MR",
            reach %in% Input_Typology$Tributary ~ "TR",
            TRUE ~ NA_character_
        )) %>%
        ungroup()
    # Step 2: summarize min/max x by event and group
    results <- toto_grouped %>%
        group_by(event, group) %>%
        summarise(
            min_x = min(x),
            max_x = max(x),
            .groups = "drop"
        )
    # This gives min/max per group

    # Step 3: Update Y_Kmin based on the group (example: set MR=6, TR=8)
    # You can map groups to new values like this:
    group_values <- c(MR = 4, TR = 5)

    toto_updated <- toto_grouped %>%
        group_by(event, group) %>%
        mutate(
            min_x = min(x),
            max_x = max(x),
            Yu_Kmin = if_else(
                x == min_x | x == max_x,
                group_values[group],
                Yu_Kmin
            )
        ) %>%
        ungroup() %>%
        select(-min_x, -max_x)

    Yu$Yu_Kmin <- toto_updated$Yu_Kmin
} else if (Experiment_id %in% c("2_WSE_KMIN_low", "2_WSE_high_KMIN_low", "2_WSE_high_KMIN_low_distributed", "2_WSE_KMIN_low_distributed")) {
    # Add boundaries values of ks for the rest
    if (Experiment_id %in% c("2_WSE_high_KMIN_low_distributed", "2_WSE_KMIN_low_distributed")) {
        Y$Kmin <- ifelse(X$reach == 1 | X$reach == 3, 34, 28)
        Yu$Yu_Kmin <- ifelse(X$reach == 1 | X$reach == 3, 8, 8)
    } else {
        Y$Kmin <- -9999
        Yu$Yu_Kmin <- -9999
    }

    all_obs_idx <- which(X$reach == 1 | X$reach == 3)

    all_slice_info <- data.frame(
        idx = all_obs_idx,
        x_specif = X$x[all_obs_idx]
    )

    positions_specified <- c(
        all_slice_info$idx[which.min(all_slice_info$x_specif)],
        all_slice_info$idx[which(all_slice_info$x_specif == 15)],
        all_slice_info$idx[which.max(all_slice_info$x_specif)]
    )

    sliced <- all_slice_info[which(all_slice_info$idx %in% positions_specified), ]
    sliced_ks <- data.frame(sliced,
        k = c(35, 28, 35)
    )

    Y$Kmin[sliced_ks$idx] <- sliced_ks$k

    # Next SU
    all_obs_idx <- which(X$reach == 2)

    all_slice_info <- data.frame(
        idx = all_obs_idx,
        x_specif = X$x[all_obs_idx]
    )

    positions_specified <- c(
        all_slice_info$idx[which.min(all_slice_info$x_specif)],
        all_slice_info$idx[which(all_slice_info$x_specif == 4)],
        all_slice_info$idx[which(all_slice_info$x_specif == 7)],
        all_slice_info$idx[which.max(all_slice_info$x_specif)]
    )

    sliced_reach_2 <- all_slice_info[which(all_slice_info$idx %in% positions_specified), ]

    sliced_ks <- data.frame(sliced_reach_2,
        k = c(25, 25, 27, 30)
    )

    Y$Kmin[sliced_ks$idx] <- sliced_ks$k

    sliced_all_reaches <- rbind(sliced, sliced_reach_2)
    Yu$Yu_Kmin[sliced_all_reaches$idx] <- 0.3
}

CalData <- cbind(X, Y, Yu)

path_experiment <- file.path(file_main_path, Experiment_id)

if (do_plot_calibration) {
    plots_CalData <- plot_CalData(
        CalData = CalData,
        scales_free = "free",
        wrap = TRUE
    )
    # Customize the plot
    plots_CalData$plot_WSE <- plots_CalData$plot_WSE +
        facet_wrap(
            ~event,
            labeller = as_labeller(
                c(
                    "1" = "Main~reach:~discharge~of~5~m^3/s",
                    "2" = "Tributary:~discharge~of~1~m^3/s"
                ),
                label_parsed
            ),
            scales = "free",
            ncol = 1
        )

    plots_CalData$plot_Kmin <- plots_CalData$plot_Kmin +
        facet_wrap(
            ~event,
            labeller = as_labeller(
                c(
                    "1" = "Main~reach:~discharge~of~5~m^3/s",
                    "2" = "Tributary:~discharge~of~1~m^3/s"
                ),
                label_parsed
            ),
            scales = "free",
            ncol = 1
        )

    plot_WSE_Thalweg <- plots_CalData$plot_WSE +
        geom_line(data = observed_data, aes(x = KP, y = Z_thalweg), color = "black")

    plots_CalData$plot_WSE_Thalweg <- plot_WSE_Thalweg

    if (!dir.exists(path_experiment)) {
        dir.create(path_experiment)
    }

    for (name in names(plots_CalData)) {
        plot <- plots_CalData[[name]]
        if (!is.null(plot)) {
            ggsave(file.path(path_experiment, paste0(name, ".png")),
                plot = plot,
                width = 30,
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
list_mod_polynomials <- list_Z_MatrixKmin <- list_Z_MatrixKflood <- list_Kmin_prior <- list_Kflood_prior <- list_Kmin_SU <- list_Kflood_SU <- list()

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
        mcmcCooking = RBaM::mcmcCooking(burn = 0.5, nSlim = 10),
        mcmcOptions = RBaM::mcmcOptions(nAdapt = 100, nCycles = 100),
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
}

# Plot DIC
if (do_plot_calibration) {
    plotDIC <- plot_DIC(path_experiment)
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


###############
# Theoretical values
####################
file_RUGFILE_synt_obs <- "data/processed_data/Synthetic_case/Rectangular_MC/TRUE_RUGFile.RUG"
Real_Ks_simulated <- read_fortran_data(
    file_path = file_RUGFILE_synt_obs,
    col_widths_RUGFile = c(1, 3, 6, 10, 10, 10, 10),
    skip = 1
)

Kmin_literature <- data.frame(
    x_start = Real_Ks_simulated$KP_start,
    x_end = Real_Ks_simulated$KP_end,
    mean_value = Real_Ks_simulated$Kmin,
    reaches_SU = ifelse(Real_Ks_simulated$id_reach == 1 | Real_Ks_simulated$id_reach == 3, "MR", "TR")
)

Kmin_segment_layer <- segment_layer_reference(
    K_literature = Kmin_literature,
    mean_col = "mean_value",
    min_col = NULL,
    max_col = NULL
)

####################
# End theoretical values
####################

# Add synthetic data
synthetic_case <- TRUE
if (synthetic_case) {
    real_synt_data <- WSE_synthetic_simplified %>%
        mutate(X1_obs = ifelse((id_reach_CAL == 1 | id_reach_CAL == 3) & id_case == "5_1",
            1,
            ifelse((id_reach_CAL == 2) & id_case == "3_2", 2, NA)
        )) %>%
        tidyr::drop_na(X1_obs)
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
        Kmin_prior = list_Kmin_prior[[id_cal_case]],
        Kflood_prior = list_Kflood_prior[[id_cal_case]],
        Kmin_SU = list_Kmin_SU[[id_cal_case]],
        Kflood_SU = list_Kflood_SU[[id_cal_case]],
        Z_MatrixKmin = list_Z_MatrixKmin[[id_cal_case]],
        Z_MatrixKflood = list_Z_MatrixKflood[[id_cal_case]],
        mod_polynomials = list_mod_polynomials[[id_cal_case]],
        Kmin_segment_layer = Kmin_segment_layer,
        Kflood_segment_layer = NULL,
        command_line_MAGE = command_line_MAGE,
        dir_workspace = dir_workspace
    )

    CalData_updated <- results_postprocess$CalData_updated
    residuals <- results_postprocess$residuals
    plot_Kmin_without_obs <- results_postprocess$plots_param$Kmin$plot_without_obs

    plot_Kmin_with_obs <- results_postprocess$plots_param$Kmin$plot_with_obs
    if (any(!is.na(CalData_updated[, c("Kmin")]))) {
        plot_Kmin_with_obs <- plot_Kmin_with_obs +
            geom_point(
                data = CalData_updated,
                aes(x = x, y = Kmin, col = "Pseudo \nobs"),
                alpha = 0.3
            ) +
            geom_errorbar(
                data = CalData_updated,
                aes(x = x, ymin = Kmin - 1.96 * Yu_Kmin, ymax = Kmin + 1.96 * Yu_Kmin, col = "Pseudo \nobs"),
                alpha = 0.3
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
                data = CalData_updated, aes(x = x, y = Kflood, col = "Pseudo \nobs"), alpha = 0.3
            ) +
            geom_errorbar(data = CalData_updated, aes(x = x, ymin = Kflood - 1.96 * Yu_Kflood, ymax = Kflood + 1.96 * Yu_Kflood, col = "Pseudo \nobs"), alpha = 0.3) +
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
                    aes(x = KP, y = WSE_real_obs, col = "synthetic data"), alpha = 0.7, shape = 2
                ) +
                scale_color_manual(
                    values =
                        c(
                            "sim" = "blue",
                            "obs" = "red",
                            "synthetic data" = "black"
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
info_events_reaches <- list(
    # 1st event: WSE.
    event_1 = list(
        type = "ZdX",
        SU1 = data.frame(
            event = c(1),
            reach = c(1, 3),
            xmin = c(0, 18),
            xmax = c(18, 25),
            tmin = c(10800, 10800),
            tmax = c(10800, 10800),
            nb_discretization = c(100, 100)
        )
    ),
    # 2nd event: WSE
    event_2 = list(
        type = "ZdX",
        SU2 = data.frame(
            event = c(2),
            reach = c(2),
            xmin = (0),
            xmax = c(20),
            tmin = c(10800),
            tmax = c(10800),
            nb_discretization = c(100)
        )
    )
)

X_pred <- grid_user(info_events_reaches)
X_pred_grid <- list()

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
        prediction_file = c("Prior", "ParamU", "Maxpost", "TotalU"),
        data = data,
        do_prediction = do_prediction,
        X_pred = X_pred,
        mod = list_mod_polynomials[[id_cal_case]],
        remant_error_list = remant_error_list,
        mcmcOptions = mcmcOptions,
        mcmcCooking = mcmcCooking,
        mcmcSummary = mcmcSummary,
        nsim_prior = 500
    )
    X_pred_grid[[id_cal_case]] <- return_prediction$X_pred_grid


    if (do_prediction) {
        cf_file <- return_prediction$cf_file

        script_path_pred <- file.path(paths$path_BaM_folder, "run_pred_BaM.sh")
        # write the bash script
        lines <- c("#!/bin/bash", "") # start with bash header

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
        # Run outside of Vscodium. To kill a job: pkill -f BaM
        # See if the runs are in parallel: pgrep -af BaM
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
        type = "dX",
        X_input = X,
        Y_observations = Y,
        Yu_observations = Yu,
        conf_level = 0.95,
        grid = X_pred_grid[[id_cal_case]],
        Input_Typology = Input_Typology,
        suffix_patterns = c("_WSE", "_Q", "_V", "_Kmin", "_Kflood"),
        desired_order = c("Total", "Parametric", "Maxpost", "Observations")
    )
}
