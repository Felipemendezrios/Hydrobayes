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
#############################################
# End load libraries
#############################################


############################################
# Module 1 : set directory paths
############################################

# Name of the experiment. All scenarios will be used the same calibration data
Experiment_id <- c(
    "Syn_case_rec_5_3_with_Ks_obs_limites_all_obs_flat_prior"
    # "Syn_case_with_manual_Ks_obs_limites_all_obs_flat_prior"
    # "Syn_case_rec_5_3_with_Ks_obs_limites_all_obs_flat_prior"
    # "Syn_case_without_ks_measures_all_obs_flat_prior"
)


# Experiments input data to be used during calibration setting
all_cal_case <- c(
    "Kmin_n_0_Kflood_n_0_flat.r",
    "Kmin_n_1_Kflood_n_0_flat.r",
    "Kmin_n_2_Kflood_n_0_flat.r",
    "Kmin_n_3_Kflood_n_0_flat.r",
    "Kmin_n_4_Kflood_n_0_flat.r"
)

# Folder related to the observations (careful with the order!)
all_events <- c(
    "5_1",
    "3_2"
)

command_line_MAGE <- ""
file_main_path <- "/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/scripts/Case_studies/Synthetic_case/Simplified_ks_rectangular_MC/Calibration_experiments"
MAGE_main_folder <- "/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/scripts/Case_studies/Synthetic_case/Simplified_ks_rectangular_MC/model_mage"
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
load("data/processed_data/Synthetic_case/WSE_synthetic_simplified_rectangle.RData")

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

# !!!!!!!!!!!!
# Ensure that order should be in coherence with mage model set up!!!!!!!!!!!!
# !!!!!!!!!!!!

## just to know the order of calibration data
Cal_measures <-
    list(
        X5_MR_1_TR =
            data.frame(
                event = CalData_event_1$event,
                reach = CalData_event_1$id_reach_CAL,
                x = CalData_event_1$KP,
                t = CalData_event_1$time
            ),
        X3_MR_2_TR =
            data.frame(
                event = CalData_event_2$event,
                reach = CalData_event_2$id_reach_CAL,
                x = CalData_event_2$KP,
                t = CalData_event_2$time
            )
    )


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

if (Experiment_id == "Syn_case_rec_5_3_with_Ks_obs_limites_all_obs_flat_prior") {
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
} else if (Experiment_id == "Syn_case_with_manual_Ks_obs_limites_all_obs_flat_prior") {
    Y$Kmin <- -9999
    Yu$Yu_Kmin <- -9999

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
        k = c(36, 34.5, 32.5)
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
        k = c(32, 31.5, 31.2, 27)
    )

    Y$Kmin[sliced_ks$idx] <- sliced_ks$k

    sliced_all_reaches <- rbind(sliced, sliced_reach_2)
    Yu$Yu_Kmin[sliced_all_reaches$idx] <- 0.3
}

CalData <- cbind(X, Y, Yu)

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
path_experiment <- file.path(file_main_path, Experiment_id)

if (!dir.exists(path_experiment)) {
    dir.create(path_experiment)
}

for (name in names(plots_CalData)) {
    plot <- plots_CalData[[name]]
    if (!is.null(plot)) { # ne sauvegarde que si ce n'est pas NULL
        ggsave(file.path(path_experiment, paste0(name, ".png")),
            plot = plot,
            width = 30,
            height = 20,
            units = "cm"
        )
    }
}
############################################
# End module 2: calibration data
############################################


############################################
# Module 4: calibration setting
############################################
