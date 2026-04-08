rm(list = ls())
graphics.off()

# Set directory (root of the repository)
dir_workspace <- "/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git"
setwd(dir_workspace)

load("data/processed_data/Rhone/all_observations_AIN.RData")
load("data/processed_data/Rhone/all_observations_Rhone.RData")

function_list <- list.files("R", full.names = TRUE)
for (i in function_list) {
    source(i)
}

key_info_event_extraction <- read.table("data/processed_data/Rhone/Boundary_conditions/Specific_boundary_conditions/Dates_start_end_extraction.csv", header = TRUE, sep = ",")


sd_WSE_fixed <- 0.025 # 2.5 cm
MAGE_main_folder <- "/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/scripts/Case_studies/Rhone/Real_Condition_model/model_mage"

Input_Typology <- list(
    Rhone = c(1, 2, 3),
    Ain = c(8, 4, 5),
    Cassier = c(6, 7)
)
Input_Model_Reach <- data.frame(
    reach = c(1, 2, 3, 4, 5, 6, 7, 8),
    KP_start = c(55900, 36250, 34500, 37491, 41211, 0, 20, 22333),
    KP_end = c(36250, 34500, 26750, 41211, 41461, 20, 40, 37491)
)

all_events <- c(
    "AIN_90",
    "RHONE_525",
    "RHONE_750"
)
mage_projet_name <- "Rhone_PCH_Ain"

#####################################################################
# Event 1: WSE Q(AIN) = 90 m3/s
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
        time = as.numeric(difftime(Date_time, init_model_date, units = c("secs"))),
        Yu_WSE = sd_WSE_fixed,
        event = 1,
        name_event = "AIN_90"
    )

check_simulation_time(
    MAGE_main_folder = MAGE_main_folder,
    mage_projet_name = mage_projet_name,
    Observations = WSE_Ain_90_raw,
    event = all_events[1] # Manual modification to check
)

#################################
set.seed(2026) #  # for reproducibility
# 1. Calibration set
WSE_Ain_90_Cal <- WSE_Ain_90_raw %>%
    group_by(id_reach_CAL) %>%
    mutate(set = "calibration") %>%
    arrange(KP + id_reach_CAL) %>% # Because order is upstream to downstream and KP is increasing
    mutate(Reach_groupped_Cal = NA_character_)

WSE_obs_1 <- assign_calibration_and_validation_data(
    Input_Typology = Input_Typology,
    Input_Model_Reach = Input_Model_Reach,
    CalData = WSE_Ain_90_Cal,
    All_observations = WSE_Ain_90_raw
)

CalData_obs_1 <- WSE_obs_1 %>% filter(set == "calibration")

#####################################################################
# Event 2: WSE Q(Rhone) = 525 m3/s
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
        time = as.numeric(difftime(Date_time, init_model_date, units = c("secs"))),
        Yu_WSE = sd_WSE_fixed,
        event = 2,
        name_event = "RHONE_525"
    )

check_simulation_time(
    MAGE_main_folder = MAGE_main_folder,
    mage_projet_name = mage_projet_name,
    Observations = WSE_Rhone_525_raw,
    event = all_events[2] # Manual modification to check
)

#################################
set.seed(2026) #  # for reproducibility
# 1. Calibration set
WSE_Rhone_525_Cal <- WSE_Rhone_525_raw %>%
    group_by(id_reach_CAL) %>%
    mutate(set = "calibration") %>%
    arrange(KP + id_reach_CAL) %>% # Because order is upstream to downstream and KP is increasing
    mutate(Reach_groupped_Cal = NA_character_)

WSE_obs_2 <- assign_calibration_and_validation_data(
    Input_Typology = Input_Typology,
    Input_Model_Reach = Input_Model_Reach,
    CalData = WSE_Rhone_525_Cal,
    All_observations = WSE_Rhone_525_raw
)

CalData_obs_2 <- WSE_obs_2 %>% filter(set == "calibration")


#####################################################################
# Event 3: WSE Q(Rhone) = 750 m3/s
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
        time = as.numeric(difftime(Date_time, init_model_date, units = c("secs"))),
        Yu_WSE = sd_WSE_fixed,
        event = 3,
        name_event = "RHONE_750"
    )

check_simulation_time(
    MAGE_main_folder = MAGE_main_folder,
    mage_projet_name = mage_projet_name,
    Observations = WSE_Rhone_750_raw,
    event = all_events[3] # Manual modification to check
)

#################################
set.seed(2026) #  # for reproducibility
# 1. Calibration set (random 70%)
WSE_Rhone_750_Cal <- WSE_Rhone_750_raw %>%
    group_by(id_reach_CAL) %>%
    mutate(set = "calibration") %>%
    arrange(KP + id_reach_CAL) %>% # Because order is upstream to downstream and KP is increasing
    mutate(Reach_groupped_Cal = NA_character_)

WSE_obs_3 <- assign_calibration_and_validation_data(
    Input_Typology = Input_Typology,
    Input_Model_Reach = Input_Model_Reach,
    CalData = WSE_Rhone_750_Cal,
    All_observations = WSE_Rhone_750_raw
)

CalData_obs_3 <- WSE_obs_3 %>% filter(set == "calibration")

# All calibration data
observed_data <- rbind(
    CalData_obs_1,
    CalData_obs_2,
    CalData_obs_3
)

##################################
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
                t = WSE_Ain_90_Cal$time
            ),
        RHONE_525 =
            data.frame(
                event = WSE_Rhone_525_Cal$event,
                reach = WSE_Rhone_525_Cal$id_reach_CAL,
                x = WSE_Rhone_525_Cal$KP,
                t = WSE_Rhone_525_Cal$time
            ),
        RHONE_750 =
            data.frame(
                event = WSE_Rhone_750_Cal$event,
                reach = WSE_Rhone_750_Cal$id_reach_CAL,
                x = WSE_Rhone_750_Cal$KP,
                t = WSE_Rhone_750_Cal$time
            )
    )

save(observed_data, file = "data/processed_data/Rhone/Ain_90_RH_525_750/observed_data.RData")
