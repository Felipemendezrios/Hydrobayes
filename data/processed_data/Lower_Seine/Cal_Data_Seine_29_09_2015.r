rm(list = ls())
graphics.off()

library(dplyr)
library(tidyr)
library(ggplot2)

# Set directory (root of the repository)
dir_workspace <- "/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git"
setwd(dir_workspace)

# Table of KP

station_KP <- data.frame(
    station = c("Poses", "Elbeuf", "Oissel", "Rouen", "Petit_Couronne", "La_Bouille", "Val_des_leux", "Duclair", "Mesnil_sous_Jumieges", "Heurteauville", "Caudebec", "Vatteville", "Aizier", "Saint_Leonard", "Tancarville"),
    KP = c(0, 16823, 27439, 41568, 51186, 58878, 65055, 78482, 85578, 97091, 111227, 116862, 123562, 130623, 137983),
    id_reach_HM = c(1, rep(2, 14))
)


all_H_obs <- get(load("data/processed_data/Lower_Seine/Observations/all_dataset_H_obs.RData"))
all_Q_obs <- get(load("data/processed_data/Lower_Seine/Observations/all_dataset_Q_obs.RData"))

all_Q_obs <- all_Q_obs %>%
    rename(
        datetime_TU = time,
        id_event = source
    )

all_H_obs <- all_H_obs %>%
    rename(
        variable = value,
        station = id_station
    )


function_list <- list.files("R", full.names = TRUE)
for (i in function_list) {
    source(i)
}

################################
# Discharge observations (Q_obs)
################################

# Metadata: link between real dates and simulation period

# all_available_events: vector with all available models
all_available_events <- unique(all_Q_obs$id_event)

# id_station: list of stations with the variable to take into account during calibration for each event
# If WSE is given, NA station is necessary
id_station_var <- list(
    "Piney_2015" = data.frame(
        station = c("Rouen", "Heurteauville", "Aizier", NA),
        var = c("Q", "Q", "Q", "WSE")
    ),
    "Druine_2016" = data.frame(
        station = c("Rouen", NA),
        var = c("Q", "WSE")
    ),
    "Druine_2017" = data.frame(
        station = c("Rouen", NA),
        var = c("Q", "WSE")
    )
)

name_events <- c(
    Piney_2015 = "09_2015",
    Druine_2016 = "2016",
    Druine_2017 = "2017"
)

check_WSE <- function(id_station_var) {
    for (event in names(id_station_var)) {
        x <- id_station_var[[event]]

        WSE <- x$var == "WSE"

        if (any(WSE) && any(!is.na(x$station[WSE]))) {
            stop(
                paste0(
                    "For event '", event,
                    "', station must be NA when var = 'WSE'."
                )
            )
        }
    }
}
check_WSE(id_station_var)

# extract_date_WSE: List of dates to extract the WSE, so time fixed
extract_date_WSE <- list(
    "Piney_2015" = c(
        as.POSIXct(
            "2015-09-29 10:00:00",
            tz = "UTC",
            format = "%Y-%m-%d %H:%M:%S"
        ),
        as.POSIXct(
            "2015-09-29 20:00:00",
            tz = "UTC",
            format = "%Y-%m-%d %H:%M:%S"
        ),
        as.POSIXct(
            "2015-09-30 06:00:00",
            tz = "UTC",
            format = "%Y-%m-%d %H:%M:%S"
        ),
        as.POSIXct(
            "2015-09-30 16:00:00",
            tz = "UTC",
            format = "%Y-%m-%d %H:%M:%S"
        )
    ),
    "Druine_2016" = NULL,
    "Druine_2017" = NULL
)

# # Extract the visualize the WSE
# toto <- all_H_obs %>%
#     filter(datetime_TU %in% extract_date_WSE[[1]]) %>%
#     arrange(datetime_TU) %>%
#     mutate(
#         kp = station_KP$KP[match(station, station_KP$station)]
#     )


# ggplot(toto %>% arrange(datetime_TU), aes(y = value)) +
#     geom_point(aes(x = kp)) +
#     facet_wrap(~datetime_TU, scales = "free")

# model_start: list of dates of each event starts in the model. NULL means that experience is not simulated
model_start <- list(
    "Piney_2015" = as.POSIXct(
        "2015-09-28 01:10:00",
        tz = "UTC",
        format = "%Y-%m-%d %H:%M:%S"
    ),
    "Druine_2016" = NULL,
    "Druine_2017" = NULL
)

# model_end: list of dates of each event ends
model_end <- list(
    "Piney_2015" = as.POSIXct(
        "2015-09-30 23:25:00",
        tz = "UTC",
        format = "%Y-%m-%d %H:%M:%S"
    ),
    "Druine_2016" = NULL,
    "Druine_2017" = NULL
)

if (any(unlist(model_end) <= unlist(model_start))) {
    stop("increasing interval values are necessary")
}
for (j in seq_along(extract_date_WSE)) {
    for (i in seq_along(extract_date_WSE[[j]])) {
        if (!between(
            extract_date_WSE[[j]][[i]],
            model_start[[j]],
            model_end[[j]]
        )) {
            stop(paste0(
                "WSE observed at ",
                extract_date_WSE[[j]][[i]],
                " is outside of the interval of the model simulation"
            ))
        }
    }
}

# Date of event to extract
key_info_event_extraction <- data.frame()
for (i in seq_along(all_available_events)) {
    event <- all_available_events[i]


    start <- model_start[[event]]
    end <- model_end[[event]]
    stations_var <- id_station_var[[event]]
    dates_WSE <- extract_date_WSE[[event]]

    measured <- data.frame()
    # Model already exists
    if (!is.null(start)) {
        # Discharge measurements
        if (exists("all_Q_obs")) {
            if (any(stations_var$var %in% "Q")) {
                measured_temps <- all_Q_obs %>%
                    filter(
                        id_event == event,
                        station %in% stations_var$station
                    ) %>%
                    group_by(station) %>%
                    summarise(
                        measured_start = min(datetime_TU),
                        measured_end = max(datetime_TU),
                        .groups = "drop"
                    ) %>%
                    left_join(
                        stations_var %>%
                            filter(!is.na(station)),
                        by = "station"
                    )
                measured <- rbind(measured, measured_temps)
            }
        }
        # WSE measurements
        if (exists("all_H_obs")) {
            if (!is.null(dates_WSE)) {
                measured_temps <-
                    lapply(dates_WSE, function(date_WSE) {
                        all_H_obs %>%
                            mutate(
                                diff_sec_temp = abs(as.numeric(
                                    difftime(datetime_TU, date_WSE, units = "secs")
                                ))
                            ) %>%
                            slice_min(diff_sec_temp, n = 1, with_ties = FALSE) %>%
                            mutate(
                                var = "WSE",
                                station = NA,
                                measured_start = date_WSE,
                                measured_end = date_WSE
                            ) %>%
                            select(-diff_sec_temp)
                    }) %>%
                    bind_rows() %>%
                    arrange(datetime_TU) %>%
                    select(station, measured_start, measured_end, var)

                measured <- rbind(measured, measured_temps)
            }
        }


        key_info_event_extraction_temp <- measured %>%
            mutate(
                id_event = event,
                model_start = start,
                model_end = end,
                diff_start = as.numeric(
                    difftime(measured_start, start, units = "secs")
                ),
                diff_end = as.numeric(
                    difftime(measured_end, start, units = "secs")
                ),
                diff_model_sec = as.numeric(
                    difftime(end, start, units = "secs")
                ),
                diff_model_fmt = sprintf(
                    "%03d:%02d:%02d:%02d",
                    diff_model_sec %/% (24 * 3600),
                    (diff_model_sec %% (24 * 3600)) %/% 3600,
                    (diff_model_sec %% 3600) %/% 60,
                    diff_model_sec %% 60
                )
            ) %>%
            relocate(
                id_event,
                model_start,
                model_end,
                diff_model_sec,
                diff_model_fmt,
                var,
                .before = everything()
            )
    } else {
        key_info_event_extraction_temp <- NULL
    }

    key_info_event_extraction <- bind_rows(
        key_info_event_extraction,
        key_info_event_extraction_temp
    )
}

################################################################
# General settings: Methodology with hydraulic model (HM)
################################################################

# Mage folder
MAGE_main_folder <- "/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/processed_data/Lower_Seine/Seine_2015"
# Name of the mage file
mage_projet_name <- "Lower_Seine"

# Event to be considered during calibration
all_events_cal <- "Piney_2015"
if (!all_events_cal %in% all_available_events) stop("Event does not exist")


# Typology: link between nom of the reach and the number of the hydraulic model
Input_Typology <- list(
    Seine = c(1, 2)
)
typology_lookup <- tibble(
    id_reach_CAL = unlist(Input_Typology),
    Reach_groupped_Cal = rep(
        names(Input_Typology),
        lengths(Input_Typology)
    )
)

# Link between reach and KP in the hydraulic model
Input_Model_Reach <- data.frame(
    reach = c(1, 2),
    KP_start = c(0, 5064),
    KP_end = c(5064, 137983)
)

######################################################
# Create data to use during calibration
######################################################

#####################################################################
# Event 1: 2015-09-28 (Piney_2015)
#####################################################################

############################
# Discharge measurements
############################
# Initial date of the modeling
init_model_date <- key_info_event_extraction %>%
    filter(
        id_event %in% all_events_cal,
        var == "Q"
    ) %>%
    select(
        id_event,
        station,
        var,
        model_start
    )


Q_obs <- all_Q_obs %>%
    # Filter the event to be considered in the calibration
    filter(id_event %in% all_events_cal) %>%
    # Filter the measures to be used as calibration data to each event
    inner_join(
        bind_rows(
            id_station_var[all_events_cal],
            .id = "id_event"
        ) %>%
            filter(var == "Q"),
        by = c("id_event", "station")
    ) %>%
    group_by(station) %>%
    left_join(
        station_KP,
        by = c("station")
    ) %>%
    rename(
        id_reach_CAL = id_reach_HM
    ) %>%
    left_join(
        init_model_date,
        by = c(
            "id_event",
            "station",
            "var"
        )
    ) %>%
    mutate(
        # id_reach_CAL = case_when(
        #     id_reach == "PCH_PGA" ~ 6,
        #     id_reach == "PGA_CAIN" ~ 4,
        #     id_reach == "CAIN_Confluence" ~ 5
        # ),
        target_datetime = datetime_TU,
        time = as.numeric(
            difftime(
                datetime_TU,
                model_start,
                units = "secs"
            )
        ),
        event = match(id_event, all_events_cal),
        name_event = unname(name_events[id_event]),
        set = "calibration"
    ) %>%
    left_join(
        typology_lookup,
        by = "id_reach_CAL"
    ) %>%
    # Average duplicated measurements
    group_by(
        station,
        id_event,
        id_reach_CAL,
        KP,
        var,
        datetime_TU
    ) %>%
    mutate(
        variable = mean(variable, na.rm = TRUE),
        uncertainty = mean(uncertainty, na.rm = TRUE),
        gauging_adcp = paste(
            unique(na.omit(gauging_adcp)),
            collapse = " / "
        )
    ) %>%
    # Garder une seule ligne par groupe
    slice(1) %>%
    ungroup() %>%
    arrange(KP + id_reach_CAL) # Because order is upstream to downstream and KP is increasing


# Check Q_measured of each event at each station or location with measurements
for (seq_event in all_available_events) {
    # the event exist in the event to be used during calibration
    if (seq_event %in% all_events_cal) {
        seq_stations <- as.vector(id_station_var[[seq_event]]$station %>% na.omit())
        for (seq_measure in seq_along(seq_stations)) {
            check_simulation_time(
                MAGE_main_folder = MAGE_main_folder,
                mage_projet_name = mage_projet_name,
                Observations = Q_obs %>% filter(station == id_station_var[[seq_event]]$station[seq_measure]),
                event = all_events_cal[match(seq_event, all_events_cal)] # Manual modification to check
            )
        }
    }
}

ggplot(Q_obs, aes(x = time, y = variable, color = id_event)) +
    geom_point() +
    geom_errorbar(aes(ymin = variable - 2 * uncertainty, ymax = variable + 2 * uncertainty)) +
    facet_wrap(~station) +
    theme_bw()

################################
# WSE observations (WSE_obs)
################################

# Metadata: link between real dates and simulation period

# Fixed value to the uncertainty of measured WSE
sd_WSE_fixed <- 0.05 # 5 cm, so 10 cm 2*sigma

############################
# WSE measurements
############################
# Initial date of the modeling
init_model_date <- key_info_event_extraction %>%
    filter(id_event %in% all_events_cal & var == "WSE") %>%
    select(
        id_event,
        station,
        var,
        model_start
    )

if (nrow(unique(key_info_event_extraction %>%
    filter(id_event %in% all_events_cal & var == "WSE") %>%
    select(
        var,
        model_start
    ))) != length(unique(key_info_event_extraction$id_event))) {
    stop("The observed WSE must have the same date for a same event, if not a new event must be done")
}


H_obs <-
    key_info_event_extraction %>%
    filter(
        var == "WSE",
        id_event %in% all_events_cal
    ) %>%
    select(id_event,
        target_datetime = measured_start
    ) %>%
    distinct() %>%
    # Create every event × target_datetime × station combination
    crossing(
        all_H_obs %>%
            distinct(station)
    ) %>%
    # Bring all observations for the corresponding station
    left_join(
        all_H_obs,
        by = "station",
        relationship = "many-to-many"
    ) %>%
    # Calculate distance between observed time and desired WSE time
    mutate(
        abs_diff = abs(
            as.numeric(
                difftime(
                    datetime_TU,
                    target_datetime,
                    units = "secs"
                )
            )
        )
    ) %>%
    # Tolerance of 10 minute of difference between target datetime and observed value
    filter(abs_diff <= 10 * 60) %>%
    # Keep ONLY the closest observation
    group_by(id_event, target_datetime, station) %>%
    slice_min(abs_diff, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    # Add station information
    left_join(
        station_KP,
        by = "station"
    ) %>%
    rename(
        id_reach_CAL = id_reach_HM
    ) %>%
    # Bring the time of model start
    left_join(
        init_model_date %>%
            distinct(id_event, model_start),
        by = "id_event"
    ) %>%
    # Add event information
    mutate(
        gauging_adcp = NA,
        var = "WSE",
        uncertainty = sd_WSE_fixed,
        station = NA,
        time = as.numeric(
            difftime(
                datetime_TU,
                model_start,
                units = "secs"
            )
        ),
        event = match(id_event, all_events_cal),
        name_event = unname(name_events[id_event])
    ) %>%
    left_join(
        typology_lookup,
        by = "id_reach_CAL"
    ) %>%
    select(-abs_diff) %>%
    select(
        datetime_TU,
        variable,
        uncertainty,
        gauging_adcp,
        id_event,
        station,
        var,
        KP,
        id_reach_CAL,
        model_start,
        target_datetime,
        time,
        event,
        name_event,
        Reach_groupped_Cal
    ) %>%
    mutate(set = "calibration") %>%
    arrange(target_datetime, KP)

ggplot(H_obs, aes(x = KP, y = variable, color = id_event)) +
    geom_point() +
    geom_errorbar(aes(ymin = variable - 2 * uncertainty, ymax = variable + 2 * uncertainty)) +
    facet_wrap(~ target_datetime + id_event) +
    theme_bw()


# All calibration data
observed_data <- rbind(
    Q_obs,
    H_obs
)

##################################
# !!!!!!!!!!!!
# Ensure that order should be in coherence with mage model set up!!!!!!!!!!!!
# !!!!!!!!!!!!

## just to know the order of calibration data
Cal_measures <-
    list(
        X90_2015 =
            data.frame(
                event = observed_data$event,
                reach = observed_data$id_reach_CAL,
                x = observed_data$KP,
                t = observed_data$time
            )
    )

observed_data <- observed_data %>%
    rename(
        "reach" = "id_reach_CAL",
        "x" = "KP",
        "t" = "time"
    )

if (any(duplicated(observed_data %>% select(t, event, x, var, datetime_TU, station)))) stop("No duplicated data is available in observed data")

save(observed_data, file = "data/processed_data/Lower_Seine/Seine_2015/observed_data.RData")
