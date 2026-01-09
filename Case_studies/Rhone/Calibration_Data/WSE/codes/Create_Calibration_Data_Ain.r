rm(list = ls())
graphics.off()
cat("\014")

library(dplyr)
library(tidyr)
library(lubridate)

setwd("/home/famendezrios/Documents/These/Cas_etude/Rhone/Calibration/Lignes_deau/")

reach_name <- "model_ST_PCH_PGA_with_observations.csv"

obs_raw <- read.csv(
    file = file.path("pre_processing", reach_name),
    header = TRUE
)

obs_raw <- obs_raw %>%
    mutate(
        # Paste the date with the time from X90_Heures
        X90_Date = dmy("12/04/2016") + hms(X90_Heures)
    )
# Check if all values are the same per profile_id
check_consistency <- obs_raw %>%
    group_by(profile_id) %>%
    summarise(
        header_val_unique = n_distinct(header_val),
        X90_Heures_unique = n_distinct(X90_Heures),
        X90_WSE_unique = n_distinct(X90_wse),
        X16_WSE_unique = n_distinct(X16_wse)
    )

# View profiles where there is inconsistency
inconsistent_profiles <- check_consistency %>%
    filter(header_val_unique > 1 | X90_WSE_unique > 1 | X16_WSE_unique > 1 | X90_Heures_unique > 1)

if (nrow(inconsistent_profiles) != 0) stop("inconsistent profiles")

obs_reach <- data.frame(obs_raw %>%
    group_by(profile_id) %>%
    summarise(
        KP = first(header_val), # assuming header_val is same for each profile
        X90_Date = first(X90_Date), # pick the first value (or use mean if needed)
        X90_WSE = first(X90_wse),
        X16_WSE = first(X16_wse)
    ))


time_duration <- as.numeric(difftime(
    obs_reach$X90_Date[length(obs_reach$X90_Date)],
    obs_reach$X90_Date[1],
    units = "secs"
))

# Add information of time to WSE measurements of 16 m3/s. In reported data, the measurement was in 2016, but it is not more information.
# We need to assume the date. Regarding data of 90 m3/s, 1 h 30 minutes were necessary to get the measurements.
# Some discharge time series were studied to defined the time period when they took the measurements.

# Q(t) at Port Galland: https://bdoh.inrae.fr/OBSERVATOIRE-DES-SEDIMENTS-DU-RHONE/PORT-GALLAND/DEB#measures-display
# Q(t) at Pont d'Ain: https://hydro.eaufrance.fr/stationhydro/V271201001/series
# Q(t) at Pont Chazey: https://hydro.eaufrance.fr/stationhydro/V294201001/series


###################

raw_qdt_station_1 <- read.csv(
    file = file.path(
        "/home/famendezrios/Documents/These/Cas_etude/Rhone/Calibration/Lignes_deau/Raw/QdT_station_Ain/",
        "QdT_pont_Ain_et_Pont_Chazey.csv"
    ),
    header = TRUE, sep = ";"
)

# V2712010 Pont_Ain
# V2942010 Chazey

qdt_station_1 <- raw_qdt_station_1 %>%
    mutate(
        STA_name = case_when(
            X.CdSiteHydro. == "V2942010" ~ "Pont_Chazey",
            X.CdSiteHydro. == "V2712010" ~ "Pont_Ain",
        )
    ) %>%
    select(
        X.DtObsHydro.,
        X.ResObsHydro.,
        STA_name
    ) %>%
    rename(
        "Date" = X.DtObsHydro.,
        "Q" = X.ResObsHydro.
    ) %>%
    slice(-1) %>%
    mutate(
        Date = as.POSIXct(Date,
            format = "%Y-%m-%d %H:%M:%S", tz = "UTC"
        ), # assuming format "01/01/2016 00:00:00"
        Q = as.numeric(Q) / 1000
    )

raw_qdt_station_2 <- read.table(
    file = file.path(
        "/home/famendezrios/Documents/These/Cas_etude/Rhone/Calibration/Lignes_deau/Raw/QdT_station_Ain/",
        "PORT-GALLAND_DEB.txt"
    ),
    header = TRUE, sep = ";", skip = 2
)

qdt_station_2 <- raw_qdt_station_2 %>%
    mutate(
        DateHeure_POSIX = dmy_hms(DateHeure),
        STA_name = "Port_Galland"
    ) %>%
    select(DateHeure_POSIX, Valeur, STA_name) %>%
    rename(
        "Date" = DateHeure_POSIX,
        "Q" = Valeur
    )

all_data_station_2016 <- rbind(qdt_station_1, qdt_station_2)

all_data_station_2016 <- all_data_station_2016 %>%
    filter(
        Date >= as.POSIXct("2016-01-01 00:00:00", tz = "UTC") &
            Date <= as.POSIXct("2016-12-31 23:59:00", tz = "UTC")
    )

library(ggplot2)

filtered_data <- all_data_station_2016 %>%
    filter(Q < 18 & Q > 14.5)

ggplot(filtered_data, aes(x = Date, y = Q, color = STA_name)) +
    geom_point(alpha = 0.7, size = 1)
labs(
    title = "Discharge over Time",
    x = "Date",
    y = "Discharge (Q)",
    color = "Station"
) +
    theme_minimal()

# Step 1: Identify the station with the fewest observations
station_counts <- filtered_data %>%
    group_by(STA_name) %>%
    summarise(n_obs = n()) %>%
    arrange(n_obs)

ref_station <- station_counts$STA_name[1] # station with fewest rows

# Step 2: Split data by station
data_list <- filtered_data %>%
    split(.$STA_name)

ref_data <- data_list[[ref_station]] # reference station data
other_stations <- setdiff(names(data_list), ref_station)

# Step 3: Check for each ref timestamp if other stations have a value within ±7 days
tolerance_days <- 2

library(purrr)

aligned <- ref_data %>%
    rowwise() %>%
    mutate(
        all_within_tolerance = all(
            map_lgl(other_stations, function(sta) {
                any(abs(data_list[[sta]]$Date - Date) <= days(tolerance_days))
            })
        )
    ) %>%
    ungroup()
# View timestamps where all stations have a value within ±1 week
aligned_filtered <- aligned %>% filter(all_within_tolerance)


# Step 1: Create intervals ±7 days around each aligned timestamp
intervals <- aligned_filtered %>%
    mutate(
        start = Date - days(7),
        end   = Date + days(7)
    ) %>%
    select(start, end)

# Step 2: Optional - combine overlapping intervals
# This gives a minimal set of continuous time ranges
intervals_combined <- intervals %>%
    arrange(start) %>%
    mutate(group = cumsum(start > lag(end, default = first(start)))) %>%
    group_by(group) %>%
    summarise(
        start = min(start),
        end   = max(end)
    ) %>%
    ungroup()

intervals_combined

# Filter filtered_data to only keep rows within any interval
filtered_aligned_data <- filtered_data %>%
    filter(
        map_lgl(Date, ~ any(.x >= intervals_combined$start & .x <= intervals_combined$end))
    )

selected_date <- as.POSIXct("2016-08-29 09:39:00", tz = "UTC")
ggplot(filtered_aligned_data, aes(x = Date, y = Q, color = STA_name)) +
    geom_point(alpha = 0.7, size = 1) +
    labs(
        title = "Discharge over Time",
        x = "Date",
        y = "Discharge (Q)",
        color = "Station"
    ) +
    theme_minimal() +
    geom_vline(xintercept = selected_date)

# The selected date is 2016-08-29 09:39:00 and the duration will be the same as the Q=90m3/s
# Number of rows
n <- nrow(obs_reach)

end_time_reach_1 <- selected_date + time_duration

# Generate evenly spaced sequence of POSIXct
obs_reach <- obs_reach %>%
    mutate(
        X16_Date = seq(from = selected_date, to = end_time_reach_1, length.out = n)
    )

ggplot(obs_reach, aes(x = X16_Date, y = X16_WSE)) +
    geom_point()
####################

obs_PCH_PGA <- data.frame(obs_reach %>%
    pivot_longer(
        cols = c(X90_WSE, X16_WSE, X90_Date, X16_Date),
        names_to = c("id_case", ".value"),
        names_pattern = "X(\\d+)_?(.*)"
    ) %>%
    rename(
        WSE = WSE,
        time = Date
    )) %>%
    mutate(
        id_reach = "PCH_PGA"
    )

ggplot(obs_PCH_PGA, aes(x = KP, y = WSE, col = id_case)) +
    geom_point()

ggplot(obs_PCH_PGA, aes(x = time, y = WSE, col = id_case)) +
    geom_point() +
    facet_wrap(~id_case, scales = "free_x")

save(obs_PCH_PGA, file = "processed/model_ST_PCH_PGA_with_observations.RData")
#######################################################################

reach_name <- "model_ST_PGA_CAIN_with_observations.csv"

obs_raw <- read.csv(
    file = file.path("pre_processing", reach_name),
    header = TRUE
)

obs_raw_clean <- obs_raw %>%
    # Convert times to hms, invalid times become NA
    mutate(time_hms = hms(X90_Heures)) %>%
    # Remove rows where parsing failed
    filter(!is.na(time_hms)) %>%
    # Combine with date
    mutate(X90_Date = dmy("12/04/2016") + time_hms) %>%
    # Drop the helper column
    select(-time_hms)

# Check if all values are the same per profile_id
check_consistency <- obs_raw_clean %>%
    group_by(profile_id) %>%
    summarise(
        header_val_unique = n_distinct(header_val),
        X90_Heures_unique = n_distinct(X90_Heures),
        X90_WSE_unique = n_distinct(X90_WSE),
        X16_WSE_unique = n_distinct(X16_wse)
    )

# View profiles where there is inconsistency
inconsistent_profiles <- check_consistency %>%
    filter(header_val_unique > 1 | X90_WSE_unique > 1 | X16_WSE_unique > 1 | X90_Heures_unique > 1)

if (nrow(inconsistent_profiles) != 0) stop("inconsistent profiles")

obs_reach <- data.frame(obs_raw_clean %>%
    group_by(profile_id) %>%
    summarise(
        KP = first(header_val), # assuming header_val is same for each profile
        X90_Date = first(X90_Date), # pick the first value (or use mean if needed)
        X90_WSE = first(X90_WSE),
        X16_WSE = first(X16_wse)
    ))


time_duration <- as.numeric(difftime(
    obs_reach$X90_Date[length(obs_reach$X90_Date)],
    obs_reach$X90_Date[1],
    units = "secs"
))
# The selected date is 2016-08-29 09:39:00 and the duration will be the same as the Q=90m3/s
# Number of rows
n <- nrow(obs_reach)

end_time_reach_2 <- end_time_reach_1 + time_duration

# Generate evenly spaced sequence of POSIXct
obs_reach <- obs_reach %>%
    mutate(
        X16_Date = seq(from = end_time_reach_1, to = end_time_reach_2, length.out = n)
    )

ggplot(obs_reach, aes(x = X16_Date, y = X16_WSE)) +
    geom_point()
####################

obs_PGA_CAIN <- data.frame(obs_reach %>%
    pivot_longer(
        cols = c(X90_WSE, X16_WSE, X90_Date, X16_Date),
        names_to = c("id_case", ".value"),
        names_pattern = "X(\\d+)_?(.*)"
    ) %>%
    rename(
        WSE = WSE,
        time = Date
    )) %>%
    mutate(
        id_reach = "PGA_CAIN"
    )

ggplot(obs_PGA_CAIN, aes(x = KP, y = WSE, col = id_case)) +
    geom_point()

ggplot(obs_PGA_CAIN, aes(x = time, y = WSE, col = id_case)) +
    geom_point() +
    facet_wrap(~id_case, scales = "free_x")

save(obs_PGA_CAIN, file = "processed/model_ST_PGA_CAIN_with_observations.RData")
#######################################################################
#######################################################################

reach_name <- "model_ST_CAIN_Confluence_with_observations.csv"

obs_raw <- read.csv(
    file = file.path("pre_processing", reach_name),
    header = TRUE
)

obs_raw_clean <- obs_raw %>%
    # Convert times to hms, invalid times become NA
    mutate(time_hms = hms(X90_Heures)) %>%
    # Remove rows where parsing failed
    filter(!is.na(time_hms)) %>%
    # Combine with date
    mutate(X90_Date = dmy("12/04/2016") + time_hms) %>%
    # Drop the helper column
    select(-time_hms)

# Check if all values are the same per profile_id
check_consistency <- obs_raw_clean %>%
    group_by(profile_id) %>%
    summarise(
        header_val_unique = n_distinct(header_val),
        X90_Heures_unique = n_distinct(X90_Heures),
        X90_WSE_unique = n_distinct(X90_wse),
        X16_WSE_unique = n_distinct(X16_wse)
    )

# View profiles where there is inconsistency
inconsistent_profiles <- check_consistency %>%
    filter(header_val_unique > 1 | X90_WSE_unique > 1 | X16_WSE_unique > 1 | X90_Heures_unique > 1)

if (nrow(inconsistent_profiles) != 0) stop("inconsistent profiles")

obs_reach <- data.frame(obs_raw_clean %>%
    group_by(profile_id) %>%
    summarise(
        KP = first(header_val), # assuming header_val is same for each profile
        X90_Date = first(X90_Date), # pick the first value (or use mean if needed)
        X90_WSE = first(X90_wse),
        X16_WSE = first(X16_wse)
    ))


time_duration <- as.numeric(difftime(
    obs_reach$X90_Date[length(obs_reach$X90_Date)],
    obs_reach$X90_Date[1],
    units = "secs"
))
# The selected date is 2016-08-29 09:39:00 and the duration will be the same as the Q=90m3/s
# Number of rows
n <- nrow(obs_reach)

end_time_reach_3 <- end_time_reach_2 + time_duration

# Generate evenly spaced sequence of POSIXct
obs_reach <- obs_reach %>%
    mutate(
        X16_Date = seq(from = end_time_reach_2, to = end_time_reach_3, length.out = n)
    )

ggplot(obs_reach, aes(x = X16_Date, y = X16_WSE)) +
    geom_point()
####################

obs_CAIN_Confluence <- data.frame(obs_reach %>%
    pivot_longer(
        cols = c(X90_WSE, X16_WSE, X90_Date, X16_Date),
        names_to = c("id_case", ".value"),
        names_pattern = "X(\\d+)_?(.*)"
    ) %>%
    rename(
        WSE = WSE,
        time = Date
    )) %>%
    mutate(
        id_reach = "CAIN_Confluence"
    )

ggplot(obs_CAIN_Confluence, aes(x = KP, y = WSE, col = id_case)) +
    geom_point()

ggplot(obs_CAIN_Confluence, aes(x = time, y = WSE, col = id_case)) +
    geom_point() +
    facet_wrap(~id_case, scales = "free_x")

save(obs_CAIN_Confluence, file = "processed/model_ST_CAIN_Confluence_with_observations.RData")
#######################################################################

all_WSE_Ain <- rbind(
    obs_PCH_PGA,
    obs_PGA_CAIN,
    obs_CAIN_Confluence
)


ggplot(all_WSE_Ain, aes(x = KP, y = WSE, col = id_reach)) +
    geom_point() +
    facet_wrap(~id_case)

ggplot(all_WSE_Ain, aes(x = time, y = WSE, col = id_reach)) +
    geom_point() +
    facet_wrap(~id_case, scales = "free_x")

save(all_WSE_Ain, file = "processed/all_observations_AIN.RData")
