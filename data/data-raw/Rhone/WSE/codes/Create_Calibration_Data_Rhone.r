rm(list = ls())
graphics.off()
cat("\014")

library(dplyr)
library(tidyr)
library(lubridate)

setwd("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/data-raw/Rhone/WSE/")

observation_data <- "Lignes_deau_Rhone.csv"

obs_raw <- read.csv(
    file = file.path("pre_processing", observation_data),
    header = TRUE,
    row.names = NULL,
    sep = ","
)


LGN_BOU_ST_data <- read.table(
    file = file.path("pre_processing", "model_ST_LGN_BOU.txt"),
    header = TRUE, sep = " "
)
LGN_BOU_ST_data <- LGN_BOU_ST_data %>%
    mutate(id_reach = "LGN_BOU")

BOU_ANT_ST_data <- read.table(
    file = file.path("pre_processing", "model_ST_BOU_ANT.txt"),
    header = TRUE, sep = " "
)

BOU_ANT_ST_data <- BOU_ANT_ST_data %>%
    mutate(id_reach = "BOU_ANT")

ANT_JNS_ST_data <- read.table(
    file = file.path("pre_processing", "model_ST_ANT_JNS.txt"),
    header = TRUE, sep = " "
)

ANT_JNS_ST_data <- ANT_JNS_ST_data %>%
    mutate(id_reach = "ANT_JNS")

all_data_model <- rbind(
    LGN_BOU_ST_data,
    BOU_ANT_ST_data,
    ANT_JNS_ST_data
)

model_profile_id_reach_KP <- all_data_model %>%
    group_by(profile_id, id_reach) %>%
    filter(label == "axe") %>%
    summarise(
        KP = first(header_val),
        Z_thalweg = first(Z),
        .groups = "drop" # This ungroups after summarising
    ) %>%
    as.data.frame()
#########################################
# Add time of the measurements. If necessary, data are available here
# Q(t) at Lagnieu: https://bdoh.inrae.fr/OBSERVATOIRE-DES-SEDIMENTS-DU-RHONE/LAGNIEU/DEB
###################

# WSE measurements have a date time underlied manually. Here the information associated to each WSE

info_date_wse_measurements <-
    data.frame(
        Q = c(
            158,
            300,
            525,
            750,
            1350
        ),
        Date_time = c(
            as.POSIXct(
                "13/02/2011 12:00:00",
                format = "%d/%m/%Y %H:%M:%S",
                tz = "UTC"
            ),
            as.POSIXct(
                "08/10/2008 14:00:00",
                format = "%d/%m/%Y %H:%M:%S",
                tz = "UTC"
            ),
            as.POSIXct(
                "24/03/2009 15:00:00",
                format = "%d/%m/%Y %H:%M:%S",
                tz = "UTC"
            ),
            as.POSIXct(
                "26/02/2010 13:00:00",
                format = "%d/%m/%Y %H:%M:%S",
                tz = "UTC"
            ),
            as.POSIXct(
                "08/12/2010 10:00:00",
                format = "%d/%m/%Y %H:%M:%S",
                tz = "UTC"
            )
        )
    )

# I'm not sure about it, so i prefer to assign the same date and time to all measurements

# speed_measurement <- 11 # km/h information obtained from WSE measurements from AIN with 90 m3/s
# distance_measurement <- 35.150 # km
# duration <- distance_measurement / speed_measurement # in hours

obs_merged <- merge(
    obs_raw,
    info_date_wse_measurements,
    by.x = "ID_Cal",
    by.y = "Q",
    all.x = TRUE
) %>% select(
    -Date
)

thalweg_model <- all_data_model[which(all_data_model$label == "axe"), ]

# Compare thalweg to verify if there is any huge difference between observed data and modeled data
ggplot(obs_merged %>% filter(ID_Cal == "Thalweg"), aes(x = KP, y = Cote, col = "obs")) +
    geom_line() +
    geom_line(data = thalweg_model, aes(x = header_val, y = Z, col = "model"))

# OK, thalweg are similar, so it does not necessary any correction

obs_filtered <- obs_merged %>% filter(ID_Cal != "Thalweg")



KP_to_remove <- obs_filtered$KP[
    which(!obs_filtered$KP %in% thalweg_model$header_val)
]

obs_assigned_KP_model <-
    obs_filtered %>%
    filter(!KP %in% KP_to_remove) %>%
    rename(
        id_case = ID_Cal,
        WSE = Cote,
        time = Date_time
    )

ggplot(
    obs_assigned_KP_model,
    aes(x = KP, y = WSE, col = id_case)
) +
    geom_point() +
    geom_line(data = thalweg_model, aes(x = header_val, y = Z, col = "model"))

if (length(which(!obs_assigned_KP_model$KP %in% thalweg_model$header_val)) != 0) stop("something is wrong")
##########################################################
all_WSE_Rhone <-
    merge(model_profile_id_reach_KP,
        obs_assigned_KP_model,
        by = "KP"
    ) %>%
    select(profile_id, KP, Z_thalweg, id_case, WSE, time, id_reach)

ggplot(
    all_WSE_Rhone,
    aes(x = KP, y = WSE, col = id_case)
) +
    geom_point() +
    geom_line(aes(x = KP, y = Z_thalweg, col = "model"))



save(all_WSE_Rhone, file = "/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/processed_data/Rhone/all_observations_Rhone.RData")
