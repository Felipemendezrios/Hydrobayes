rm(list = ls())
cat("\014")

load("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/processed_data/Rhone/Boundary_conditions/All_boundary_conditions/metadata.RData")


metadata_modified <- data.frame(
    files = file.path("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/processed_data/Rhone/Boundary_conditions/All_boundary_conditions", paste0(metadata$export_name, ".RData")),
    export_name = metadata$export_name
)

setwd("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/processed_data/Rhone/Boundary_conditions/Specific_boundary_conditions/")

library(data.table)
library(ggplot2)
library(dplyr)
library(lubridate)

load_RData_as <- function(path, new_name, envir = .GlobalEnv) {
    e <- new.env()
    obj <- load(path, envir = e)
    assign(new_name, e[[obj]], envir = envir)
}
extract_data <- function(
    path_data,
    range_dates_measures,
    out_path,
    export_name) {
    # Load and rename safely
    load_RData_as(
        path_data,
        new_name = "Data",
        envir = environment()
    )
    # Build logical keep vector
    keep <- Reduce(`|`, lapply(seq_len(nrow(range_dates_measures)), function(i) {
        between(
            Data$Date,
            range_dates_measures$start[i],
            range_dates_measures$end[i]
        )
    }))

    Data_subset <- Data[keep, ]

    interval_id <- name_id <- rep(NA_integer_, nrow(Data))
    for (i in seq_len(nrow(range_dates_measures))) {
        idx <- between(
            Data$Date,
            range_dates_measures$start[i],
            range_dates_measures$end[i]
        )

        interval_id[idx] <- i
        name_id[idx] <- range_dates_measures$id_case[i]
    }
    Data_subset$id_time_interval <- interval_id[keep]
    Data_subset$name_id <- name_id[keep]

    out_path_RData <- file.path(out_path, paste0(export_name, "_subset.RData"))
    out_path_csv <- file.path(out_path, paste0(export_name, "_subset.csv"))

    # Convert the discharge in numeric format
    Data_subset$Q <- as.numeric(gsub(",", ".", Data_subset$Q))

    # Add interval id to the ranges
    range_dates_measures <- range_dates_measures %>%
        mutate(id_time_interval = row_number())

    Data_subset <- Data_subset %>%
        left_join(
            range_dates_measures %>% select(id_time_interval, start),
            by = "id_time_interval"
        ) %>%
        group_by(id_time_interval) %>%
        arrange(Date) %>%
        mutate(
            t_minutes = as.numeric(difftime(Date, start, units = "mins"))
        ) %>%
        ungroup() %>%
        select(Date, t_minutes, Q, id_time_interval, name_id)

    save(Data_subset, file = out_path_RData)
    write.table(
        file = out_path_csv,
        Data_subset, row.names = FALSE, col.names = TRUE, sep = ";"
    )
    return(invisible(out_path))
}



######
# Load measurements to get the range of dates

# AIN
load("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/processed_data/Rhone/all_observations_AIN.RData")
load("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/processed_data/Rhone/all_observations_Rhone.RData")

range_dates_measures_AIN <- all_WSE_Ain %>%
    group_by(id_case) %>%
    summarise(min(time), max(time)) %>%
    mutate(id_case = paste0("AIN_", id_case))

range_dates_measures_RHONE <- all_WSE_Rhone %>%
    group_by(id_case) %>%
    summarise(min(time), max(time)) %>%
    mutate(id_case = paste0("RHONE_", id_case))

range_dates_measures <- rbind(range_dates_measures_AIN, range_dates_measures_RHONE)

# Add extension of the time series to warn up the model
extension_time <- days(5)
# Range of dates searched
range_dates_measures <- range_dates_measures %>%
    rename(
        start = `min(time)`,
        end = `max(time)`
    ) %>%
    mutate(
        # id_case = as.numeric(id_case),
        start = start - extension_time,
        end = end + extension_time,
    ) %>%
    arrange(id_case)

if (any(range_dates_measures$end <= range_dates_measures$start)) {
    stop("increasing interval values")
}

results_files <- vector("character", nrow(metadata_modified))

for (i in seq_len(nrow(metadata_modified))) {
    path_data <- metadata_modified$files[i]
    export_name <- metadata_modified$export_name[i]
    results_files[i] <- extract_data(
        path_data = path_data,
        range_dates_measures = range_dates_measures,
        out_path = getwd(),
        export_name = export_name
    )
}
