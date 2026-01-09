rm(list = ls())
cat("\014")

load("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/processed_data/Rhone/Boundary_conditions/All_boundary_conditions/metadata.RData")


metadata_modified <- data.frame(
    files = file.path("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/processed_data/Rhone/Boundary_conditions/All_boundary_conditions/", paste0(metadata$export_name, ".RData")),
    export_name = metadata$export_name
)

setwd("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/processed_data/Rhone/Boundary_conditions/Specific_boundary_conditions/")

library(data.table)
library(ggplot2)
library(dplyr)

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

    interval_id <- rep(NA_integer_, nrow(Data))
    for (i in seq_len(nrow(range_dates_measures))) {
        idx <- between(
            Data$Date,
            range_dates_measures$start[i],
            range_dates_measures$end[i]
        )

        interval_id[idx] <- i
    }
    Data_subset$id_time_interval <- interval_id[keep]

    out_path_RData <- file.path(out_path_directory, paste0(export_name, "_subset.RData"))
    out_path_csv <- file.path(out_path_directory, paste0(export_name, "_subset.csv"))
    if (dir.exists(out_path_directory)) {
        unlink(file.path(out_path_RData), recursive = TRUE, force = TRUE)
        unlink(file.path(out_path_csv), recursive = TRUE, force = TRUE)
    } else {
        dir.create(out_path_directory, recursive = TRUE)
    }

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
        select(Date, t_minutes, Q, id_time_interval)

    save(Data_subset, file = out_path_RData)
    write.table(
        file = out_path_csv,
        Data_subset, row.names = FALSE, col.names = TRUE, sep = ";"
    )
    return(invisible(out_path))
}


# Range of dates searched
range_dates_measures <- data.frame(
    start = c(
        as.POSIXct("09/04/2016 00:00:00", format = "%d/%m/%Y %H:%M:%S", tz = "UTC"),
        as.POSIXct("26/03/2023 00:00:00", format = "%d/%m/%Y %H:%M:%S", tz = "UTC")
    ),
    end = c(
        as.POSIXct("15/04/2016 23:00:00", format = "%d/%m/%Y %H:%M:%S", tz = "UTC"),
        as.POSIXct("31/03/2023 23:00:00", format = "%d/%m/%Y %H:%M:%S", tz = "UTC")
    )
)

range_dates_measures <- data.frame(
    start = c(
        as.POSIXct("22/03/2023 09:00:00", format = "%d/%m/%Y %H:%M:%S", tz = "UTC")
    ),
    end = c(
        as.POSIXct("31/08/2023 22:00:00", format = "%d/%m/%Y %H:%M:%S", tz = "UTC")
    )
)

if (any(range_dates_measures$end <= range_dates_measures$start)) {
    stop("increasing interval values")
}


out_path_directory <- paste0(
    nrow(range_dates_measures),
    "_event_", substr(
        range_dates_measures[1, 1],
        1, 10
    ),
    "_to_",
    substr(
        range_dates_measures[nrow(range_dates_measures), 2],
        1, 10
    )
)

results_files <- vector("character", nrow(metadata_modified))

for (i in seq_len(nrow(metadata_modified))) {
    path_data <- metadata_modified$files[i]
    export_name <- metadata_modified$export_name[i]
    results_files[i] <- extract_data(
        path_data = path_data,
        range_dates_measures = range_dates_measures,
        out_path = out_path_directory,
        export_name = export_name
    )
}
