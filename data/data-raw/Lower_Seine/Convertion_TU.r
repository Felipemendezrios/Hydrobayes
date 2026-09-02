rm(list = ls())
graphics.off()

library(ggplot2)
library(dplyr)

setwd("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data")

# Date of reference
date_ref <- as.POSIXct(
    "2015-09-28 01:10:00",
    format = "%Y-%m-%d %H:%M:%S",
    tz = "UTC"
)

date_ref2 <- as.POSIXct(
    "2015-09-30 23:55:00",
    format = "%Y-%m-%d %H:%M:%S",
    tz = "UTC"
)

diff <- difftime(date_ref2, date_ref, units = "secs")

all_categorize <- c("BC", "Observations_H", "Observations_Q")

# categorize = 'Observations_Q'
for (categorize in all_categorize) {
    if (categorize == "BC" | categorize == "Observations_H") {
        all_H_obs <- data.frame()
        if (categorize == "BC") {
            # Boundary conditions
            folder_read <- "data-raw/Lower_Seine/Boundary_conditions_model"
            folder_write <- "processed_data/Lower_Seine/boundaries_conditions"
            file_read <- "Tancarville_TU+1.txt"
            file_write_model <- "Tancarville_TU_model.txt"
            file_write_all <- "Tancarville_TU_all.txt"
        } else if (categorize == "Observations_H") {
            # Stage record
            folder_read <- "data-raw/Lower_Seine/Observations/Stage_records_TU+1"
            folder_write <- "processed_data/Lower_Seine/Observations"

            file_read <- list.files(folder_read)
        }

        # Counter
        idx <- 1
        for (station in file_read) {
            var_obs <- read.table(
                file.path(folder_read, station),
                header = FALSE,
                sep = ";",
                col.names = c("datetime_TU_1", "value")
            )

            var_obs$datetime_TU_1 <- as.POSIXct(
                var_obs$datetime_TU_1,
                format = "%Y-%m-%d %H:%M:%S",
                tz = "UTC"
            )

            # Put data in TU
            var_obs$datetime_TU <- var_obs$datetime_TU_1 - 3600

            var_obs$diff_sec <- as.numeric(
                difftime(var_obs$datetime_TU, date_ref, units = "secs")
            )
            ggplot(var_obs[which(var_obs$diff_sec >= -10000), ][1:100, ], aes(y = value)) +
                geom_line(aes(x = datetime_TU_1, color = "TU+1")) +
                geom_line(aes(x = datetime_TU, color = "TU")) +
                geom_vline(xintercept = date_ref) +
                theme_bw()

            var_obs <- var_obs %>%
                select(datetime_TU, diff_sec, value)
            if (categorize == "BC") {
                write.table(var_obs, file = file.path(folder_write, file_write_all[idx]), sep = ";", row.names = FALSE)

                var_obs_extract <- var_obs %>%
                    filter(diff_sec >= 0 & diff_sec <= diff) %>%
                    select(diff_sec, value)

                write.table(var_obs_extract, file = file.path(folder_write, file_write_model[idx]), sep = ";", row.names = FALSE)
            } else {
                var_obs <- var_obs %>%
                    mutate(id_station = sub("\\.txt$", "", station))
                all_H_obs <- rbind(all_H_obs, var_obs)
            }

            # Update counter
            idx <- idx + 1
        }
        if (categorize == "Observations_H") {
            save(all_H_obs, file = file.path(folder_write, "all_dataset_H_obs.RData"))
        }
    } else if (categorize == "Observations_Q") {
        main_folder_read <- "data-raw/Lower_Seine/Observations/Discharge"
        folder_write <- "processed_data/Lower_Seine/Observations"

        all_folder_Q <- list.dirs(main_folder_read, recursive = FALSE)

        all_obs_Q <- list()
        for (folder_read in all_folder_Q) {
            file_read <- list.files(folder_read, full.names = TRUE)
            for (id_gauging in file_read) {
                Q_gaugings <- get(load(id_gauging))

                # Index name to identify the author of the campaign
                source <- basename(folder_read)
                # Modify the lists
                Q_gaugings <- lapply(Q_gaugings, function(df) {
                    df$source <- source

                    if (source == "Piney_2015") {
                        # Correct time zone to stick to TU
                        # Measurements of Piney need to be correct because they are TU+1
                        df$time <- df$time - 3600
                        df$station <- sub(".*_", "", df$gauging_adcp)
                    } else {
                        df$station <- "Rouen"
                    }
                    df
                })

                # Add this list to the main list
                all_obs_Q[[length(all_obs_Q) + 1]] <- Q_gaugings
            }
        }
        all_obs_Q_df <- dplyr::bind_rows(
            lapply(all_obs_Q, dplyr::bind_rows)
        )
        save(all_obs_Q_df, file = file.path(folder_write, "all_dataset_Q_obs.RData"))
    } else {
        stop("Categorize does not exist")
    }
}
