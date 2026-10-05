rm(list = ls()) # Clean workspace

library(stringr)
library(ggplot2)


setwd("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/data-raw/Lower_Seine/Observations/Discharge/Piney_2015/")

path_raw_data_oursin_estimation <- "Intermedial_results/1_files_generated_python"

path_intermedial_files <- "Intermedial_results/2_files_generated_r"

# Pretreatment of raw data gauging

GaugingsFileName <- list.files(path_raw_data_oursin_estimation, pattern = ".csv")
ColumnNameTemp <- str_split(GaugingsFileName, pattern = "_")
ResultNameTemp <- str_split(GaugingsFileName, pattern = ".csv")

Gaugingdata <- c()
Gaugingdata_raw <- c()
Gaugingdata_filtered <- c()
gauging_data_merge_temp <- c()
gauging_data_2015 <- c()
u2_mb <- 0 # u_mb = keep default value
u2_cov <- 0 # u_cov = 3.6% (already solved in QRevInt for negative discharge values)



# Add constant component to the uncertainty calculation
# Apply first filter: transect with less than 1 minute to cross the section (too quick for, higher uncertainties, worst discharge estimation)
# Correct UTC variation during discharge measurements: M9 Aizier -3600 RG Aizier + 3600
# Second filter: filter manually some transect which are not coherent at all

for (i in 1:length(GaugingsFileName)) {
   # setwd(path_raw_data_oursin_estimation)

   raw_data <- read.csv(file.path(path_raw_data_oursin_estimation, GaugingsFileName[i]), header = T, sep = ";")
   data_temp <- raw_data[c(-1)]
   ColumnName <- ColumnNameTemp[[i]][c(1, 2)]
   ColumnNameF <- paste0(ColumnName, collapse = "_")
   ResultName <- ResultNameTemp[[i]][1]

   for (j in 1:3) {
      data_temp[, j] <- as.POSIXct(strptime(as.character(data_temp[, j]),
         format = "%d/%m/%Y %H:%M", tz = "GMT"
      ))
   }

   data_temp <- data.frame(data_temp,
      u_s = (data_temp$oursin / 2), # Verify if QRevInt returns the uncertainty at 95% confidence interval
      delta_t = (as.numeric(difftime(data_temp$t_end, data_temp$t_start, units = "secs"))),
      delta_Qt = rep(NA, nrow(data_temp)),
      delta_tt = rep(NA, nrow(data_temp)),
      u_ns = rep(NA, nrow(data_temp)),
      # u_oursin_modifie = rep(0, nrow(data_temp)), # Remove if u2_cov and u2_mb = 0
      u_total = rep(NA, nrow(data_temp))
   )

   data_temp <- data_temp[order(data_temp$t_start), ]

   id_Q_neg <- which(data_temp$discharge < 0)
   id_Q_pos <- which(data_temp$discharge >= 0)
   # data_temp$u_oursin_modifie[id_Q_neg] <- sqrt(data_temp$u_s[id_Q_neg] + u2_cov + u2_mb)
   # data_temp$u_oursin_modifie[id_Q_pos] <- sqrt(data_temp$u_s[id_Q_pos] + u2_mb)

   id_Qt <- which(colnames(data_temp) == c("delta_Qt"))
   id_tt <- which(colnames(data_temp) == c("delta_tt"))
   id_uq <- which(colnames(data_temp) == c("u_ns")) # Call the good name of the column
   id_u_q_i <- which(colnames(data_temp) == c("u_total"))

   for (j in 1:nrow(data_temp)) {
      if (j == 1) {
         delta_Qt <- data_temp$discharge[j + 1] - data_temp$discharge[j]
         delta_tt <- as.numeric(difftime(data_temp$t_moy[j + 1], data_temp$t_moy[j], units = "secs"))
         u_ns <- (1 / (2 * sqrt(3) * data_temp$discharge[j])) * (data_temp$delta_t[j] * (delta_Qt / delta_tt)) # k=1 standard deviation

         data_temp[j, id_Qt] <- delta_Qt
         data_temp[j, id_tt] <- delta_tt
         data_temp[j, id_uq] <- u_ns
      } else if (j == nrow(data_temp)) {
         delta_Qt <- data_temp$discharge[j] - data_temp$discharge[j - 1]
         delta_tt <- as.numeric(difftime(data_temp$t_moy[j], data_temp$t_moy[j - 1], units = "secs"))
         u_ns <- (1 / (2 * sqrt(3) * data_temp$discharge[j])) * (data_temp$delta_t[j] * (delta_Qt / delta_tt))

         data_temp[j, id_Qt] <- delta_Qt
         data_temp[j, id_tt] <- delta_tt
         data_temp[j, id_uq] <- u_ns
      } else {
         data_q <- t(data.frame(
            "1" = data_temp$discharge[j - 1],
            "2" = data_temp$discharge[j],
            "3" = data_temp$discharge[j + 1]
         ))

         data_t <- t(data.frame(
            "1" = data_temp$t_moy[j - 1],
            "2" = data_temp$t_moy[j],
            "3" = data_temp$t_moy[j + 1]
         ))

         data_reg <- data.frame(
            "time" = data_t,
            "discharge" = data_q
         )
         data_reg$time <- as.POSIXct(strptime(as.character(data_reg$time),
            format = "%Y-%m-%d %H:%M", tz = "GMT"
         ))

         reg <- lm(discharge ~ time, data_reg)
         q_reg <- fitted(reg)
         t_reg <- data_reg$time

         delta_Qt <- q_reg[3] - q_reg[1]
         delta_tt <- as.numeric(difftime(t_reg[3], t_reg[1], units = "secs"))
         u_ns <- (1 / (2 * sqrt(3)) * data_temp$discharge[j]) * (data_temp$delta_t[j] * (delta_Qt / delta_tt))

         data_temp[j, id_Qt] <- delta_Qt
         data_temp[j, id_tt] <- delta_tt
         data_temp[j, id_uq] <- u_ns
      }
   }

   u_total <- sqrt(
      # ((data_temp$u_oursin_modifie / 100)^2) *
      (data_temp$u_s^2 + data_temp$u_ns^2) * (data_temp$discharge^2)
   )
   data_temp[, id_u_q_i] <- u_total

   # ######
   # ### Simplified calculation
   # ######

   # # Get from Qrevint the value from some ADCP transect and get a mean value
   # U_typical <- rep(c(0.1, 0.1, 0.1), 2) # m2
   # Aire_mouille <- rep(c(4800, 2120, 1380), 2) # m2
   # min_velocity <- c(0.05) # m/s

   # Q_min_measured <- round(Aire_mouille[i] * min_velocity / 5) * 5
   # print(Q_min_measured)
   # u_q_i_simplified <- sqrt(((U_typical[i])^2) * (data_temp$discharge^2) + Q_min_measured^2) / 2
   # data_temp[, id_u_q_i] <- u_q_i_simplified


   Gaugingdata_station_raw <- data.frame(
      data_temp$t_moy,
      data_temp$discharge,
      data_temp$u_total
   )

   colnames(Gaugingdata_station_raw) <- c(
      "Time",
      paste0("q_", ColumnNameF),
      paste0("u_q_", ColumnNameF)
   )


   if (length(which(data_temp$delta_t <= 60)) != 0) {
      data_tempF <- data_temp[c(-which(data_temp$delta_t <= 60)), ]
   } else {
      data_tempF <- data_temp
   }

   Gaugingdata_station <- data.frame(
      data_tempF$t_moy,
      data_tempF$discharge,
      data_tempF$u_total
   )

   colnames(Gaugingdata_station) <- c(
      "Time",
      paste0("q_", ColumnNameF),
      paste0("u_q_", ColumnNameF)
   )

   plot_qt <- ggplot() +
      geom_point(aes(x = data_tempF$t_moy, y = data_tempF$discharge)) +
      geom_errorbar(aes(
         x = data_tempF$t_moy,
         ymin = data_tempF$discharge - 2 * data_tempF$u_total,
         ymax = data_tempF$discharge + 2 * data_tempF$u_total
      )) +
      labs(
         x = "Time (hours)",
         y = "Q (m3/s)",
         title = ColumnNameF
      ) +
      theme_classic() +
      theme(
         plot.title = element_text(size = 14, hjust = 0.5),
         legend.position = "right",
         legend.title = element_text(size = 12),
         legend.text = element_text(size = 12),
         axis.text = element_text(size = 12),
         axis.title = element_text(size = 12)
      )

   ggsave(file.path(path_intermedial_files, "1_filter_time_gauging", paste0(ColumnNameF, "_first_filter.png")))

   # setwd(file.path(path_intermedial_files,'filter_time_gauging'))
   save(Gaugingdata_station, file = file.path(path_intermedial_files, "1_filter_time_gauging", paste0(ResultName, ".RData")))
   Gaugingdata[[i]] <- Gaugingdata_station
   Gaugingdata_raw[[i]] <- Gaugingdata_station_raw

   # Independent analyze station by station:
   # Remove manually incoherent transects and put all data in TU+1

   if (ColumnNameF == "M9_Aizier") {
      gauging_filtered <- Gaugingdata_station[-c(17, 19, 23), ]
      gauging_filtered$Time <- gauging_filtered$Time - 3600 # Correction of time zone
   } else if (ColumnNameF == "M9_Heurteauville") {
      gauging_filtered <- Gaugingdata_station[-c(64), ]
   } else if (ColumnNameF == "M9_Rouen") {
      gauging_filtered <- Gaugingdata_station[-c(6), ]
   } else if (ColumnNameF == "RG_Aizier") {
      gauging_filtered <- Gaugingdata_station[-c(22, 24, 25), ]
      gauging_filtered$Time <- gauging_filtered$Time + 3600 # Correction of time zone
   } else if (ColumnNameF == "RG_Heurteauville") {
      gauging_filtered <- Gaugingdata_station[, ]
   } else if (ColumnNameF == "RG_Rouen") {
      gauging_filtered <- Gaugingdata_station[, ]
   }

   save(gauging_filtered, file = file.path(path_intermedial_files, "2_filter_manually_correct_TU", paste0(ResultName, ".RData")))

   Gaugingdata_filtered[[i]] <- gauging_filtered

   plot_qt_f <- ggplot() +
      geom_point(aes(x = gauging_filtered[, 1], y = gauging_filtered[, 2])) +
      geom_errorbar(aes(
         x = gauging_filtered[, 1],
         ymin = gauging_filtered[, 2] - 2 * gauging_filtered[, 3],
         ymax = gauging_filtered[, 2] + 2 * gauging_filtered[, 3]
      )) +
      labs(
         x = "Time (hours)",
         y = "Q (m3/s)",
         title = ColumnNameF
      ) +
      theme_classic() +
      theme(
         plot.title = element_text(size = 14, hjust = 0.5),
         legend.position = "right",
         legend.title = element_text(size = 12),
         legend.text = element_text(size = 12),
         axis.text = element_text(size = 12),
         axis.title = element_text(size = 12)
      )

   ggsave(file.path(path_intermedial_files, "2_filter_manually_correct_TU", paste0(ColumnNameF, ".png")))
}


for (i in 1:length(Gaugingdata_filtered)) {
   temp <- Gaugingdata_filtered[[i]]
   colnames(temp) <- c("time", "variable", "uncertainty")
   # , "uncertainty_simplified")
   gauging_data_merge_temp[[i]] <- cbind(temp, gauging_adcp = rep(colnames(Gaugingdata_filtered[[i]][2]), nrow(Gaugingdata_filtered[[i]])))
}

for (i in 1:(length(gauging_data_merge_temp) / 2)) {
   gauging_data_2015[[i]] <- rbind(gauging_data_merge_temp[[i]], gauging_data_merge_temp[[i + 3]])
}


save(Gaugingdata_raw, file = file.path(path_intermedial_files, "Gaugings_M9_RD_all_raw_data.RData"))

# save(Gaugingdata_filtered, file = file.path("Final_results", "Gaugings_M9_RD_data.RData"))

save(gauging_data_2015, file = file.path("Final_results", "gauging_data_2015.RData"))

# Folder with processed data
# save(Gaugingdata_filtered, file = file.path("/home/famendezrios/Documents/papiers/Qmec/Data/Processed_data/Lower_Seine/Gaugings/Piney_2018/", "Gaugings_M9_RD_data.RData"))

save(gauging_data_2015, file = file.path("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/data-raw/Lower_Seine/Observations/Discharge/Piney_2015/gauging_data_2015.RData"))

# toto <- gauging_data_2015
# ggplot(toto[[3]], aes(x = time, y = variable, color = gauging_adcp)) +
#    geom_point() +
#    geom_errorbar(aes(
#       ymin = variable - 2 * uncertainty,
#       ymax = variable + 2 * uncertainty
#    )) +
#    labs(
#       x = "Time (hours)",
#       y = "Q (m3/s)"
#       # title =
#    ) +
#    theme_classic() +
#    theme(
#       plot.title = element_text(size = 14, hjust = 0.5),
#       legend.position = "right",
#       legend.title = element_text(size = 12),
#       legend.text = element_text(size = 12),
#       axis.text = element_text(size = 12),
#       axis.title = element_text(size = 12)
#    )
