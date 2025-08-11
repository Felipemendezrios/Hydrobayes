cat("\014")
rm(list = ls())

library(stringr)
library(dplyr)

path <- "/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/data-raw/HHLab/compound_single_friction"
path_results <- "/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/processed_data/HHLab/compound_single_friction/"

cases <- list.files(path)[str_detect(list.files(path), "case")]

string_flow_depth_MC <- "Flowdepth_MC_Qf_08_at_y_star_1p2_1p5_1p8.csv"
string_flow_depth_FP <- "Flowdepth_FP_Qf_08_at_y_star_0p3_and_0p7.csv"
string_velocity <- "Mean_Velocity_Qf_08_at_z_star_0p94_x_star_various.csv"

So <- 1.1e-3
delta_h_floodplain <- 0.117
channel_long <- 18

all_data_calibration <- list(
    WSE = NA,
    Discharge = NA,
    Velocity = NA
)
WS_profiles <- list()
j <- 1
for (i in cases) {
    path_relative <- file.path(path, i)

    # Flow depth

    # Main channel (MC)
    h_MC <- read.table(file.path(path_relative, string_flow_depth_MC),
        header = TRUE,
        colClasses = c("numeric", "numeric", "numeric"),
        skip = 1,
        sep = ","
    )
    colnames(h_MC) <- c("y_star_m", "x", "h_mm")

    h_mean_MC_by_x_star <- data.frame(h_MC %>%
        group_by(x) %>%
        summarise(
            h_mean = mean(h_mm / 1000, na.rm = TRUE),
            ID = "MC"
        ))

    # Floodplain (FP)
    h_FP <- read.table(file.path(path_relative, string_flow_depth_FP),
        header = TRUE,
        colClasses = c("numeric", "numeric", "numeric"),
        skip = 1,
        sep = ","
    )
    colnames(h_FP) <- c("y_star_m", "x", "h_mm")

    h_mean_FP_by_x_star <- data.frame(
        h_FP %>%
            group_by(x) %>%
            summarise(h_mean = mean(h_mm / 1000, na.rm = TRUE)),
        ID = "FP"
    )

    h_mean_MC_by_x_star <- h_mean_MC_by_x_star[order(h_mean_MC_by_x_star$x), ]
    h_mean_FP_by_x_star <- h_mean_FP_by_x_star[order(h_mean_FP_by_x_star$x), ]

    # Correct water depth to water level from last point recorded
    z <- 1
    h_mean_FP_by_x_star$z_riverbed <- h_mean_MC_by_x_star$z_riverbed <- NA

    for (p in h_mean_FP_by_x_star$x) {
        h_mean_FP_by_x_star$z_riverbed[z] <- So * (channel_long - p)
        z <- z + 1
    }

    z <- 1

    for (p in h_mean_MC_by_x_star$x) {
        h_mean_MC_by_x_star$z_riverbed[z] <- So * (channel_long - p)
        z <- z + 1
    }

    h_mean_MC_by_x_star$z_mean <- h_mean_MC_by_x_star$z_riverbed + h_mean_MC_by_x_star$h_mean
    h_mean_FP_by_x_star$z_mean <- h_mean_FP_by_x_star$z_riverbed + h_mean_FP_by_x_star$h_mean + delta_h_floodplain

    model <- lm(formula = z_mean ~ x, data = h_mean_MC_by_x_star)
    residuals <- resid(model) # or model$residuals
    sd_residuals <- sd(residuals)

    h_mean_MC_by_x_star$Yu <- sd_residuals

    model <- lm(formula = z_mean ~ x, data = h_mean_FP_by_x_star)
    residuals <- resid(model) # or model$residuals
    sd_residuals <- sd(residuals)

    h_mean_FP_by_x_star$Yu <- sd_residuals

    data_WS_profile <- rbind(h_mean_MC_by_x_star, h_mean_FP_by_x_star)

    all_data_calibration$WSE <- data_WS_profile
    save(all_data_calibration, file = file.path(path_results, "data_HHLab_uniform_case.RData"))
}
