cat("\014")
rm(list = ls())

library(stringr)
library(dplyr)
library(readODS)


path <- "/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/data-raw/HHLab/compound_transition_friction"
path_results <- "/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/processed_data/HHLab/compound_transition_friction/"

string_depth <- "CWM_CMW_lignes_d_eau_modif.ods"

So <- 1.05e-3
delta_h_floodplain <- 0.115
channel_long <- 18

all_data_calibration <- list(
    WSE = NA,
    Discharge = NA,
    Velocity = NA
)
WS_profiles <- list()

# relative zero
relative_zero <- data.frame(case = c("CWM", "CMW"), x_0 = c(9050, 9150))

# Flow depth
raw_data <- read_ods(file.path(path, string_depth))



# Merge raw_data with relative_zero using str_detect
raw_data_with_zero <- raw_data %>%
    mutate(
        # Find the matching x_0 for each ID_experiment
        x_0 = relative_zero$x_0[
            sapply(ID_experiment, function(id) {
                which(str_detect(id, relative_zero$case))
            })
        ]
    )

# Adjust x (mm) by adding x_0, except for "CWMQ12 floodplain"
raw_data_adjusted <- raw_data_with_zero %>%
    mutate(
        `x (mm)` = ifelse(ID_experiment == "CWMQ12" & ID_channel == "Floodplain", `x (mm)`, `x (mm)` + x_0)
    ) %>%
    select(-x_0) %>% # Remove the x_0 column if you don't need it
    rename(x = "x (mm)", H = "H_m (mm)")

library(ggplot2)
ggplot(raw_data_adjusted, aes(x = x, y = H, col = ID_experiment)) +
    geom_point() +
    facet_wrap(~ID_channel)

data_H_mean_separated <- raw_data_adjusted %>%
    group_by(ID_experiment, x, ID_channel) %>%
    summarise(
        h_mean = mean(H) / 1000,
        x = x / 1000,
        .groups = "drop"
    )

ggplot(data_H_mean_separated, aes(x = x, y = h_mean, col = ID_experiment)) +
    geom_point()

# Correct water depth to water level from last point recorded

data_H_mean_separated$z_riverbed <- NA

cases <- unique(data_H_mean_separated$ID_experiment)

for (i in cases) {
    mask <- which(data_H_mean_separated$ID_experiment == i)
    mask_data <- data_H_mean_separated[mask, ]
    z <- 1
    for (p in mask_data$x) {
        id <- mask[z]
        data_H_mean_separated$z_riverbed[id] <- So * (channel_long - p)
        z <- z + 1
    }
}

ggplot(data_H_mean_separated, aes(x = x, y = h_mean, col = ID_channel)) +
    geom_point() +
    geom_line(aes(y = z_riverbed)) +
    facet_wrap(~ID_experiment)

data_H_mean_separated$h_from_riverbed <- data_H_mean_separated$h_mean


data_H_mean_separated$z_mean <- data_H_mean_separated$z_riverbed + data_H_mean_separated$h_from_riverbed

data_WS_profile <- data.frame(data_H_mean_separated)

model <- lm(formula = h_from_riverbed ~ x, data = data_WS_profile)
residuals <- resid(model) # or model$residuals
sd_residuals <- sd(residuals)

data_WS_profile$Yu <- sd_residuals

all_data_calibration$WSE <- data_WS_profile
save(all_data_calibration, file = file.path(path_results, "data_HHLab_uniform_case.RData"))
