cat("\014")
list <- ls()

library(stringr)
library(dplyr)

cas <- "smooth_bed" # rough_bed smooth_bed

path <- file.path("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/data-raw/HHLab/", cas)
path_results <- file.path("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/processed_data/HHLab/", cas)

cases <- list.files(path)[str_detect(list.files(path), "case")]

string_h_1500 <- "_flow_depth_y_1500_mm.txt"
string_h_500 <- "_flow_depth_y_500_mm.txt"
string_parameters <- "_mixing_layer_mean_flow_parameters.csv"

So <- 1.04e-3

all_data_calibration <- list(
    WSE = NA,
    Discharge = NA,
    Velocity = NA
)
WS_profiles <- list()
j <- 1
for (i in cases) {
    path_relative <- file.path(path, i)

    h500 <- read.table(file.path(path_relative, paste0(i, string_h_500)),
        header = TRUE,
        colClasses = c("numeric", "numeric")
    )
    colnames(h500) <- c("x", "h_500")

    h1500 <- read.table(file.path(path_relative, paste0(i, string_h_1500)), header = TRUE)
    colnames(h1500) <- c("x", "h_1500")

    param <- read.table(file.path(path_relative, paste0(i, string_parameters)), sep = ";", row.names = NULL)

    if (!identical(h1500$x, h500$x)) stop("Measures are not taken at the same distances or a data is missing")

    WS_profile_temp <- merge(h500, h1500, by = "x")
    WS_profile_temp[, -1] <- WS_profile_temp[, -1] / 1000 # Pass to meter unity
    WS_profile_temp <- WS_profile_temp %>%
        mutate(h_mean = (h_500 + h_1500) / 2)

    WS_profile <- WS_profile_temp[c("x", "h_mean")]

    str_colonnes <- c("x (m)", "h (mm)")
    param_extraited <- t(param[which(param[, 1] %in% str_colonnes), ])
    param_extraited <- data.frame(param_extraited[-1, ])
    colnames(param_extraited) <- c("x", "h_mean")
    param_extraited$h_mean <- as.numeric(param_extraited$h_mean) / 1000 # Pass to meter unity

    param_extraited <- param_extraited %>%
        mutate(across(everything(), as.numeric))
    #     mutate(ks_mean = sqrt(8 * 9.81 / cf_mean) * 1 / (h_mean**(1 / 6)))

    data_WS_profile <- full_join(WS_profile, param_extraited)
    data_WS_profile <- data_WS_profile[order(data_WS_profile$x), ]

    # Correct water depth to water level from last point recorded
    z <- 1
    data_WS_profile$z_riverbed <- NA

    for (p in data_WS_profile$x) {
        data_WS_profile$z_riverbed[z] <- So * (data_WS_profile$x[nrow(data_WS_profile)] - p)
        z <- z + 1
    }

    data_WS_profile$z_mean <- data_WS_profile$z_riverbed + data_WS_profile$h_mean

    model <- lm(formula = z_mean ~ x, data = data_WS_profile)
    residuals <- resid(model) # or model$residuals
    sd_residuals <- sd(residuals)

    data_WS_profile$Yu <- sd_residuals
    WS_profiles[[j]] <- data_WS_profile

    names(WS_profiles)[j] <- i

    all_data_calibration$WSE <- data_WS_profile

    save(all_data_calibration, file = file.path(path_results, paste0("data_HHLab_", i, ".RData")))

    if (i == last(cases)) {
        all_data_calibration$WSE <- WS_profiles
        save(all_data_calibration, file = file.path(path_results, paste0("data_HHLab_all_cases.RData")))
    }
    j <- j + 1
}
