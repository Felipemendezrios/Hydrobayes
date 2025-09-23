cat("\014")
rm(list = ls())

library(stringr)
library(dplyr)
library(lubridate)

path <- "/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/data-raw/Rhone/case_1_river_linear_extent"
path_results <- "/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/processed_data/Rhone/case_1_river_linear_extent"


# Common files for all experiments
riverbed <- read.table(file.path(path, "Thalweg.csv"), header = TRUE, sep = ";", dec = ",")
colnames(riverbed) <- c("KP", "Z_riverbed")
ZdQ_aval <- read.table(file.path(path, "CL_aval_ZdQ.csv"), header = TRUE, sep = ";", dec = ",")
colnames(ZdQ_aval) <- c("ref", "Q_LGN", "Z_JNS")

cases <- list.files(path)[str_detect(list.files(path), "case")]

ZdX_csv <- "ZdX.csv"
QdT_csv <- "QdT.csv"

all_data_calibration <- list(
    WSE = NA,
    Discharge = NA,
    Velocity = NA
)
WS_profiles <- list()
QdT_upstream <- list()
j <- 1
customized_names <- paste0("Q_", str_extract(cases, "(?<=Q_)\\d+(?=_Date)"))

for (i in cases) {
    path_relative <- file.path(path, i)
    # WSE measurements
    ZdX <- read.table(file.path(path_relative, ZdX_csv),
        header = TRUE,
        sep = ";",
        dec = ","
    )
    data_WS_profile <- ZdX[, -c(1, 2)]
    colnames(data_WS_profile) <- c("x", "z_observed")
    # Upstream discharge time series
    data_QdT_upstream <- read.table(file.path(path_relative, QdT_csv),
        header = TRUE,
        sep = ";",
        dec = ","
    )
    data_QdT_upstream <- data_QdT_upstream %>%
        mutate(
            Date = dmy_hms(paste(Date, Heure)),
            Lagnieu = as.numeric(Lagnieu),
            Ain = as.numeric(Ain),
        ) %>%
        rename(LGN = Lagnieu) %>%
        select(-2)

    model <- lm(formula = z_observed ~ x, data = data_WS_profile)
    residuals <- resid(model) # or model$residuals
    sd_residuals <- sd(residuals)

    data_WS_profile$Yu <- sd_residuals
    WS_profiles[[j]] <- data_WS_profile
    QdT_upstream[[j]] <- data_QdT_upstream

    names(WS_profiles)[j] <- names(QdT_upstream)[j] <- customized_names[j]

    all_data_calibration$WSE <- data_WS_profile

    save(all_data_calibration, file = file.path(path_results, paste0("data_WSE_Rhone_", customized_names[j], ".RData")))

    if (i == last(cases)) {
        all_data_calibration$WSE <- WS_profiles
        save(all_data_calibration, file = file.path(path_results, paste0("data_WSE_Rhone_all_cases.RData")))
        save(QdT_upstream, file = file.path(path_results, paste0("data_Q_upstream_Rhone_all_cases.RData")))
        save(riverbed, file = file.path(path_results, paste0("data_Thalweg_Rhone_all_cases.RData")))
        save(ZdQ_aval, file = file.path(path_results, paste0("data_Downstream_boundary_condition_Rhone_all_cases.RData")))
    }
    j <- j + 1
}
