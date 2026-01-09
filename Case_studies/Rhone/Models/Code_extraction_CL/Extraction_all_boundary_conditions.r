rm(list = ls())
cat("\014")

setwd("/home/famendezrios/Documents/These/Cas_etude/Rhone/modelisation/Code_extraction_CL/All_boundary_conditions/")

library(data.table)
library(ggplot2)


extract_all_data <- function(
    data_source,
    path,
    export_name,
    out_dir,
    tz = "UTC") {
    if (data_source == "Hydroportail") {
        # Read data
        Data <- data.frame(fread(path, header = TRUE, skip = 1))
        Data <- Data[, c(2, 4)]
        colnames(Data) <- c("Date", "Q")
    } else if (data_source == "OSR") {
        Data <- data.frame(fread(path, header = TRUE, skip = 2))
        Data <- Data[, c(1, 2)]
        colnames(Data) <- c("Date", "Q")
    } else {
        stop("data_source must be either OSR or Hydroportail")
    }
    if (!inherits(Data$Date, "POSIXct")) {
        Data$Date <- as.POSIXct(
            Data$Date,
            format = "%d/%m/%Y %H:%M:%S",
            tz = tz
        )
    }
    # Ensure output directory exists
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    # Save
    out_file <- file.path(out_dir, paste0(export_name, ".RData"))
    save(Data, file = out_file)

    return(invisible(out_file))
}



metadata <- data.frame(
    export_name = c(
        "Ain_Pont_d_Ain",
        "Ain_Port_Galland",
        "Ain_Chazey_sur_Ain",
        "Bourbre_Tignieu_jameyzieu",
        "Rhone_Anthon",
        "Rhone_Jons",
        "Rhone_Lagnieu"
    ),
    path_file = c(
        "/home/famendezrios/Documents/These/Cas_etude/Rhone/Calibration/Condition_limite_BDOH_hydroportail/Ain_Pont_Ain_HYDROPORTAIL/Data_Pont_Ain_1980_2025.csv",
        "/home/famendezrios/Documents/These/Cas_etude/Rhone/Calibration/Condition_limite_BDOH_hydroportail/Ain_Port_Galland_OSR/PORT-GALLAND_DEB_1984_2025.txt",
        "/home/famendezrios/Documents/These/Cas_etude/Rhone/Calibration/Condition_limite_BDOH_hydroportail/Ain_Chazey_sur_Ain_HYDROPORTAIL/Chazey_data_1958_2025.csv",
        "/home/famendezrios/Documents/These/Cas_etude/Rhone/Calibration/Condition_limite_BDOH_hydroportail/Bourbre_Tignieu_jameyzieu_OSR/BOURBRE_DEB_1981_2025.txt",
        "/home/famendezrios/Documents/These/Cas_etude/Rhone/Calibration/Condition_limite_BDOH_hydroportail/Rhone_Anthon_HYDROPORTAIL/anthon_data_2023_2025.csv",
        "/home/famendezrios/Documents/These/Cas_etude/Rhone/Calibration/Condition_limite_BDOH_hydroportail/Rhone_Jons_OSR/JONS_DEB_1993_2025.txt",
        "/home/famendezrios/Documents/These/Cas_etude/Rhone/Calibration/Condition_limite_BDOH_hydroportail/Rhone_Lagnieu_OSR/LAGNIEU_DEB_1991_2025.txt"
    ),
    source = c(
        "Hydroportail",
        "OSR",
        "Hydroportail",
        "OSR",
        "Hydroportail",
        "OSR",
        "OSR"
    )
)


saved_files <- mapply(
    extract_all_data,
    data_source = metadata$source,
    path        = metadata$path_file,
    export_name = metadata$export_name,
    MoreArgs    = list(out_dir = getwd()),
    SIMPLIFY    = FALSE
)



save(metadata, file = file.path(getwd(), "metadata.RData"))
