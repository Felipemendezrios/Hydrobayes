# DIC
cat("\014")
rm(list = ls())

library(ggplot2)

workspace <- "/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/Case_studies/HHLab"
setwd(workspace)

# Here, define calibration type
# 'Calibration_time_series'
# 'Calibration_water_surface_profiles'
path_results <- "Calibration_water_surface_profiles"

# Choose DIC estimation : DIC1, DIC2, DIC3, D(maxpost)
DIC_criterion <- "DIC3"

all_files <- list.files(
    file.path(
        path_results
    ),
    recursive = TRUE
)

DIC_path_Files <- all_files[grep(all_files, pattern = "Results_DIC.txt")]
if (length(DIC_path_Files) == 0) stop("Any Results_DIC.txt file was found")
DIC_results <- c()

for (i in 1:length(DIC_path_Files)) {
    DIC_by_degree <- read.table(file.path(path_results, DIC_path_Files[i]), col.names = c("Criteria", "Value"))
    # Match the criterion chosen
    DIC_by_degree <- DIC_by_degree[which(DIC_by_degree[, 1] == DIC_criterion), ]
    # Assign the polynomial degree
    extraction <- strsplit(DIC_path_Files[i], "/")[[1]][1]
    # Extraire le chiffre après le underscore
    chiffre <- sub(".*_(\\d+)", "\\1", extraction)

    # Convertir en numérique (optionnel)
    DIC_by_degree$Degree <- as.numeric(chiffre)

    DIC_results <- rbind(DIC_results, DIC_by_degree)
}
min_local <- DIC_results[which.min(DIC_results$Value), ]

DIC_plot <- ggplot(DIC_results, aes(x = factor(Degree), y = Value, col = factor(Criteria))) +
    geom_point(size = 3) +
    geom_point(data = min_local, aes(x = , factor(Degree), y = Value), col = "blue", size = 3) +
    annotate("segment",
        x = factor(min_local$Degree), y = min_local$Value * 1.005, xend = factor(min_local$Degree), yend = min_local$Value * 1.002,
        linewidth = 2, linejoin = "mitre",
        arrow = arrow(type = "closed", length = unit(0.01, "npc"))
    ) +
    theme_bw() +
    labs(
        x = "Polynomial degree",
        y = "Value",
        col = "Criterion",
        title = "DIC criterion :",
        subtitle = "Comparison of several polynomial degrees"
    ) +
    theme(
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)
    )


ggsave(
    filename = file.path(path_results, "DIC_estimation.png"),
    plot = DIC_plot,
    dpi = 300, width = 1400, height = 1400,
    units = "px"
)
