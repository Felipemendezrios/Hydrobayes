rm(list = ls())


WSE_Durance <- read.table("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/data-raw/Durance_Thiercelin/Ref_data.csv", sep = ";", header = TRUE)

colnames(WSE_Durance) <- c("KP", "WSE")

WSE_Durance$Yu_WSE <- 0.01
WSE_Durance$id_reach_CAL <- 1
library(ggplot2)
ggplot(
    WSE_Durance,
    aes(
        x = KP,
        y = WSE,
        ymin = WSE - 2 * Yu_WSE,
        ymax = WSE + 2 * Yu_WSE
    )
) +
    geom_point(size = 0.5) +
    geom_errorbar() +
    theme_bw()

WSE_Durance$time <- 43200
WSE_Durance$id_case <- 1

save(WSE_Durance, file = "/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/processed_data/Durance_Thiercelin/WSE_Durance.RData")
