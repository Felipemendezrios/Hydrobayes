rm(list = ls())
library(dplyr)

WSE_Durance <- read.table("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/data-raw/Durance_Thiercelin/Ref_data.csv", sep = ",", header = TRUE)

colnames(WSE_Durance) <- c("KP", "WSE")

WSE_Durance$Yu_WSE <- 0.05
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

# Add thalweg
file_bathy <- "/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/data-raw/Durance_Thiercelin/_PAMHYR_/PMRiver/default-mage/net/Reach_002.ST"
# Lagnieu reach
bathymetry <- readLines(file_bathy)
bathymetry <- bathymetry[-1]

# Remove headers and get xyz values

fields <- strsplit(trimws(bathymetry), "\\s+")
n_fields <- lengths(fields)

is_header <- n_fields >= 5

is_end <- sapply(fields, function(x) {
    length(x) >= 3 && all(x[1:3] == "999.9990")
})

is_xyz <- (n_fields %in% c(3, 4)) & !is_end

block_id <- cumsum(is_header) # each header starts a new block
block_id[is_end] <- NA # remove end marker lines

# keep only headers
header_lines <- fields[is_header]

# get first and before-last values
header_info <- t(sapply(header_lines, function(x) {
    first_val <- x[1]
    before_last <- x[length(x) - 1]
    c(profile_id = first_val, header_val = before_last)
}))

header_info <- as.data.frame(header_info, stringsAsFactors = FALSE)

xyz_list <- lapply(which(is_xyz), function(i) {
    nums <- as.numeric(fields[[i]][1:3])
    lbl <- if (length(fields[[i]]) >= 4) fields[[i]][4] else NA
    blk <- block_id[i] # which header/block this line belongs to
    c(nums, label = lbl, block = blk)
})

xyz_df <- as.data.frame(do.call(rbind, xyz_list), stringsAsFactors = FALSE)
colnames(xyz_df) <- c("X", "Y", "Z", "label", "block")
xyz_df$profile_id <- header_info$profile_id[as.numeric(xyz_df$block)]
xyz_df$header_val <- header_info$header_val[as.numeric(xyz_df$block)]

# optional: remove temporary block column
xyz_df$block <- NULL

# Convert numeric columns
xyz_df$X <- as.numeric(xyz_df$X)
xyz_df$Y <- as.numeric(xyz_df$Y)
xyz_df$Z <- as.numeric(xyz_df$Z)
xyz_df$profile_id <- as.numeric(xyz_df$profile_id)
xyz_df$header_val <- as.numeric(xyz_df$header_val)
xyz_df$id_reach <- "Reach_1"

Thalweg_data <- xyz_df %>%
    group_by(profile_id, id_reach) %>%
    summarise(
        KP = first(header_val),
        Z_thalweg = min(Z),
        .groups = "drop"
    )


save(Thalweg_data, file = "/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/processed_data/Durance_Thiercelin/Thalweg_Durance.RData")
