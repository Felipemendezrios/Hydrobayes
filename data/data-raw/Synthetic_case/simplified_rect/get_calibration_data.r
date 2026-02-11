rm(list = ls())
graphics.off()

library(dplyr)

setwd("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/data-raw/Synthetic_case/simplified_rect/")

case_1 <- c("3_2")
case_2 <- c("5_1")

list_files_case_1 <- list.files(case_1, full.names = TRUE)
list_files_case_2 <- list.files(case_2, full.names = TRUE)

read_bief_file <- function(file) {
    # Lire toutes les lignes
    lines <- readLines(file)

    # Extraire numéro du bief
    bief <- as.numeric(
        sub(".*Numéro du bief *: *([0-9]+).*", "\\1", lines[1])
    )

    # Extraire l'heure
    heure <- as.numeric(
        sub(".*date \\(heures\\): *([0-9\\.]+).*", "\\1", lines[1])
    )

    # Lire les données numériques (à partir de la ligne 3)
    data <- read.table(
        file,
        skip = 2,
        col.names = c("Pk", "Cotes")
    )

    # Ajouter les colonnes
    data$bief <- bief
    data$heure <- heure

    return(data)
}

data_list_case_1 <- lapply(list_files_case_1, read_bief_file)
data_list_case_2 <- lapply(list_files_case_2, read_bief_file)


# Fusionner en un seul data.frame
data_case_1 <- do.call(rbind, data_list_case_1) %>% mutate(id_case = case_1)
data_case_2 <- do.call(rbind, data_list_case_2) %>% mutate(id_case = case_2)

WSE_synthetic_temp <- rbind(data_case_1, data_case_2) %>%
    arrange(id_case, Pk, bief, heure) %>%
    mutate(heure = heure * 3600) %>%
    rename(
        "KP" = "Pk",
        "WSE_real_obs" = "Cotes",
        "id_reach_CAL" = "bief",
        "time" = "heure"
    ) %>%
    filter(time == 10800)

# Perturb the simulated data to used it afterward as observations
set.seed(2026) # reproducibility
WSE_synthetic_temp$Yu_WSE <- 0.001
WSE_synthetic_temp$WSE <-
    rnorm(
        n = nrow(WSE_synthetic_temp),
        mean = WSE_synthetic_temp$WSE_real_obs,
        sd = WSE_synthetic_temp$Yu_WSE
    )



# Get thalweg
##################################
## MAIN REACH UPSTREAM
##################################
file_bathy <- "/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/scripts/Case_studies/Synthetic_case/Simplified_ks_rectangular_MC/model_mage/5_1/net/"
Bathy_ST1 <- readLines(file.path(file_bathy, "Reach_001.ST"))
Bathy_ST1 <- Bathy_ST1[-1]

# Remove headers and get xyz values
fields <- strsplit(trimws(Bathy_ST1), "\\s+")
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
    c(profile_id = first_val, KP = before_last)
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
xyz_df$KP <- header_info$KP[as.numeric(xyz_df$block)]

# optional: remove temporary block column
xyz_df$block <- NULL

# Convert numeric columns
xyz_df$X <- as.numeric(xyz_df$X)
xyz_df$Y <- as.numeric(xyz_df$Y)
xyz_df$Z <- as.numeric(xyz_df$Z)
xyz_df$profile_id <- as.numeric(xyz_df$profile_id)
xyz_df$KP <- as.numeric(xyz_df$KP)

table_Z_thalweg <- data.frame(
    xyz_df %>%
        group_by(KP) %>%
        filter(Z == min(Z)) %>%
        summarise(
            Z_thalweg = first(Z)
        ) %>%
        mutate(id_reach_CAL = 1)
)

##################################
## TRIBUTARY
##################################
Bathy_ST2 <- readLines(file.path(file_bathy, "Reach_002.ST"))
Bathy_ST2 <- Bathy_ST2[-1]

# Remove headers and get xyz values
fields <- strsplit(trimws(Bathy_ST2), "\\s+")
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
    c(profile_id = first_val, KP = before_last)
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
xyz_df$KP <- header_info$KP[as.numeric(xyz_df$block)]

# optional: remove temporary block column
xyz_df$block <- NULL

# Convert numeric columns
xyz_df$X <- as.numeric(xyz_df$X)
xyz_df$Y <- as.numeric(xyz_df$Y)
xyz_df$Z <- as.numeric(xyz_df$Z)
xyz_df$profile_id <- as.numeric(xyz_df$profile_id)
xyz_df$KP <- as.numeric(xyz_df$KP)

table_Z_thalweg <- rbind(table_Z_thalweg, data.frame(
    xyz_df %>%
        group_by(KP) %>%
        filter(Z == min(Z)) %>%
        summarise(
            Z_thalweg = first(Z)
        ) %>%
        mutate(id_reach_CAL = 2)
))

##################################
## MAIN CHANNEL DOWNSTREAM
##################################
Bathy_ST3 <- readLines(file.path(file_bathy, "Reach_003.ST"))
Bathy_ST3 <- Bathy_ST3[-1]

# Remove headers and get xyz values
fields <- strsplit(trimws(Bathy_ST3), "\\s+")
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
    c(profile_id = first_val, KP = before_last)
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
xyz_df$KP <- header_info$KP[as.numeric(xyz_df$block)]

# optional: remove temporary block column
xyz_df$block <- NULL

# Convert numeric columns
xyz_df$X <- as.numeric(xyz_df$X)
xyz_df$Y <- as.numeric(xyz_df$Y)
xyz_df$Z <- as.numeric(xyz_df$Z)
xyz_df$profile_id <- as.numeric(xyz_df$profile_id)
xyz_df$KP <- as.numeric(xyz_df$KP)

table_Z_thalweg <- rbind(table_Z_thalweg, data.frame(
    xyz_df %>%
        group_by(KP) %>%
        filter(Z == min(Z)) %>%
        summarise(
            Z_thalweg = first(Z)
        ) %>%
        mutate(id_reach_CAL = 3)
))

########################
# Add the thalweg to the data


WSE_synthetic_simplified <- WSE_synthetic_temp %>% left_join(table_Z_thalweg, by = c("id_reach_CAL", "KP"))


library(ggplot2)
wse_obs_sim <- ggplot(WSE_synthetic_simplified, aes(x = KP, y = WSE)) +
    # Calibration data
    geom_point(aes(color = "Calibration Data")) +
    # Real observations
    geom_point(aes(x = KP, y = WSE_real_obs, color = "real obs")) +
    # Error bars
    geom_errorbar(aes(
        ymin = WSE - 1.96 * Yu_WSE, ymax = WSE + 1.96 * Yu_WSE,
        color = "Calibration Data"
    )) +
    scale_color_manual(
        values = c(
            "Calibration Data" = "#1b9e77",
            "real obs" = "black",
            "riverbed" = "brown"
        ),
        breaks = c("real obs", "Calibration Data", "riverbed"),
        labels = c(
            "Real observations",
            "Calibration Data",
            "Riverbed"
        )
    ) +
    facet_wrap(~ id_case + id_reach_CAL,
        scales = "free_x",
        labeller = labeller(
            id_case = c(
                "100_70" = "MC 100 TR 70",
                "120_40" = "MC 120 TR 40"
            ),
            reach = c(
                "1" = "Upstream main reach",
                "2" = "Upstream tributary reach",
                "3" = "Downstream main reach"
            )
        )
    ) +
    labs(x = "Kilometer point (km)", y = "Water surface elevation (meters)", color = "Reach") +
    theme_bw()
wse_obs_sim

wse_obs_sim +
    # Riverbed
    geom_line(aes(y = Z_thalweg, color = "riverbed"))

# Save data frame
save(WSE_synthetic_simplified, file = "/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/processed_data/Synthetic_case/WSE_synthetic_simplified_rectangle.RData")
