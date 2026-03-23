rm(list = ls())
graphics.off()

library(dplyr)
library(fuzzyjoin)

setwd("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/data-raw/Synthetic_case/Rhone/Calibration/extraction_sim")

case_1 <- c("Rhone_synt_Ev_1")
case_2 <- c("Rhone_synt_Ev_2")

list_files_case_1 <- list.files(case_1, full.names = TRUE)
list_files_case_2 <- list.files(case_2, full.names = TRUE)

read_bief_file <- function(file) {
    # Lire toutes les lignes
    lines <- readLines(file)

    reach_numbers <- numeric()
    times <- numeric()
    Pks <- numeric()
    Cotes <- numeric()

    current_reach <- NA
    current_time <- NA

    for (line in lines) {
        # Extract reach number and time from comment lines
        if (grepl("^\\*", line)) {
            # Extract reach number
            match_reach <- regexec("Numéro du bief :\\s*(\\d+)", line, perl = TRUE)
            if (match_reach[[1]][1] != -1) {
                current_reach <- as.numeric(regmatches(line, match_reach)[[1]][2])
            }

            # Extract time
            match_time <- regexec("date \\(heures\\):\\s*([0-9.]+)", line, perl = TRUE)
            if (match_time[[1]][1] != -1) {
                current_time <- as.numeric(regmatches(line, match_time)[[1]][2])
            }
        }
        # Extract Pk and Cote values from data lines
        else if (grepl("[0-9.E+]+", line)) {
            parts <- strsplit(line, "\\s+")[[1]]
            parts <- parts[parts != ""] # Remove empty strings
            if (length(parts) >= 2) {
                Pk <- as.numeric(parts[1])
                Cote <- as.numeric(parts[2])

                # Append to vectors
                reach_numbers <- c(reach_numbers, current_reach)
                times <- c(times, current_time)
                Pks <- c(Pks, Pk)
                Cotes <- c(Cotes, Cote)
            }
        }
    }

    result_df <- data.frame(
        reach_number = reach_numbers,
        hours = times,
        KP = Pks,
        WSE = Cotes
    )

    return(result_df)
}

data_list_case_1 <- lapply(list_files_case_1, read_bief_file)
data_list_case_2 <- lapply(list_files_case_2, read_bief_file)


# Fusionner en un seul data.frame
data_case_1 <- do.call(rbind, data_list_case_1) %>% mutate(id_case = case_1)
data_case_2 <- do.call(rbind, data_list_case_2) %>% mutate(id_case = case_2)

WSE_synthetic_base <- rbind(data_case_1, data_case_2) %>%
    arrange(id_case, KP, reach_number, hours) %>%
    mutate(hours = hours * 3600) %>%
    rename(
        "KP" = "KP",
        "WSE_real_obs" = "WSE",
        "id_reach_CAL" = "reach_number",
        "time" = "hours"
    ) %>%
    filter(time == 10800)

# Perturb the simulated data to used it afterward as observations

# Get thalweg
##################################
## MAIN REACH UPSTREAM
##################################
file_bathy <- "/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/data-raw/Synthetic_case/Rhone/PamHYR_synthetic/_PAMHYR_/Rhone_synt_Ev_1/default-mage/net/"
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

#############################################################################
# Low uncertainty
#############################################################################

# Perturb the simulated data to used it afterward as observations
set.seed(2026) # reproducibility

WSE_synthetic_low_unc <- WSE_synthetic_base

WSE_synthetic_low_unc$Yu_WSE <- 0.025
WSE_synthetic_low_unc$WSE <-
    rnorm(
        n = nrow(WSE_synthetic_low_unc),
        mean = WSE_synthetic_low_unc$WSE_real_obs,
        sd = WSE_synthetic_low_unc$Yu_WSE
    )

library(ggplot2)

ggplot(
    WSE_synthetic_low_unc,
    aes(x = KP)
) +
    geom_point(aes(y = WSE_real_obs, col = "real")) +
    geom_point(aes(y = WSE, col = "perturbed")) +
    geom_errorbar(aes(ymin = WSE - 1.96 * Yu_WSE, ymax = WSE + 1.96 * Yu_WSE, col = "perturbed")) +
    theme_bw() +
    facet_wrap(~ id_case + id_reach_CAL)


WSE_synthetic_simplified <- WSE_synthetic_low_unc %>%
    difference_left_join(
        table_Z_thalweg,
        by = c("id_reach_CAL", "KP"),
        max_dist = 0.02
    ) %>%
    filter(id_reach_CAL.x == id_reach_CAL.y) %>%
    select(-c("id_reach_CAL.y", "KP.x")) %>%
    rename(
        id_reach_CAL = id_reach_CAL.x,
        KP = "KP.y"
    )

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
save(WSE_synthetic_simplified, file = "/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/processed_data/Synthetic_case/Rhone/Low_uncertainty/WSE_synthetic_Rhone.RData")
