rm(list = ls())
graphics.off()

library(dplyr)
library(ggplot2)
#################################################
# Extraction of parameter of the rectangular geometry :
# width and delta_z or depth
#################################################
get_ST_fortran <- function(
    path_ST, skip = 1) {
    bathy <- readLines(file.path(path_ST))
    bathy <- bathy[-c(1:skip)]

    # Remove headers and get xyz values

    fields <- strsplit(trimws(bathy), "\\s+")
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
    return(list(
        xyz_df = xyz_df,
        header_info = header_info
    ))
}

setwd("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/data-raw/Synthetic_case/Rhone/PamHyr/_PAMHYR_/Rhone_PCH_Hautes_eaux/default-mage/net/")

all_ST <- list.files()

Groupped_reaches <- list(
    Rhone = data.frame(
        ST = c("Reach_002.ST", "Reach_003.ST", "Reach_004.ST"),
        nb_reach = c(1, 2, 3)
    ),
    AIN = data.frame(
        ST = c("Reach_005.ST", "Reach_006.ST"),
        nb_reach = c(4, 5)
    )
)

results_geom_rect <- list()

# Check
for (group in names(Groupped_reaches)) {
    # group = names(Groupped_reaches)[2]
    files <- Groupped_reaches[[group]]$ST
    reaches <- Groupped_reaches[[group]]$nb_reach
    missing_files <- files[!files %in% all_ST]

    if (length(missing_files) > 0) {
        stop(paste("Les fichiers suivants n'existent pas :", paste(missing_files, collapse = ", ")))
    }

    results_st <- lapply(files, get_ST_fortran)

    all_xyz_df <- do.call(rbind, lapply(results_st, function(x) x$xyz_df))

    # Concaténation des header_info
    all_header_info <- do.call(rbind, lapply(results_st, function(x) x$header_info))


    save(all_xyz_df, file = paste0("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/data-raw/Synthetic_case/Rhone/rectangular_geometry/ST_read_XYZ_", group, ".RData"))
    save(all_header_info, file = paste0("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/data-raw/Synthetic_case/Rhone/rectangular_geometry/ST_read_headers", group, ".RData"))



    all_xyz_df_based <- data.frame(
        all_xyz_df %>%
            group_by(profile_id, KP) %>%
            mutate(
                # Point R (rg)
                X_rg = X[which(label == "rg")],
                Y_rg = Y[which(label == "rg")],
                Z_rg = Z[which(label == "rg")],

                # Point rd pour vecteur directeur
                X_rd = X[which(label == "rd")],
                Y_rd = Y[which(label == "rd")],
                Z_rd = Z[which(label == "rd")],

                # Vecteur directeur
                X_vect_dir = X_rd - X_rg,
                Y_vect_dir = Y_rd - Y_rg,
                Z_vect_dir = Z_rd - Z_rg,

                # Norme du vecteur
                norm_vect = sqrt(X_vect_dir^2 + Y_vect_dir^2 + Z_vect_dir^2),

                # Distance projetée sur la droite
                distance_proj = (
                    (X - X_rg) * X_vect_dir +
                        (Y - Y_rg) * Y_vect_dir +
                        (Z - Z_rg) * Z_vect_dir
                ) / norm_vect,
                # Find the row number of the first 'rg' in the group
                first_rg = which(label == "rg")[1],
                # Find the row number of the last 'rd' in the group
                last_rd = max(which(label == "rd")),
                # Create the new label
                new_label = case_when(
                    row_number() < first_rg ~ "flood_rg",
                    row_number() >= first_rg & row_number() <= last_rd ~ "main_channel",
                    row_number() > last_rd ~ "flood_rd",
                    TRUE ~ label
                )
            ) %>%
            select(-first_rg, -last_rd) # Remove temporary columns
    ) %>%
        ungroup()

    summarise_values <- all_xyz_df_based %>%
        group_by(profile_id, KP, new_label) %>%
        summarise(
            width = median(distance_proj, na.rm = TRUE)
        )

    # Calculate the median of width by new_label (across all data)
    median_widths <- summarise_values %>%
        group_by(new_label) %>%
        summarise(median_width = median(width))

    # Create the plot with vertical lines
    ggplot(summarise_values, aes(x = KP, y = width, color = new_label)) +
        geom_point(size = 3, alpha = 0.7) +
        # Add a vertical line for each median
        geom_hline(
            data = median_widths,
            aes(yintercept = median_width, color = new_label),
            linetype = "dashed",
            linewidth = 0.8,
            show.legend = FALSE
        ) +
        labs(
            title = "Variation Width from rg",
            y = "Width",
            x = "KP",
            color = "Label"
        ) +
        scale_color_manual(
            values = c("flood_rd" = "blue", "flood_rg" = "green", "main_channel" = "red")
        ) +
        theme_bw()

    ggsave(file = paste0("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/data-raw/Synthetic_case/Rhone/rectangular_geometry/variation_width_", group, ".png"))

    # Calculate the median of riverbed by new_label (across all data)

    all_xyz_df_based_riverbed <- all_xyz_df_based %>%
        group_by(profile_id, KP) %>%
        mutate(
            dz = Z - min(Z),
            new_label_grouped = case_when(
                new_label %in% c("flood_rd", "flood_rg") ~ "floodplain",
                TRUE ~ new_label
            )
        )

    summarise_riverbed <- all_xyz_df_based_riverbed %>%
        group_by(profile_id, KP, new_label_grouped) %>%
        summarise(
            riverbed = median(dz, na.rm = TRUE)
        )

    median_dz <- summarise_riverbed %>%
        group_by(new_label_grouped) %>%
        summarise(median_dz = median(riverbed))



    # Create the plot with vertical lines
    ggplot(summarise_riverbed, aes(y = riverbed, x = KP, color = new_label_grouped)) +
        geom_point(size = 3, alpha = 0.7) +
        # Add a vertical line for each median
        geom_hline(
            data = median_dz,
            aes(yintercept = median_dz, color = new_label_grouped),
            linetype = "dashed",
            linewidth = 0.8,
            show.legend = FALSE
        ) +
        labs(
            title = "Variation riverbed from thalweg",
            y = "Elevation",
            x = "KP",
            color = "Label"
        ) +
        scale_color_manual(
            values = c("floodplain" = "brown", "main_channel" = "yellow")
        ) +
        theme_bw()

    ggsave(file = paste0("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/data-raw/Synthetic_case/Rhone/rectangular_geometry/variation_dz_", group, ".png"))

    results_geom_rect[[group]] <- list(
        width = data.frame(median_widths),
        riverbed = median_dz
    )
}


#################################################
# Extraction of the slope:
# Extract from a simulation of the real model with a high discharge to obtain the energy slope which could be assumed as similar to the slope of the friction and slope of the bathymetry (manning assumption)
#################################################

setwd("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/data-raw/Synthetic_case/Rhone/PamHyr/_PAMHYR_/Rhone_PCH_Hautes_eaux/default-mage/")

library(readr)

# Read file
lines <- readLines("Mage_fin.ini")

# Keep only numeric data lines
data_lines <- lines[grepl("^\\s*\\d", lines)]

# Split by multiple spaces
split_lines <- strsplit(trimws(data_lines), "\\s+")

# Function to fix rows where ! or # is separated
fix_line <- function(x) {
    if (length(x) == 25) {
        # Find position of ! or #
        pos <- which(x %in% c("!", "#"))

        if (length(pos) == 1 && pos > 1) {
            # Attach symbol to previous value
            x[pos - 1] <- paste0(x[pos - 1], x[pos])
            x <- x[-pos]
        }
    }

    return(x)
}

# Apply correction
split_fixed <- lapply(split_lines, fix_line)

# Check column count
table(sapply(split_fixed, length))

df_sim_high_Q <- as.data.frame(do.call(rbind, split_fixed),
    stringsAsFactors = FALSE
)

# Assign column names manually
colnames(df_sim_high_Q) <- c(
    "IB", "IS", "Debit", "Cote", "Pm", "Z_Fond", "Z_moyen", "Z_Berge",
    "H", "Largeur", "Fr", "V", "Q_Lat", "Qdev", "K_min", "K_maj",
    "Surface", "Perimetre", "Q_mineur", "Q_moyen",
    "V_mineur", "V_moyen", "Debitance", "Boussinesq"
)

# Remove symbol and convert to numeric
df_sim_high_Q$H <- as.numeric(gsub("[#!]", "", df_sim_high_Q$H))

# Convert everything to numeric
df_sim_high_Q[] <- lapply(df_sim_high_Q, as.numeric)



# Check
for (group in names(Groupped_reaches)) {
    # group <- names(Groupped_reaches)[1]
    reaches <- Groupped_reaches[[group]]$nb_reach

    results_geom_rect[[group]]$So <- abs(as.vector(coef(lm(Cote ~ Pm, data = df_sim_high_Q %>% filter(IB %in% reaches)))["Pm"]))
}



save(results_geom_rect, file = "/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/data-raw/Synthetic_case/Rhone/rectangular_geometry/values_width_dz.RData")
