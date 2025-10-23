Cross_sections_interpolation_and_export_rectangular_channel <- function(
    So,
    Pk,
    positions_banks,
    borders_heigh,
    path_export) {
    PK <- sort(Pk)
    ST <- data.frame(PK = PK)
    ST$z <- NA

    j <- 1
    for (i in PK) {
        ST$z[j] <- So * (PK[length(PK)] - i)
        j <- j + 1
    }

    con <- file(path_export, open = "w")

    format_str_entete <- "%6d%6d%6d%6d%13.8f  %-20s"
    format_str_xyz <- "%10.8f %10.8f %10.8f %-3s"

    # Use a for loop to write lines
    for (i in 1:nrow(ST)) {
        # Appliquer ligne par ligne
        entete_line <- sprintf(
            format_str_entete,
            as.integer(i),
            as.integer(0),
            as.integer(0),
            as.integer(6),
            as.numeric(ST$PK[i]),
            as.character("")
        )
        # Écrire dans un fichier
        writeLines(entete_line, con)

        XYZ_df <- data.frame(
            x = c(rep(ST$PK[i], 6), 999.999),
            y = c(rep(positions_banks[1], 3), rep(positions_banks[2], 3), 999.999),
            z = c(
                rep(borders_heigh + ST$z[i], 2),
                rep(ST$z[i], 2),
                rep(borders_heigh + ST$z[i], 2),
                999.999
            ),
            char = c("", "rg", rep("", 2), "rd", "", "")
        )

        # Générer toutes les lignes formatées
        XYZ <- apply(XYZ_df, 1, function(row) {
            sprintf(
                format_str_xyz,
                as.numeric(row["x"]),
                as.numeric(row["y"]),
                as.numeric(row["z"]),
                as.character(row["char"])
            )
        })
        # Écrire dans un fichier
        writeLines(XYZ, con)
    }

    # Close the connection
    close(con)
}



# assign_reach <- function(specific_points_Model, covariate_value) {
#     reaches <- rep(NA, length(covariate_value))

#     if (nrow(specific_points_Model == 1)) {
#         all_set <- dplyr::between(covariate_value, specific_points_Model$KP_start, specific_points_Model$KP_end)
#         if (!all(all_set)) stop("specific_points_Model must contain all covariate_value")
#         reaches[which(all_set)] <- specific_points_Model$reach[1]
#     } else {
#         duplicated_numbers <- covariate_value[duplicated(covariate_value)]

#         if ((nrow(specific_points_Model) - 1) != length(duplicated_numbers)) stop("Number of rows of specific_points_Model minus 1 must coincide with the number of repeated values in covariate_value")

#         for (i in 1:nrow(specific_points_Model)) {
#             idx_position_reaches <- which(dplyr::between(covariate_value, specific_points_Model$KP_start[i], specific_points_Model$KP_end[i]))

#             idx_duplicated_values <- which(covariate_value[idx_position_reaches] %in% duplicated_numbers)

#             if (length(idx_duplicated_values) != 0) {
#                 if (length(idx_duplicated_values) != 2 & length(idx_duplicated_values) != 4) stop("Duplicated values are not found, it must be two or four")

#                 # First duplicated values
#                 if (i == 1 & length(idx_duplicated_values) == 2) {
#                     reaches[idx_position_reaches] <- specific_points_Model$reach[i]

#                     reaches[idx_duplicated_values[2]] <- specific_points_Model$reach[i + 1]
#                 } else if (length(idx_duplicated_values) == 4) {
#                     idx_position_reaches <- idx_position_reaches[-idx_duplicated_values[c(1, 2)]]
#                     reaches[idx_position_reaches] <- specific_points_Model$reach[i]
#                     reaches[idx_duplicated_values[4]] <- specific_points_Model$reach[i + 1]
#                 } else if (length(idx_duplicated_values) == 2) {
#                     idx_position_reaches <- idx_position_reaches[-idx_duplicated_values[c(1, 2)]]
#                     reaches[idx_position_reaches] <- specific_points_Model$reach[i]
#                 }
#             } else {
#                 reaches[idx_position_reaches] <- specific_points_Model$reach[i]
#             }
#         }
#     }
#     grid_covariant_discretized <- data.frame(Reach = reaches, Covariate = covariate_value)
#     return(grid_covariant_discretized)
# }


assign_reach_from_a_grid <- function(metadata_specific_points, grid) {
    reaches <- rep(NA, length(grid))

    if (nrow(metadata_specific_points) == 1) {
        all_set <- dplyr::between(grid, metadata_specific_points$KP_start, metadata_specific_points$KP_end)
        if (!all(all_set)) stop("metadata_specific_points must contain all grid")
        reaches[which(all_set)] <- metadata_specific_points$reach[1]
    } else {
        duplicated_numbers <- grid[duplicated(grid)]
        if ((nrow(metadata_specific_points) - 1) != length(duplicated_numbers)) stop("Number of rows of specific_points_Model minus 1 must coincide with the number of repeated values in the grid")

        for (i in 1:nrow(metadata_specific_points)) {
            idx_position_reaches <- which(dplyr::between(grid, metadata_specific_points$KP_start[i], metadata_specific_points$KP_end[i]))

            idx_duplicated_values <- which(grid[idx_position_reaches] %in% duplicated_numbers)

            if (length(idx_duplicated_values) != 0) {
                if (length(idx_duplicated_values) != 2 & length(idx_duplicated_values) != 4) stop("Duplicated values are not found, the number of duplicated values must be two or four")

                # First duplicated values
                if (i == 1 & length(idx_duplicated_values) == 2) {
                    reaches[idx_position_reaches] <- metadata_specific_points$reach[i]

                    reaches[idx_position_reaches[idx_duplicated_values[2]]] <- metadata_specific_points$reach[i + 1]
                } else if (length(idx_duplicated_values) == 4) {
                    idx_position_reaches <- idx_position_reaches[-idx_duplicated_values[c(1, 2)]]

                    reaches[idx_position_reaches] <- metadata_specific_points$reach[i]

                    reaches[idx_position_reaches[(idx_duplicated_values[4] - 2)]] <- metadata_specific_points$reach[i + 1] # -2 because two values have been removed previously
                } else if (length(idx_duplicated_values) == 2) {
                    idx_position_reaches <- idx_position_reaches[-idx_duplicated_values[c(1, 2)]]
                    reaches[idx_position_reaches] <- metadata_specific_points$reach[i]
                }
            } else {
                stop("It should never get here. Something is wrong")
                reaches[idx_position_reaches] <- metadata_specific_points$reach[i]
            }
        }
    }

    if (any(is.na(reaches))) {
        return(stop("Some values in the vector reaches have not been assigned"))
    }

    return(reaches)
}
