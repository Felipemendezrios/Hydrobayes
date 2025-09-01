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
