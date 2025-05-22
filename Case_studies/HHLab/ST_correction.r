So <- 1.04e-3

PK <- c(
    0,
    0.06,
    seq(0.5040, 16.5040, by = 0.1),
    0.65,
    1.65,
    3.65,
    5.65,
    7.65,
    9.65,
    11.65,
    15.65
)
PK <- sort(PK)
ST <- data.frame(PK = PK)
ST$z <- NA
j <- 1
for (i in PK) {
    ST$z[j] <- So * (PK[length(PK)] - i)
    j <- j + 1
}

fin_line <- c(rep("999.9990", 3))

b <- c(0, 2)
z_paroi <- 1

file_path <- "/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/Case_studies/HHLab//MAGE_PaMHYR2/_PAMHYR_/smooth_bed/default-mage/net/Reach_002.ST"

con <- file(file_path, open = "w")

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
        y = c(rep(0, 3), rep(2, 3), 999.999),
        z = c(
            rep(1 + ST$z[i], 2),
            rep(ST$z[i], 2),
            rep(1 + ST$z[i], 2),
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
