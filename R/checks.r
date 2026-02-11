check_mage_folder <- function(name_folder) {
    if (!dir.exists(name_folder)) stop("Check the path to the Mage main folder. It does not exist")
}



check_simulation_time <- function(
    MAGE_main_folder,
    mage_projet_name,
    Observations,
    event) {
    PARFile <- readLines(file.path(MAGE_main_folder, event, paste0(mage_projet_name, ".PAR")), n = 7)
    end_par_file <- PARFile[3]
    timestep_bin_par_file <- PARFile[7]

    # Get the mage time in seconds
    duration <- sub("^final_time\\s+", "", end_par_file)
    parts <- as.numeric(strsplit(duration, ":")[[1]])

    end_time_sim <- parts[1] * 86400 + # days
        parts[2] * 3600 + # hours
        parts[3] * 60 + # minutes
        parts[4] # seconds

    if (nrow(Observations %>% filter(Observations$time > end_time_sim)) != 0) stop("There are some observed data beyond the simulation time period")
}
