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

check_type_pred <- function(type) {
    if (!any(type %in% c("ZdX", "QdXT", "VdXT", "ZdT"))) stop(paste0(type, " must be either ZdX, or QdXT, or VdXT, or ZdT"))
}

check_dX <- function(df_event) {
    if (!identical(df_event$tmin, df_event$tmax)) stop("tmin and tmax must be identical")
}

check_dT <- function(df_event) {
    if (!identical(df_event$xmin, df_event$xmax)) stop("xmin and xmax must be identical")
}

check_dXT <- function(df_event) {
    if (!identical(df_event$xmin, df_event$xmax) | !identical(df_event$tmin, df_event$tmax)) stop("Both xmin, xmax and tmin, tmax must be identical")
}

check_calibration_case <- function(path) {
    if (!dir.exists(path)) stop(paste0("Calibration case is not performed yet: ", basename(path)))

    if (!file.exists(file.path(
        path,
        "BaM",
        "Results_Cooking.txt"
    ))) {
        stop("Results_Cooking.txt file is missing. Calibration is not performed yet")
    }
}

check_calibration_done <- function(path) {
    if (!file.exists(file.path(path, "Results_Cooking.txt"))) stop("MCMC is still running or calculation is not going to the end. Verify if calibration is already finished or verify that calibration has not error messages. Please put final_results = FALSE as input")
}

check_suffix_pred <- function(file) {
    # Check if all lists are empty
    if (length(file) == 0) stop(paste0("No files found for prediction of : ", basename(file)))
}
