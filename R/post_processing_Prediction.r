postprocess_prediction <- function(
    paths,
    mod_polynomials,
    Kmin_prior,
    Kflood_prior) {
    # Check if Results_Cooking.txt file exists
    check_calibration_done(path = paths$path_temp_plots)
    message("Processing: ", basename(dirname(paths$path_temp_plots)))

    # Count number of parameters different to FIX distribution
    n_param_Kmin_to_estimate <- sum(
        get_prior_distribution(Kmin_prior) != "FIX"
    )

    n_param_Kflood_to_estimate <- sum(
        get_prior_distribution(Kflood_prior) != "FIX"
    )
}
