# ============================================================
# POSTPROCESS CALIBRATION WORKFLOW
# ============================================================


# ============================================================
# 1. LOAD EXPERIMENT
# ============================================================

load_experiment <- function(file_main_path, cal_case, path_experiment, all_events) {
    path_input <- file.path(file_main_path, "Experiments_Input_Data", cal_case)

    if (!file.exists(path_input)) {
        stop("Experiment input file does not exist: ", cal_case)
    }

    source(path_input)

    path_polynomial <- file.path(
        path_experiment,
        sub("\\.r$", "", cal_case)
    )

    path_temp_plots <- file.path(path_polynomial, "BaM")

    path_post <- file.path(path_temp_plots, "post_traitement")

    path_post_data <- file.path(path_post, "RData")

    dir.create(path_post, showWarnings = FALSE)
    dir.create(path_post_data, showWarnings = FALSE)

    path_model_mage <- paste0(
        path_polynomial,
        "/model_mage/",
        all_events,
        "/"
    )

    return(list(
        path_temp_plots = path_temp_plots,
        path_post = path_post,
        path_post_data = path_post_data,
        path_model_mage = path_model_mage
    ))
}


# ============================================================
# 2. MCMC DIAGNOSTICS
# ============================================================

plot_mcmc_diagnostics <- function(
    path_temp_plots,
    path_post,
    final_calibration) {
    file <- if (final_calibration) {
        "Results_Cooking.txt"
    } else {
        "Results_MCMC.txt"
    }


    fullpath <- file.path(path_temp_plots, file)


    if (!file.exists(fullpath)) {
        stop("MCMC file not found: ", fullpath)
    }


    mcmc <- RBaM::readMCMC(fullpath)


    trace <- patchwork::wrap_plots(
        RBaM::tracePlot(mcmc),
        ncol = 3
    )


    ggplot2::ggsave(
        file.path(
            path_post,
            paste0("MCMC_", tools::file_path_sans_ext(file), ".png")
        ),
        trace,
        width = 20,
        height = 20,
        units = "cm"
    )


    density <- patchwork::wrap_plots(
        RBaM::densityPlot(mcmc),
        ncol = 3
    )


    ggplot2::ggsave(
        file.path(
            path_post,
            paste0("density_", tools::file_path_sans_ext(file), ".png")
        ),
        density,
        width = 20,
        height = 20,
        units = "cm"
    )


    return(mcmc)
}



# ============================================================
# 3. EXTRACT MAP
# ============================================================

extract_MAP_parameters <- function(
    path_temp_plots,
    final_calibration,
    Kmin_prior,
    Kflood_prior) {
    Fix_Kmin <- sum(
        sapply(Kmin_prior, function(x) x$prior$dist == "FIX")
    )

    Fix_Kflood <- sum(
        sapply(Kflood_prior, function(x) x$prior$dist == "FIX")
    )

    if (final_calibration) {
        summary <- read.table(
            file.path(
                path_temp_plots,
                "Results_Summary.txt"
            ),
            header = TRUE
        )

        summary_zoom <- summary[c(11, 16), ]

        MAP <- as.numeric(
            summary_zoom[
                2,
                1:(length(Kmin_prior) - Fix_Kmin + length(Kflood_prior) - Fix_Kflood)
            ]
        )
    } else {
        res <- read.table(
            file.path(
                path_temp_plots,
                "Results_MCMC.txt"
            ),
            header = TRUE
        )


        MAP <- as.numeric(
            res[
                which.max(res$LogPost),
                1:(length(Kmin_prior) - Fix_Kmin + length(Kflood_prior) - Fix_Kflood)
            ]
        )
    }

    return(MAP)
}



# ============================================================
# 4. COMPUTE K
# ============================================================

compute_K <- function(
    Z_MatrixK,
    mcmc,
    Kmin_prior,
    Kflood_prior,
    K_SU,
    MAP,
    n_param_Kmin_to_estimate,
    n_param_Kflood_to_estimate,
    do_main_channel) {
    SU <- do.call(
        rbind,
        lapply(names(K_SU), function(river) {
            data.frame(
                KP = K_SU[[river]][[1]]$KP,
                reach = K_SU[[river]][[1]]$reach,
                id_river = river
            )
        })
    )
    if (do_main_channel) {
        if (n_param_Kmin_to_estimate == 0) {
            return(stop("Any parameter is estimated in the main channel. Please verify the parameters in the main channel"))
        }
    } else {
        if (n_param_Kflood_to_estimate == 0) {
            return(
                stop("Any parameter is estimated in the floodplain. Please verify the parameters in the floodplain")
            )
        }
    }

    res <- K_plot(
        matrix_spatialisation = Z_MatrixK,
        mcmc = mcmc,
        n_param_Kmin = length(Kmin_prior),
        n_param_Kflood = length(Kflood_prior),
        SU_reaches = SU,
        MAP_param_vector = MAP,
        main_channel = do_main_channel
    )


    return(list(
        df_MAP = res[[1]],
        plot = res[[2]],
        envelope = res[[3]]
    ))
}



# ============================================================
# 5. COMPUTE RESIDUALS
# ============================================================

compute_residuals <- function(
    path_temp_plots) {
    CalData <- read.table(
        file.path(
            path_temp_plots,
            "CalibrationData.txt"
        ),
        header = TRUE
    )

    residuals <- read.table(
        file.path(
            path_temp_plots,
            "Results_Residuals.txt"
        ),
        header = TRUE
    )

    # Colonnes X_obs
    X_obs_cols <- grep("^X[0-9]+_obs$", colnames(residuals), value = TRUE)

    # Colonnes Y_sim
    Y_sim_cols <- grep("^Y[0-9]+_sim$", colnames(residuals), value = TRUE)

    # Colonnes Y_res
    Y_res_cols <- grep("^Y[0-9]+_res$", colnames(residuals), value = TRUE)

    # Colonnes Y_stdres
    Y_stdres_cols <- grep("^Y[0-9]+_stdres$", colnames(residuals), value = TRUE)

    res <- data.frame(
        residuals[, c(
            X_obs_cols,
            Y_sim_cols,
            Y_res_cols,
            Y_stdres_cols
        )],
        CalData[, c(
            grep("^Yu_", colnames(CalData), value = TRUE)
        )]
    )
    return(res)
}


# ============================================================
# 6. PLOT K with referenced values
# ============================================================

plot_K_and_ref <- function(
    K_results,
    K_segment_layer,
    path_post,
    path_post_data) {
    final_plot <-
        K_results$plot +

        K_segment_layer +

        ggplot2::scale_linetype_manual(
            name = "Reference values",
            values = c(
                "mean" = "dashed",
                "min" = "dashed",
                "max" = "dashed"
            ),
            labels = c("synthetic data")
        )
    return(final_plot)
}



# ============================================================
# MAIN FUNCTION
# ============================================================

postprocess_calibration <- function(
    all_cal_case,
    file_main_path,
    path_experiment,
    all_events,
    final_calibration = TRUE,
    list_Kmin_prior,
    list_Kflood_prior,
    list_Kmin_SU,
    list_Kflood_SU,
    list_Z_MatrixK,
    list_Z_MatrixKflood,
    Kmin_segment_layer = NULL,
    Kflood_segment_layer = NULL) {
    for (id_cal_case in 1:length(all_cal_case)) {
        message("Processing: ", all_cal_case[[id_cal_case]])

        # Load experiments
        exp <- load_experiment(
            file_main_path = file_main_path,
            cal_case = all_cal_case[[id_cal_case]],
            path_experiment = path_experiment,
            all_events = all_events
        )

        # MCMC
        mcmc <- plot_mcmc_diagnostics(
            path_temp_plots = exp$path_temp_plots,
            path_post = exp$path_post,
            final_calibration = final_calibration
        )

        # Count number of parameters different to FIX distribution

        n_param_Kmin_to_estimate <- sum(
            sapply(list_Kmin_prior[[id_cal_case]], function(x) x$prior$dist != "FIX")
        )

        n_param_Kflood_to_estimate <- sum(
            sapply(list_Kflood_prior[[id_cal_case]], function(x) x$prior$dist != "FIX")
        )

        # MAP
        MAP <- extract_MAP_parameters(
            path_temp_plots = exp$path_temp_plots,
            final_calibration = final_calibration,
            Kmin_prior = list_Kmin_prior[[id_cal_case]],
            Kflood_prior = list_Kflood_prior[[id_cal_case]]
        )

        # Kmin
        if (n_param_Kmin_to_estimate) {
            Kmin <- compute_K(
                Z_MatrixK = list_Z_MatrixKmin[[id_cal_case]],
                mcmc = mcmc,
                Kmin_prior = list_Kmin_prior[[id_cal_case]],
                Kflood_prior = list_Kflood_prior[[id_cal_case]],
                K_SU = list_Kmin_SU[[id_cal_case]],
                MAP = MAP,
                do_main_channel = TRUE,
                n_param_Kmin_to_estimate = n_param_Kmin_to_estimate,
                n_param_Kflood_to_estimate = n_param_Kflood_to_estimate
            )

            if (!is.null(Kmin_segment_layer)) {
                final_plot_Kmin <- plot_K_and_ref(
                    K_results = Kmin,
                    Kmin_segment_layer = Kmin_segment_layer,
                    path_post = exp$path_post,
                    path_post_data = exp$path_post_data
                )
            } else {
                final_plot_Kmin <- Kmin$plot
            }

            ggplot2::ggsave(
                file.path(exp$path_post, "Kmin.png"),
                final_plot_Kmin,
                width = 20,
                height = 20,
                units = "cm"
            )

            save(
                final_plot_Kmin,
                file = file.path(
                    exp$path_post_data,
                    "Plot_friction_estimation_plot_Kmin_plot.RData"
                )
            )
        } else {
            Kmin <- NULL
            warning("Any parameter is estimated in the main channel. Kmin returned is NULL")
        }
        # Kflood
        if (n_param_Kflood_to_estimate != 0) {
            Kflood <- compute_K(
                Z_MatrixK = list_Z_MatrixKflood[[id_cal_case]],
                mcmc = mcmc,
                Kmin_prior = list_Kmin_prior[[id_cal_case]],
                Kflood_prior = list_Kflood_prior[[id_cal_case]],
                K_SU = list_Kflood_SU[[id_cal_case]],
                MAP = MAP,
                do_main_channel = FALSE,
                n_param_Kmin_to_estimate = n_param_Kmin_to_estimate,
                n_param_Kflood_to_estimate = n_param_Kflood_to_estimate
            )

            if (!is.null(Kflood_segment_layer)) {
                final_plot_Kflood <- plot_K_and_ref(
                    K_results = Kflood,
                    K_segment_layer = Kflood_segment_layer,
                    path_post = exp$path_post,
                    path_post_data = exp$path_post_data
                )
            } else {
                final_plot_Kflood <- Kflood$plot
            }

            ggplot2::ggsave(
                file.path(exp$path_post, "Kflood.png"),
                final_plot_Kflood,
                width = 20,
                height = 20,
                units = "cm"
            )

            save(
                final_plot_Kflood,
                file = file.path(
                    exp$path_post_data,
                    "Plot_friction_estimation_plot_Kflood_plot.RData"
                )
            )
        } else {
            Kflood <- NULL
            warning("Any parameter is estimated in the floodplain. Kflood returned is NULL")
        }

        # Residuals
        residuals <- compute_residuals(
            path_temp_plots = exp$path_temp_plots
        )

        save(
            residuals,
            file = file.path(
                exp$path_post_data,
                "Data_Z_sim_vs_obs.RData"
            )
        )
    }
}
