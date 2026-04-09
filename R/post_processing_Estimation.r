extract_MAP_parameters <- function(
    path_BaM_folder,
    final_calibration,
    Kmin_prior,
    Kflood_prior) {
    Fix_Kmin <- sum(
        get_prior_distribution(Kmin_prior) == "FIX"
    )

    Fix_Kflood <- sum(
        get_prior_distribution(Kflood_prior) == "FIX"
    )

    if (final_calibration) {
        summary <- read.table(
            file.path(
                path_BaM_folder,
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
                path_BaM_folder,
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

compute_K <- function(
    Z_MatrixK,
    mcmc,
    Kmin_prior,
    Kflood_prior,
    K_SU,
    MAP,
    do_main_channel) {
    ## missing a level to search [[1]] ? fixed? only SU_1
    SU <- do.call(
        rbind,
        lapply(names(K_SU), function(typology) {
            do.call(
                rbind,
                lapply(names(K_SU[[typology]]), function(SU) {
                    data.frame(
                        KP = K_SU[[typology]][[SU]]$KP,
                        reach = K_SU[[typology]][[SU]]$reach,
                        typology = typology,
                        id_reach_SU = K_SU[[typology]][[SU]]$id_reach_SU
                    )
                })
            )
        })
    ) %>% arrange(reach)

    res <- K_plot(
        matrix_spatialisation = Z_MatrixK,
        mcmc = mcmc,
        SU_reaches = SU,
        MAP_param_vector = MAP,
        dist_prior_Kmin = get_prior_distribution(Kmin_prior),
        init_guess_prior_Kmin = get_all_init_prior_theta(Kmin_prior),
        dist_prior_Kflood = get_prior_distribution(Kflood_prior),
        init_guess_Kflood_Kmin = get_all_init_prior_theta(Kflood_prior),
        main_channel = do_main_channel
    )

    return(list(
        df_MAP = res[[1]],
        plot = res[[2]],
        envelope = res[[3]]
    ))
}

compute_residuals <- function(
    path_BaM_folder,
    final_calibration,
    path_model_HM_events = NULL,
    Z_MatrixKmin = NULL,
    Z_MatrixKflood = NULL,
    SU_Kmin = NULL,
    SU_Kflood = NULL,
    MAP = NULL,
    mod_polynomials = NULL,
    dir_workspace = NULL,
    X = NULL,
    Y = NULL,
    command_line_MAGE = "") {
    if (final_calibration) {
        residuals <- read.table(
            file.path(
                path_BaM_folder,
                "Results_Residuals.txt"
            ),
            header = TRUE
        )
    } else {
        # I need to pass information of MAP run model and get results if final_calibration !=TRUE
        if (is.null(path_model_HM_events)) {
            stop("If there are not the final calibration results, path_model_HM_events must be given as input")
        }
        # Residuals file is not ready, so I need to create by myself with MAP estimation from sampled data while MCMC is turning
        # Get data format from mage results
        temporal_dir <- tempdir()

        files <- list.files(temporal_dir, full.names = TRUE)

        exclude_file <- file.path(temporal_dir, "vscode-R")

        # List all files in the temp directory
        files <- list.files(temporal_dir, full.names = TRUE)
        # Exclude the specific file
        files_to_delete <- setdiff(files, exclude_file)

        # Remove remaining files if any
        if (length(files_to_delete) > 0) {
            unlink(files_to_delete, recursive = TRUE)
        }

        file.copy(
            from = path_model_HM_events,
            to = temporal_dir, recursive = TRUE
        )
        temp_path <- file.path(temporal_dir, basename(path_model_HM_events))

        path_RUGFile <- file.path(
            temp_path,
            list.files(temp_path)[grep(list.files(temp_path), pattern = ".RUG")]
        )
        MAP_RUGFile <- lapply(path_RUGFile, function(file) {
            read_fortran_data(
                file_path = file,
                col_widths_RUGFile = c(1, 3, 6, 10, 10, 10, 10),
                skip = 1
            )
        })

        param_MAP_values <- get_param_vector_MAP_values(
            SU_Kmin = SU_Kmin,
            SU_Kflood = SU_Kflood,
            MAP = MAP
        )

        values_Kmin_RUGFile <- Z_MatrixKmin %*% param_MAP_values[1:ncol(Z_MatrixKmin)]
        values_Kflood_RUGFile <- Z_MatrixKflood %*% param_MAP_values[(ncol(Z_MatrixKmin) + 1):(ncol(Z_MatrixKmin) + ncol(Z_MatrixKflood))]

        for (i in seq_along(MAP_RUGFile)) {
            MAP_RUGFile[[i]]$Kflood <- values_Kflood_RUGFile
            MAP_RUGFile[[i]]$Kmin <- values_Kmin_RUGFile

            write_RUGFile(
                RUG_path = path_RUGFile[i],
                RUGFile_data = MAP_RUGFile[[i]],
                RUG_format = "%1s%3d      %10.3f%10.3f%10.2f%10.2f"
            )
        }

        REPFile <- list.files(temp_path)[grep(list.files(temp_path), pattern = ".REP")]

        REPFile <- str_remove(REPFile, pattern = ".REP")
        if (length(unique(REPFile)) != 1) {
            stop("names of the REPfile must be identical in all the mage projects")
        }
        REPFile <- unique(REPFile)

        for (id_temp_path in temp_path) {
            setwd(id_temp_path)
            system2(MAGE_executable,
                args = c(
                    file.path(REPFile),
                    command_line_MAGE
                ),
                wait = TRUE
            )
        }

        runModel(
            workspace = temporal_dir,
            mod = mod_polynomials,
            X = X,
            stout = NULL
        )

        sim <- read.table(file.path(temporal_dir, "Y.txt"), header = TRUE)
        colnames(sim) <- paste0("Y", 1:ncol(sim), "_sim")

        # The standard deviation from a single run stdres = std_calibration_data
        u_res <- data.frame(matrix(NA, nrow = nrow(sim), ncol = ncol(sim)))
        colnames(u_res) <- paste0("Y", 1:ncol(u_res), "_stdres")

        # Get observations
        obs <- convert_9999_to_NA(Y)
        res <- obs - sim
        colnames(res) <- paste0("Y", 1:ncol(res), "_res")

        residuals <- cbind(sim, res, u_res)
        setwd(dir_workspace)
    }

    # Colonnes Y_sim
    Y_sim_cols <- grep("^Y[0-9]+_sim$", colnames(residuals), value = TRUE)

    # Colonnes Y_res
    Y_res_cols <- grep("^Y[0-9]+_res$", colnames(residuals), value = TRUE)

    # Colonnes Y_stdres
    Y_stdres_cols <- grep("^Y[0-9]+_stdres$", colnames(residuals), value = TRUE)

    res <- residuals[, c(
        Y_sim_cols,
        Y_res_cols,
        Y_stdres_cols
    )]

    return(res)
}

# MAIN FUNCTION
postprocess_calibration <- function(
    paths,
    X_input,
    Y_observations,
    Yu_observations,
    type = "dx",
    final_calibration = TRUE,
    Kmin_prior,
    Kflood_prior,
    Kmin_SU,
    Kflood_SU,
    Z_MatrixKmin,
    Z_MatrixKflood,
    Key_Info_Typology_Model_Reach,
    summary_SU_Kmin,
    summary_SU_Kflood,
    mod_polynomials = NULL,
    Kmin_segment_layer = NULL,
    Kflood_segment_layer = NULL,
    dir_workspace = NULL,
    command_line_MAGE = "") {
    message("Processing: ", basename(dirname(paths$path_BaM_folder)))

    # MCMC
    mcmc <- plot_mcmc_diagnostics(
        path_BaM_folder = paths$path_BaM_folder,
        path_plot_folder = paths$path_plot_folder,
        final_calibration = final_calibration
    )

    # MAP
    MAP <- extract_MAP_parameters(
        path_BaM_folder = paths$path_BaM_folder,
        final_calibration = final_calibration,
        Kmin_prior = Kmin_prior,
        Kflood_prior = Kflood_prior
    )

    CalData <- convert_9999_to_NA(cbind(X_input, Y_observations, Yu_observations))

    summary_HM <- do.call(
        rbind,
        lapply(names(Key_Info_Typology_Model_Reach), function(typology) {
            nodes_HM <- Key_Info_Typology_Model_Reach[[typology]]$Model_Reach_nodes
            reach_HM <- unique(Key_Info_Typology_Model_Reach[[typology]]$reach)
            # each segment is from MR_nodes[i] to MR_nodes[i+1]
            data.frame(
                typology = typology,
                id_reach_HM = reach_HM
            )
        })
    )
    CalData_HM <- CalData %>%
        left_join(summary_HM, by = c("reach" = "id_reach_HM"))

    CalData_updated <- CalData_HM %>%
        left_join(summary_SU_Kmin, by = "typology", relationship = "many-to-many") %>%
        filter(x >= KP_min_SU & x <= KP_max_SU) %>%
        distinct(event, reach, x, .keep_all = TRUE) %>%
        select(-c(KP_max_SU, KP_min_SU)) %>%
        rename(
            id_SU_Kmin = id_SU,
            id_reach_SU_Kmin = id_reach_SU,
        )

    CalData_updated <- CalData_updated %>%
        left_join(summary_SU_Kflood, by = "typology", relationship = "many-to-many") %>%
        filter(x >= KP_min_SU & x <= KP_max_SU) %>%
        distinct(event, reach, x, .keep_all = TRUE) %>%
        select(-c(KP_max_SU, KP_min_SU)) %>%
        rename(
            id_SU_Kflood = id_SU,
            id_reach_SU_Kflood = id_reach_SU,
        )

    # Kmin
    Kmin <- compute_K(
        Z_MatrixK = Z_MatrixKmin,
        mcmc = mcmc,
        Kmin_prior = Kmin_prior,
        Kflood_prior = Kflood_prior,
        K_SU = Kmin_SU,
        MAP = MAP,
        do_main_channel = TRUE
    )

    ls_spatial_friction_Kmin <- list(
        df_envelope = Kmin[[3]],
        df_MAP = Kmin[[1]]
    )
    save(ls_spatial_friction_Kmin,
        file = file.path(paths$path_RData, "Data_friction_estimation_ls_spatial_friction_Kmin.RData")
    )

    if (!is.null(Kmin_segment_layer)) {
        final_plot_Kmin <- plot_K_and_ref(
            K_results = Kmin,
            K_segment_layer = Kmin_segment_layer,
            path_plot_folder = paths$path_plot_folder,
            path_RData = paths$path_RData
        )
    } else {
        final_plot_Kmin <- Kmin$plot
    }

    ggplot2::ggsave(
        file.path(paths$path_plot_folder, "Kmin.png"),
        final_plot_Kmin,
        width = 20,
        height = 20,
        units = "cm"
    )

    save(
        final_plot_Kmin,
        file = file.path(
            paths$path_RData,
            "Plot_friction_estimation_plot_Kmin_plot.RData"
        )
    )

    # Kflood

    Kflood <- compute_K(
        Z_MatrixK = Z_MatrixKflood,
        mcmc = mcmc,
        Kmin_prior = Kmin_prior,
        Kflood_prior = Kflood_prior,
        K_SU = Kflood_SU,
        MAP = MAP,
        do_main_channel = FALSE
    )

    ls_spatial_friction_Kflood <- list(
        df_envelope = Kflood[[3]],
        df_MAP = Kflood[[1]]
    )
    save(ls_spatial_friction_Kflood,
        file = file.path(paths$path_RData, "Data_friction_estimation_ls_spatial_friction_Kflood.RData")
    )

    if (!is.null(Kflood_segment_layer)) {
        final_plot_Kflood <- plot_K_and_ref(
            K_results = Kflood,
            K_segment_layer = Kflood_segment_layer,
            path_plot_folder = paths$path_plot_folder,
            path_RData = paths$path_RData
        )
    } else {
        final_plot_Kflood <- Kflood$plot
    }

    ggplot2::ggsave(
        file.path(paths$path_plot_folder, "Kflood.png"),
        final_plot_Kflood,
        width = 20,
        height = 20,
        units = "cm"
    )

    save(
        final_plot_Kflood,
        file = file.path(
            paths$path_RData,
            "Plot_friction_estimation_plot_Kflood_plot.RData"
        )
    )


    # Residuals
    residuals <- compute_residuals(
        path_BaM_folder = paths$path_BaM_folder,
        final_calibration = final_calibration,
        path_model_HM_events = paths$path_model_HM_events,
        Z_MatrixKmin = Z_MatrixKmin,
        Z_MatrixKflood = Z_MatrixKflood,
        SU_Kmin = Kmin_SU,
        SU_Kflood = Kflood_SU,
        MAP = MAP,
        mod_polynomials = mod_polynomials,
        dir_workspace = dir_workspace,
        X = X_input,
        Y = Y_observations,
        command_line_MAGE = command_line_MAGE
    )

    colnames(X_input) <- paste0("X", 1:ncol(Y_observations), "_obs")
    colnames(Y_observations) <- paste0("Y", 1:ncol(Y_observations), "_obs")
    colnames(Yu_observations) <- paste0("Yu", 1:ncol(Yu_observations), "_obs")
    obs_sim_residuals <- cbind(
        X_input,
        Y_observations,
        Yu_observations,
        residuals
    )
    save(
        obs_sim_residuals,
        file = file.path(
            paths$path_RData,
            "Data_Z_sim_vs_obs.RData"
        )
    )

    # Plot for each output variable
    plots <- plot_obs_sim_MAP(
        all_obs_simulations = obs_sim_residuals,
        type = type
    )

    return(
        list(
            data_param = list(
                Kmin = list(
                    df_MAP = Kmin[[1]],
                    envelop = Kmin[[3]]
                ),
                Kflood = list(
                    df_MAP = Kflood[[1]],
                    envelop = Kflood[[3]]
                )
            ),
            residuals = residuals,
            plots_param = list(
                Kmin = list(
                    plot_without_obs = Kmin[[2]],
                    plot_with_obs = final_plot_Kmin
                ),
                Kflood = list(
                    plot_without_obs = Kflood[[2]],
                    plot_with_obs = final_plot_Kflood
                )
            ),
            plots_MAP_output_variables = plots,
            CalData_updated = CalData_updated
        )
    )
}
