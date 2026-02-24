get_all_data_pred <- function(
    X_grid,
    path,
    patterns = c("_WSE", "_Q", "_V", "_Kmin", "_Kflood")) {
    # List files for each pattern
    file_lists <- lapply(patterns, function(p) {
        list.files(
            path = path,
            pattern = paste0("(", p, "\\.env$|Maxpost", p, ".spag$)"),
            full.names = TRUE
        )
    })
    # Apply function to each element of file_lists
    invisible(lapply(seq_along(file_lists), function(i) {
        check_suffix_pred(file_lists[[i]], i)
    }))

    n_elements <- sum(lengths(file_lists))
    # Initialize an empty data frame
    all_data_pred <- data.frame(
        X_grid[rep(seq_len(nrow(X_grid)), times = n_elements), ],
        value = NA_real_,
        q2.5 = NA_real_,
        q97.5 = NA_real_,
        variable = NA_character_,
        id_pred = NA_character_,
        stringsAsFactors = FALSE
    )
    # Index de remplissage
    row_index <- 1
    # Loop over each pattern and its files
    for (i in seq_along(file_lists)) {
        if (length(file_lists[[i]]) > 0) {
            # Read all files for this pattern
            for (file in file_lists[[i]]) {
                # Extract the prefix using the known pattern
                filename <- basename(file)
                prefix <- sub(paste0(patterns[i], "\\..*"), "", filename)

                if (prefix == "Maxpost") {
                    # Read the file
                    df_temp <- read.table(file, header = FALSE, sep = "", stringsAsFactors = FALSE)

                    df <- data.frame(
                        value = df_temp[, 1],
                        q2.5 = NA_real_,
                        q97.5 = NA_real_,
                        id_pred = prefix
                    )
                } else if (prefix %in% c("ParamU", "TotalU")) {
                    prefix_customized <- ifelse(prefix == "TotalU", "Total", "Parametric")
                    # Read the file
                    df_temp <- read.table(file, header = TRUE, sep = "", stringsAsFactors = FALSE)
                    df <- df_temp %>%
                        dplyr::select(Median, q2.5, q97.5) %>%
                        dplyr::rename(value = Median) %>%
                        dplyr::mutate(id_pred = prefix_customized)
                } else {
                    stop("It must never arrive here")
                }
                n_rows <- nrow(df)
                idx <- row_index:(row_index + n_rows - 1)

                all_data_pred$value[idx] <- df$value
                all_data_pred$q2.5[idx] <- df$q2.5
                all_data_pred$q97.5[idx] <- df$q97.5
                # Add a column for the pattern (without the leading underscore)
                all_data_pred$variable[idx] <- sub("^_", "", patterns[i])
                all_data_pred$id_pred[idx] <- df$id_pred

                row_index <- row_index + n_rows
            }
        }
    }
    return(all_data_pred = all_data_pred)
}

get_all_CalData_pred <- function(
    X_input,
    Y_observations,
    Yu_observations,
    suffix_patterns,
    conf_level = 0.95) {
    Y_observations <- convert_9999_to_NA(Y_observations)
    Yu_observations <- convert_9999_to_NA(Yu_observations)
    # colnames(X_input) <- colnames(grid)

    z_val <- qnorm(1 - (1 - conf_level) / 2)

    all_CalData <- c()
    for (i in seq_along(Y_observations)) {
        CalData <- X_input %>%
            mutate(
                "q2.5" = (Y_observations[, i] - z_val * Yu_observations[, i]),
                "q97.5" = (Y_observations[, i] + z_val * Yu_observations[, i]),
                value = Y_observations[, i],
                variable = suffix_patterns[i],
                id_pred = "Observations"
            )
        all_CalData <- rbind(all_CalData, CalData)
    }
    return(all_CalData = all_CalData)
}


postprocess_prediction <- function(
    paths,
    type = "dX", # dX ou dT
    X_input,
    Y_observations,
    Yu_observations,
    conf_level = 0.95,
    grid,
    Input_Typology,
    suffix_patterns = c("_WSE", "_Q", "_V", "_Kmin", "_Kflood"),
    desired_order = c("Total", "Parametric", "Maxpost", "Observations")) {
    # Check if Results_Cooking.txt file exists
    check_calibration_done(path = paths$path_BaM_folder)
    if (!type %in% c("dX", "dT")) stop(paste0("Type must be either dX or dT. You tapped : ", type))
    message("Processing: ", basename(dirname(paths$path_BaM_folder)))

    all_data_pred <- get_all_data_pred(
        X_grid = grid,
        path = paths$path_BaM_folder,
        patterns = suffix_patterns
    )

    clean_suffix_patterns <- sub("^_", "", suffix_patterns)
    all_CalData_pred <- get_all_CalData_pred(
        X_input = X_input,
        Y_observations = Y_observations,
        Yu_observations = Yu_observations,
        suffix_patterns = clean_suffix_patterns,
        conf_level = conf_level
    )

    colnames(all_CalData_pred)[1:ncol(grid)] <- colnames(grid)

    all_data <- rbind(all_data_pred, all_CalData_pred)
    all_data$SU <- NA
    for (i in seq_along(Input_Typology)) {
        all_data$SU[all_data$reach %in% Input_Typology[[i]]] <- names(Input_Typology)[i]
    }

    # Get actual levels present in the data
    present_levels <- intersect(
        desired_order,
        unique(all_data$id_pred)
    )
    # Apply factor levels only for those that exist
    all_data$id_pred <- factor(
        all_data$id_pred,
        levels = present_levels
    )

    # Plots
    if (type == "dX") {
        all_data$xaxis <- all_data$x
        xlabel <- "Streamwise position (meters)"
    } else {
        all_data$xaxis <- all_data$t / 3600
        xlabel <- "Time (hours)"
    }


    for (i in seq_along(Y_observations)) {
        # Flag to indicate if calibration data is presented of each output variable to plot total or parametric uncertainty. Particular case for Kmin and Kflood, they are not really structural error model, but only parametric.
        any_obs_Y <- any(!is.na(convert_9999_to_NA(Y_observations[i])) & !colnames(Y_observations)[i] %in% c("Kmin", "Kflood"))

        all_data_output <- all_data %>%
            # Get information of each output
            filter(variable == clean_suffix_patterns[i]) %>%
            rename(
                sim_value = value,
                min = q2.5,
                max = q97.5
            )

        # Plot observation and simulation with uncertainties
        plot_unc_by_SU <- plot_obs_sim_unc(
            data_output_var = all_data_output,
            any_obs_Y = any_obs_Y,
            wrap = "event_SU"
        )

        plot_unc_by_SU <- plot_unc_by_SU +
            labs(
                title = paste0("Comparison of simulations and observations by event and SU of:\n", clean_suffix_patterns[i]),
                x = xlabel,
                y = clean_suffix_patterns[i]
            )

        # Save
        ggsave(
            file.path(
                paths$path_plot_folder,
                paste0("sim_vs_obs_by_SU_", clean_suffix_patterns[i], "_with_uncertainties.png")
            ),
            plot_unc_by_SU,
            width = 17,
            height = 22,
            units = "cm",
            dpi = 300
        )

        save(plot_unc,
            file = file.path(
                paths$path_RData,
                paste0("sim_obs_plot_by_SU_", clean_suffix_patterns[i], ".RData")
            )
        )

        plot_unc_by_reach <- plot_obs_sim_unc(
            data_output_var = all_data_output,
            any_obs_Y = any_obs_Y,
            wrap = "event_reach_HM"
        )
        plot_unc_by_reach <- plot_unc_by_reach +
            labs(
                title = paste0("Comparison of simulations and observations by event and reach of :\n", clean_suffix_patterns[i]),
                x = xlabel,
                y = clean_suffix_patterns[i]
            )
        ggsave(
            file.path(
                paths$path_plot_folder,
                paste0("sim_vs_obs_by_reaches_", clean_suffix_patterns[i], "_with_uncertainties.png")
            ),
            plot_unc_by_reach,
            width = 17,
            height = 26,
            units = "cm",
            dpi = 300
        )

        save(plot_unc_by_reach,
            file = file.path(
                paths$path_RData,
                paste0("plot_unc_by_reach_", clean_suffix_patterns[i], ".RData")
            )
        )
        ############################
        # Residuals with uncertainty
        ############################

        all_data_res <- all_data %>%
            left_join(
                all_data %>%
                    filter(id_pred == "Maxpost") %>%
                    select(event, reach, x, t, variable,
                        Maxpost = value
                    ),
                by = c("event", "reach", "x", "t", "variable")
            ) %>%
            mutate(
                sim_value = value - Maxpost,
                min = q2.5 - Maxpost,
                max = q97.5 - Maxpost
            )

        var_output_data_res <- all_data_res %>%
            # Get information of each output
            filter(variable == clean_suffix_patterns[i]) %>%
            # Remove Maxpost because it is the referent
            filter(id_pred != "Maxpost")

        # Plot residuals with uncertainties
        plot_unc_res_by_SU <- plot_obs_sim_unc(
            data_output_var = var_output_data_res,
            any_obs_Y = any_obs_Y,
            wrap = "event_SU"
        )

        plot_unc_res_by_SU <- plot_unc_res_by_SU +
            geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
            labs(
                title = paste0(
                    "Residuals with uncertainties by event and SU of :\n",
                    clean_suffix_patterns[i]
                ),
                x = xlabel,
                y = clean_suffix_patterns[i],
                color = "Residuals \n(obs-sim)"
            )

        # Save plot and data
        ggsave(
            filename = file.path(
                paths$path_plot_folder,
                paste0("residual_by_SU_", clean_suffix_patterns[i], "_with_uncertainties.png")
            ),
            plot = plot_unc_res_by_SU,
            dpi = 300,
            width = 20,
            height = 17,
            units = "cm"
        )
        save(plot_unc_res_by_SU,
            file = file.path(
                paths$path_RData,
                paste0("plot_unc_res_by_SU_", clean_suffix_patterns[i], "_with_uncertainties.RData")
            )
        )

        plot_unc_res_by_reach <- plot_obs_sim_unc(
            data_output_var = var_output_data_res,
            any_obs_Y = any_obs_Y,
            wrap = "event_reach_HM"
        )
        plot_unc_res_by_reach <- plot_unc_res_by_reach +
            geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
            labs(
                title = paste0(
                    "Residuals with uncertainties by event and reach of :\n",
                    clean_suffix_patterns[i]
                ),
                x = xlabel,
                y = clean_suffix_patterns[i],
                color = "Residuals \n(obs-sim)"
            )
        # Save plot and data
        ggsave(
            filename = file.path(
                paths$path_plot_folder,
                paste0("residual_by_reach_", clean_suffix_patterns[i], "_with_uncertainties.png")
            ),
            plot = plot_unc_res_by_reach,
            dpi = 300,
            width = 20,
            height = 17,
            units = "cm"
        )
        save(plot_unc_res_by_reach,
            file = file.path(
                paths$path_RData,
                paste0("plot_unc_res_by_reach_", clean_suffix_patterns[i], "_with_uncertainties.RData")
            )
        )
    }
}
