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
                } else if (prefix %in% c("ParamU", "TotalU", "Prior")) {
                    prefix_customized <- ifelse(prefix == "TotalU", "Total", ifelse(prefix == "ParamU", "Parametric", "Prior"))
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
                value = Y_observations[, i],
                "q2.5" = (Y_observations[, i] - z_val * Yu_observations[, i]),
                "q97.5" = (Y_observations[, i] + z_val * Yu_observations[, i]),
                variable = suffix_patterns[i],
                id_pred = "Observations"
            )
        all_CalData <- rbind(all_CalData, CalData)
    }
    all_CalData <- all_CalData %>%
        filter(
            !is.na(value),
            !is.na(q2.5),
            !is.na(q97.5)
        )
    return(all_CalData = all_CalData)
}


postprocess_prediction <- function(
    paths,
    X_input,
    Y_observations,
    Yu_observations,
    date_ref,
    conf_level = 0.95,
    summary_SU_Kmin,
    summary_SU_Kflood,
    grid,
    Input_Typology,
    suffix_patterns = c("_WSE", "_Q", "_V", "_Kmin", "_Kflood"),
    desired_order = c("Total", "Parametric", "Maxpost", "Observations")) {
    # Check if Results_Cooking.txt file exists
    check_calibration_done(path = paths$path_BaM_folder)

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

    all_data <- rbind(all_data_pred, all_CalData_pred)
    all_data$typology <- NA
    for (i in seq_along(Input_Typology)) {
        all_data$typology[all_data$reach %in% Input_Typology[[i]]] <- names(Input_Typology)[i]
    }

    all_data <- all_data %>%
        left_join(summary_SU_Kmin, by = "typology", relationship = "many-to-many") %>%
        filter(x >= KP_min_SU & x <= KP_max_SU) %>%
        distinct(event, reach, x, t, variable, id_pred, .keep_all = TRUE) %>%
        select(-c(KP_max_SU, KP_min_SU)) %>%
        rename(
            id_SU_Kmin = id_SU,
            id_reach_SU_Kmin = id_reach_SU,
        ) %>%
        left_join(summary_SU_Kflood, by = "typology", relationship = "many-to-many") %>%
        filter(x >= KP_min_SU & x <= KP_max_SU) %>%
        distinct(event, reach, x, t, variable, id_pred, .keep_all = TRUE) %>%
        select(-c(KP_max_SU, KP_min_SU)) %>%
        rename(
            id_SU_Kflood = id_SU,
            id_reach_SU_Kflood = id_reach_SU,
        )

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

    all_data_output <- all_data %>%
        # Get information of each output
        filter(variable %in% clean_suffix_patterns) %>%
        rename(
            min = q2.5,
            max = q97.5
        )

    # In the data frame all_data_output, two columns have similar names but two different significations.
    # In fact, var columns indicate the variable that the numerical grid is associated,
    # By contrast, variable column is the predicted variable depending on the id_pred, it exists in all the numerical grid.
    # So if I need to ensure that var == variable to filter the numerical grid corresponding to the expected variable.
    # Ex. The first values of the numerical grid are related to WSE (var), then 5 values of Q.
    # Predicted values has WSE for all 10 values, same for Q. Then, for plotting, I need to filter data to compare the first 5 values of the numerical grid with the 5 predicted values and same for Q.

    # Plot observation and simulation with uncertainties
    plot_unc_by_SU <- plot_obs_sim_unc_lastest(
        data_input = all_data_output,
        wrap = "event_SU",
        ncol = 2
    )

    # Plot observation and simulation with uncertainties (event_reach_HM)
    plot_unc_by_HM <- plot_obs_sim_unc_lastest(
        data_input = all_data_output,
        wrap = "event_reach_HM",
        ncol = 2
    )

    ############################
    # Residuals with uncertainty
    ############################

    all_data_res <-
        all_data %>%
        left_join(
            all_data %>%
                filter(id_pred == "Maxpost") %>%
                select(event, reach, x, t, date_time_format, var, target_datetime, id_campaign, variable,
                    Maxpost = value
                ),
            by = c("event", "reach", "x", "t", "date_time_format", "var", "target_datetime", "id_campaign", "variable")
        ) %>%
        mutate(
            value = value - Maxpost,
            min = q2.5 - Maxpost,
            max = q97.5 - Maxpost
        ) %>% # Remove Maxpost because it is the referent
        filter(id_pred != "Maxpost")

    # Plot residuals with uncertainties by SU
    plot_unc_res_by_SU <- plot_obs_sim_unc_lastest(
        data_input = all_data_res,
        wrap = "event_SU",
        ncol = 2
    )

    # Plot residuals with uncertainties by HM
    plot_unc_res_by_reach <- plot_obs_sim_unc_lastest(
        data_input = all_data_res,
        wrap = "event_reach_HM",
        ncol = 2
    )

    for (i_save in seq_along(plot_unc_by_SU[[1]])) {
        # Save sim vs obs by SU
        ggsave(
            file.path(
                paths$path_plot_folder,
                paste0("sim_vs_obs_by_SU_", plot_unc_by_SU[[2]][i_save], "_with_uncertainties.png")
            ),
            plot_unc_by_SU[[1]][[i_save]],
            width = 17,
            height = 22,
            units = "cm",
            dpi = 300
        )

        # Save sim vs obs by HM
        ggsave(
            file.path(
                paths$path_plot_folder,
                paste0("sim_vs_obs_by_reaches_", plot_unc_by_HM[[2]][i_save], "_with_uncertainties.png")
            ),
            plot_unc_by_HM[[1]][[i_save]],
            width = 17,
            height = 26,
            units = "cm",
            dpi = 300
        )

        # res by SU
        plot_unc_res_by_SU[[1]][[i_save]] <- plot_unc_res_by_SU[[1]][[i_save]] +
            geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
            labs(
                title = paste0(
                    "Residuals with uncertainties by event and SU of :\n",
                    plot_unc_res_by_SU[[2]][i_save]
                ),
                color = "Residuals \n(obs-sim)"
            )

        # Save res by SU
        ggsave(
            filename = file.path(
                paths$path_plot_folder,
                paste0("residual_by_SU_", plot_unc_res_by_SU[[2]][i_save], "_with_uncertainties.png")
            ),
            plot = plot_unc_res_by_SU[[1]][[i_save]],
            dpi = 300,
            width = 20,
            height = 17,
            units = "cm"
        )

        # res by HM
        plot_unc_res_by_SU[[1]][[i_save]] <- plot_unc_res_by_SU[[1]][[i_save]] +
            geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
            labs(
                title = paste0(
                    "Residuals with uncertainties by event and reach of :\n",
                    plot_unc_res_by_reach[[2]][i_save]
                ),
                color = "Residuals \n(obs-sim)"
            )
        # Save res by  HM
        ggsave(
            filename = file.path(
                paths$path_plot_folder,
                paste0("residual_by_reach_", plot_unc_res_by_reach[[2]][i_save], "_with_uncertainties.png")
            ),
            plot = plot_unc_res_by_SU[[1]][[i_save]],
            dpi = 300,
            width = 20,
            height = 17,
            units = "cm"
        )
    }
    # Save RData
    save(plot_unc_by_SU,
        file = file.path(
            paths$path_RData,
            paste0("sim_obs_plot_by_SU.RData")
        )
    )
    save(plot_unc_by_HM,
        file = file.path(
            paths$path_RData,
            paste0("plot_unc_by_reach.RData")
        )
    )
    save(plot_unc_res_by_SU,
        file = file.path(
            paths$path_RData,
            paste0("plot_unc_res_by_SU.RData")
        )
    )
    save(plot_unc_res_by_reach,
        file = file.path(
            paths$path_RData,
            paste0("plot_unc_res_by_reach_with_uncertainties.RData")
        )
    )
}
