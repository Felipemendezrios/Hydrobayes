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
    check_calibration_done(path = paths$path_temp_plots)
    if (!type %in% c("dX", "dT")) stop(paste0("Type must be either dX or dT. You tapped : ", type))
    message("Processing: ", basename(dirname(paths$path_temp_plots)))

    all_data_pred <- get_all_data_pred(
        X_grid = grid,
        path = paths$path_temp_plots,
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
        all_data_output <- all_data %>%
            # Get information of each output
            filter(variable == clean_suffix_patterns[i])


        all_data_unc_obs <- all_data_output %>%
            filter(!id_pred %in% c("Maxpost", "Prior", "Observations"))

        # Plot simulation uncertainties ribbons
        sim_obs_plot <-
            ggplot(
                data = all_data_unc_obs,
                aes(
                    x = xaxis,
                    ymin = q2.5,
                    ymax = q97.5,
                ), alpha = 0.65
            ) +
            geom_ribbon(
                aes(fill = id_pred),
                alpha = 0.65
            ) +
            facet_wrap(~ event + SU,
                scales = "free",
                ncol = 1
            ) +
            theme_bw() +
            labs(
                x = xlabel,
                y = clean_suffix_patterns[i],
                title = paste0("Comparison of simulations and observations by event and SU of :\n", clean_suffix_patterns[i]),
                fill = "Uncertainties",
                color = "Data"
            ) +
            theme(
                strip.text = element_text(size = 12),
                axis.text.x = element_text(angle = 45, hjust = 1),
                plot.title = element_text(hjust = 0.5),
                legend.title = element_text(hjust = 0.5),
                legend.position = "bottom"
            )
        # Check if Maxpost is available
        if (any(levels(all_data_output$id_pred) == "Maxpost")) {
            sim_obs_plot <- sim_obs_plot +
                geom_line(
                    data = all_data_output %>% filter(id_pred == "Maxpost"),
                    aes(y = value, color = id_pred)
                )
        }
        # Plot observations
        sim_obs_plot <- sim_obs_plot +
            # Add observations
            geom_point(
                data = all_data_output %>% filter(id_pred == "Observations"),
                aes(y = value, col = id_pred)
            ) +
            geom_errorbar(
                data = all_data_output %>% filter(id_pred == "Observations"),
                aes(
                    col = id_pred
                ),
                na.rm = TRUE
            ) +
            scale_fill_manual(
                values = c(
                    "Parametric" = "pink",
                    "Total" = "red"
                )
            ) +
            scale_color_manual(
                values = c(
                    "Observations" = "black",
                    "Maxpost" = "blue"
                )
            )

        # Save
        ggsave(
            file.path(
                paths$path_post,
                paste0("sim_vs_obs_by_SU_", clean_suffix_patterns[i], "_with_uncertainties.png")
            ),
            sim_obs_plot,
            width = 17,
            height = 22,
            units = "cm",
            dpi = 300
        )

        save(sim_obs_plot,
            file = file.path(
                paths$path_post_data,
                paste0("sim_obs_plot_by_SU_", clean_suffix_patterns[i], ".RData")
            )
        )

        # Plot by reach (mage)
        sim_obs_plot_by_reaches <- sim_obs_plot +
            facet_wrap(
                ~ event + reach,
                scales = "free",
                ncol = 2
            )

        ggsave(
            file.path(
                paths$path_post,
                paste0("sim_vs_obs_by_reaches_mage_", clean_suffix_patterns[i], "_with_uncertainties.png")
            ),
            sim_obs_plot_by_reaches,
            width = 17,
            height = 26,
            units = "cm",
            dpi = 300
        )

        save(sim_obs_plot_by_reaches,
            file = file.path(
                paths$path_post_data,
                paste0("sim_obs_plot_by_reaches_mage_", clean_suffix_patterns[i], ".RData")
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
                residual = value - Maxpost,
                residual_low = q2.5 - Maxpost,
                residual_high = q97.5 - Maxpost
            )

        var_output_data_res <- all_data_res %>%
            # Get information of each output
            filter(variable == clean_suffix_patterns[i])

        all_unc_var_output_res <- var_output_data_res %>%
            filter(!id_pred %in% c("Maxpost", "Prior", "Observations"))
        # Set residual plot
        res_sim_obs_plot <-
            ggplot(
                data = all_unc_var_output_res,
                aes(
                    x = xaxis,
                    y = residual,
                    ymin = residual_low,
                    ymax = residual_high
                )
            ) +
            # Residual uncertainty bands as ribbon
            geom_ribbon(
                aes(fill = id_pred),
                alpha = 0.65
            ) +
            facet_wrap(~ event + SU,
                scales = "free",
                ncol = 1
            ) +
            theme_bw() +
            labs(
                x = xlabel,
                y = clean_suffix_patterns[i],
                title = paste0("Residuals with uncertainties by event and SU of :\n", clean_suffix_patterns[i]),
                fill = "Uncertainties",
                color = "Residuals \n(obs-sim)"
            ) +
            theme(
                strip.text = element_text(size = 12),
                axis.text.x = element_text(angle = 45, hjust = 1),
                plot.title = element_text(hjust = 0.5),
                legend.title = element_text(hjust = 0.5),
                legend.position = "bottom"
            )

        # Plot observations
        res_sim_obs_plot <- res_sim_obs_plot +
            geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
            # Add observations
            geom_point(
                data = var_output_data_res %>% filter(id_pred == "Observations"),
                aes(col = id_pred)
            ) +
            geom_errorbar(
                data = var_output_data_res %>% filter(id_pred == "Observations"),
                aes(
                    col = id_pred
                ),
                na.rm = TRUE
            ) +
            scale_fill_manual(
                values = c(
                    "Parametric" = "pink",
                    "Total" = "red"
                )
            ) +
            scale_color_manual(
                values = c(
                    "Observations" = "black"
                )
            )

        # Save plot and data
        ggsave(
            filename = file.path(paths$path_post, paste0("residual_", clean_suffix_patterns[i], "_with_uncertainties.png")),
            plot = res_sim_obs_plot,
            dpi = 300,
            width = 20,
            height = 17,
            units = "cm"
        )
        save(res_sim_obs_plot,
            file = file.path(paths$path_post_data, paste0("res_sim_obs_plot_", clean_suffix_patterns[i], "_with_uncertainties.RData"))
        )
    }
}
