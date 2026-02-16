plot_CalData <- function(
    CalData,
    scales_free = "free_y",
    wrap = TRUE) {
    ##########################################
    # Plot WSE (Output number 1)
    ##########################################
    if (any(CalData$WSE != -9999)) {
        CalData$WSE <- convert_9999_to_NA(CalData$WSE)
        plot_WSE <-
            ggplot(data = CalData, aes(
                x = x,
                y = WSE,
                col = factor(reach)
            ))

        if (any(CalData$Yu_WSE != 0)) {
            plot_WSE <- plot_WSE +
                geom_errorbar(aes(
                    ymin = WSE - qnorm(0.975) * Yu_WSE,
                    ymax = WSE + qnorm(0.975) * Yu_WSE
                ))
        }

        plot_WSE <- plot_WSE +
            geom_point() +
            labs(
                x = "Streamwise position (m)",
                y = "Water surface elevation (m)",
                col = "Reaches in the \nhydraulic model"
            ) +
            theme_bw() +
            theme(
                plot.title = element_text(hjust = 0.5),
                legend.title = element_text(hjust = 0.5)
            )

        if (wrap) {
            plot_WSE <- plot_WSE +
                facet_wrap(~event, scales = scales_free, ncol = 1)
        }

        plot_WSE <- plot_WSE +
            theme(
                strip.text = element_text(size = 12, face = "bold"),
                axis.title.x = element_text(size = 13), # x-axis title
                axis.title.y = element_text(size = 13), # y-axis title
                axis.text.x = element_text(size = 12, hjust = 1), # x-axis ticks
                axis.text.y = element_text(size = 12), # y-axis ticks
                legend.title = element_text(size = 14, face = "bold"),
                legend.text = element_text(size = 12)
            )
    } else {
        plot_WSE <- NULL
    }
    ##########################################
    # Plot Q (Output number 2)
    ##########################################
    if (any(CalData$Q != -9999)) {
        CalData$Q <- convert_9999_to_NA(CalData$Q)
        plot_Q <-
            ggplot(data = CalData, aes(
                x = x,
                y = Q,
                col = factor(reach)
            ))

        if (any(CalData$Yu_Q != 0)) {
            plot_Q <- plot_Q +
                geom_errorbar(aes(
                    ymin = Q - qnorm(0.975) * Yu_Q * Q / 100,
                    ymax = Q + qnorm(0.975) * Yu_Q * Q / 100
                ))
        }

        plot_Q <- plot_Q +
            geom_point() +
            labs(
                x = "Streamwise position (m)",
                y = expression("Discharge (m"^
                    {
                        3
                    } * "/s)"),
                col = "Reaches in the \nhydraulic model"
            ) +
            theme_bw() +
            theme(
                plot.title = element_text(hjust = 0.5),
                legend.title = element_text(hjust = 0.5)
            )
        if (wrap) {
            plot_Q <- plot_Q +
                facet_wrap(~event, scales = scales_free, ncol = 1)
        }

        plot_Q <- plot_Q +
            theme(
                strip.text = element_text(size = 12, face = "bold"),
                axis.title.x = element_text(size = 13), # x-axis title
                axis.title.y = element_text(size = 13), # y-axis title
                axis.text.x = element_text(size = 12, hjust = 1), # x-axis ticks
                axis.text.y = element_text(size = 12), # y-axis ticks
                legend.title = element_text(size = 14, face = "bold"),
                legend.text = element_text(size = 12)
            )
    } else {
        plot_Q <- NULL
    }



    ##########################################
    # Plot V (Output number 3)
    ##########################################
    if (any(CalData$V != -9999)) {
        CalData$V <- convert_9999_to_NA(CalData$V)
        plot_V <-
            ggplot(data = CalData, aes(
                x = x,
                y = V,
                col = factor(reach)
            ))

        if (any(CalData$Yu_V != 0)) {
            plot_V <- plot_V +
                geom_errorbar(aes(
                    ymin = Q - qnorm(0.975) * Yu_V * Q / 100,
                    ymax = Q + qnorm(0.975) * Yu_V * Q / 100
                ))
        }

        plot_V <- plot_V +
            geom_point() +
            labs(
                x = "Streamwise position (m)",
                y = "Velocity (m/s)",
                col = "Reaches in the \nhydraulic model"
            ) +
            theme_bw() +
            theme(
                plot.title = element_text(hjust = 0.5),
                legend.title = element_text(hjust = 0.5)
            )
        if (wrap) {
            plot_V <- plot_V +
                facet_wrap(~event, scales = scales_free, ncol = 1)
        }

        plot_V <- plot_V +
            theme(
                strip.text = element_text(size = 12, face = "bold"),
                axis.title.x = element_text(size = 13), # x-axis title
                axis.title.y = element_text(size = 13), # y-axis title
                axis.text.x = element_text(size = 12, hjust = 1), # x-axis ticks
                axis.text.y = element_text(size = 12), # y-axis ticks
                legend.title = element_text(size = 14, face = "bold"),
                legend.text = element_text(size = 12)
            )
    } else {
        plot_V <- NULL
    }



    ##########################################
    # Plot Kmin (Output number 4)
    ##########################################
    if (any(CalData$Kmin != -9999)) {
        CalData$Kmin <- convert_9999_to_NA(CalData$Kmin)
        plot_Kmin <-
            ggplot(data = CalData, aes(
                x = x,
                y = Kmin,
                col = factor(reach)
            ))

        if (any(CalData$Yu_Kmin != 0)) {
            plot_Kmin <- plot_Kmin +
                geom_errorbar(aes(
                    ymin = Kmin - qnorm(0.975) * Yu_Kmin,
                    ymax = Kmin + qnorm(0.975) * Yu_Kmin
                ))
        }

        plot_Kmin <- plot_Kmin +
            geom_point() +
            labs(
                x = "Streamwise position (m)",
                y = expression("Friction coefficient (m"^
                    {
                        1 / 3
                    } * "/s)"),
                col = "Reaches in the \nhydraulic model"
            ) +
            theme_bw() +
            theme(
                plot.title = element_text(hjust = 0.5),
                legend.title = element_text(hjust = 0.5)
            )

        if (wrap) {
            plot_Kmin <- plot_Kmin +
                facet_wrap(~event, scales = scales_free, ncol = 1)
        }

        plot_Kmin <- plot_Kmin +
            theme(
                strip.text = element_text(size = 12, face = "bold"),
                axis.title.x = element_text(size = 13), # x-axis title
                axis.title.y = element_text(size = 13), # y-axis title
                axis.text.x = element_text(size = 12, hjust = 1), # x-axis ticks
                axis.text.y = element_text(size = 12), # y-axis ticks
                legend.title = element_text(size = 14, face = "bold"),
                legend.text = element_text(size = 12)
            )
    } else {
        plot_Kmin <- NULL
    }


    ##########################################
    # Plot Kflood(Output number 5)
    ##########################################
    if (any(CalData$Kflood != -9999)) {
        CalData$Kflood <- convert_9999_to_NA(CalData$Kflood)
        plot_Kflood <-
            ggplot(data = CalData, aes(
                x = x,
                y = Kflood,
                col = factor(reach)
            ))

        if (any(CalData$Yu_Kflood != 0)) {
            plot_Kflood <- plot_Kflood +
                geom_errorbar(aes(
                    ymin = Kflood - qnorm(0.975) * Yu_Kflood,
                    ymax = Kflood + qnorm(0.975) * Yu_Kflood
                ))
        }

        plot_Kflood <- plot_Kflood +
            geom_point() +
            labs(
                x = "Streamwise position (m)",
                y = expression("Friction coefficient (m"^
                    {
                        1 / 3
                    } * "/s)"),
                col = "Reaches in the \nhydraulic model"
            ) +
            theme_bw() +
            theme(
                plot.title = element_text(hjust = 0.5),
                legend.title = element_text(hjust = 0.5)
            )
        if (wrap) {
            plot_Kflood <- plot_Kflood +
                facet_wrap(~event, scales = scales_free, ncol = 1)
        }

        plot_Kflood <- plot_Kflood +
            theme(
                strip.text = element_text(size = 12, face = "bold"),
                axis.title.x = element_text(size = 13), # x-axis title
                axis.title.y = element_text(size = 13), # y-axis title
                axis.text.x = element_text(size = 12, hjust = 1), # x-axis ticks
                axis.text.y = element_text(size = 12), # y-axis ticks
                legend.title = element_text(size = 14, face = "bold"),
                legend.text = element_text(size = 12)
            )
    } else {
        plot_Kflood <- NULL
    }


    return(list(
        plot_WSE = plot_WSE,
        plot_Q = plot_Q,
        plot_V = plot_V,
        plot_Kmin = plot_Kmin,
        plot_Kflood = plot_Kflood
    ))
}


# Plot spatialisation of k in main channel or floodplain
K_plot <- function(
    matrix_spatialisation,
    mcmc,
    n_param_Kmin,
    n_param_Kflood,
    SU_reaches,
    MAP_param_vector,
    main_channel = TRUE) {
    if (main_channel) {
        indx <- 1:n_param_Kmin
        label_title <- "Friction coefficient estimation in the main channel \nwith parametric uncertainty"
    } else {
        indx <- (n_param_Kmin + 1):(n_param_Kmin + n_param_Kflood)
        label_title <- "Friction coefficient estimation in the floodplain \nwith parametric uncertainty"
    }

    k_estimated_MCMC <- as.data.frame(as.matrix(matrix_spatialisation) %*% as.matrix(t(mcmc[, indx])))

    k_estimated_MCMC$KP <- SU_reaches$KP
    k_estimated_MCMC$reaches_nb <- SU_reaches$reach
    k_estimated_MCMC$reaches_SU <- SU_reaches$id_river

    # Convert to long format
    df_MCMC_sampling <- tidyr::pivot_longer(
        k_estimated_MCMC,
        cols = -c(reaches_nb, reaches_SU, KP),
        values_to = "Value"
    ) %>%
        select(reaches_nb, reaches_SU, KP, Value) %>%
        mutate(ID = "MCMC Sampling")

    k_estimated_MAP <- as.matrix(matrix_spatialisation) %*% MAP_param_vector[indx]

    df_MAP <- data.frame(
        reaches_nb = SU_reaches$reach,
        reaches_SU = SU_reaches$id_river,
        KP = SU_reaches$KP,
        Value = k_estimated_MAP,
        ID = "MAP"
    )

    # Get 95% uncertainty for envelope curve : create ribbon data from MCMC
    df_envelope <- df_MCMC_sampling %>%
        filter(ID == "MCMC Sampling") %>%
        group_by(KP, reaches_SU, reaches_nb) %>%
        summarise(
            ymin = quantile(Value, probs = 0.025, na.rm = TRUE),
            ymax = quantile(Value, probs = 0.975, na.rm = TRUE),
            ID = "Parametric\nuncertainty", # so we can map to fill
            .groups = "drop"
        )

    K_plot <-
        ggplot() +
        geom_ribbon(
            data = df_envelope,
            aes(x = KP, ymin = ymin, ymax = ymax, fill = ID)
        ) +
        geom_line(
            data = df_MAP,
            aes(x = KP, y = Value, color = ID)
        ) +
        scale_fill_manual(values = c("Parametric\nuncertainty" = "pink")) +
        scale_color_manual(values = c("MAP" = "black")) +
        labs(
            title = label_title,
            x = "Streamwise position (meters)",
            y = expression("Friction coefficient (m"^
                {
                    1 / 3
                } * "/s)"),
            fill = "95% credibility\ninterval",
            color = NULL
        ) +
        theme_bw() +
        theme(
            plot.title = element_text(hjust = 0.5),
            legend.title = element_text(hjust = 0.5)
        ) +
        facet_wrap(~reaches_SU, scales = "free", ncol = 1)


    return(list(df_MAP, K_plot, df_envelope))
}

# Plot the simulation of ZQdX
ZQdX_sim <- function(
    sim_event_mm_m3_s,
    Q_input = NULL,
    Qplot = TRUE,
    title_label,
    ylabel) {
    if (Qplot == TRUE && is.null(Q_input)) {
        stop("Q_input must be provided if Qplot is TRUE")
    }

    # Replace -9999 and -1e9 with NA
    sim_event_mm_m3_s[sim_event_mm_m3_s == -9999] <- NA
    sim_event_mm_m3_s[sim_event_mm_m3_s == -1e9] <- NA

    ZQplot <- ggplot(sim_event_mm_m3_s, aes(x = X3_obs)) +
        labs(title = title_label, x = "Streamwise position (meters)", y = ylabel) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))
    if (Qplot == TRUE) {
        ZQplot <- ZQplot +
            geom_point(aes(y = Y2_sim, col = "sim")) +
            scale_color_manual(values = c("sim" = "blue", "boundary\nconditions" = "red")) +
            facet_wrap(~X1_obs, ncol = 2)

        # Create a data frame for the horizontal lines
        hline_data <- data.frame(
            X1_obs = 1:length(Q_input),
            yintercept = Q_input,
            linetype = 1:length(Q_input), # or use named linetypes like c("solid", "dashed", "dotted")
            col = "boundary\nconditions" # same color for all lines
        )

        ZQplot <- ZQplot +
            geom_hline(
                data = hline_data,
                aes(yintercept = yintercept, linetype = factor(linetype), color = col),
                show.legend = TRUE
            )

        ZQplot <- ZQplot +
            scale_linetype_manual(values = c("1" = "dashed", "2" = "solid", "3" = "dotted", "4" = "twodash")) +
            labs(colour = NULL, linetype = "boundary\ncondition\nevents")
    } else {
        ZQplot <- ZQplot +
            geom_point(aes(y = Y1_sim, col = "sim")) +
            geom_point(aes(y = Y1_obs, col = "obs"))

        if (any(sim_event_mm_m3_s$Yu_z != 0)) {
            ZQplot <- ZQplot +
                geom_errorbar(aes(
                    ymin = Y1_obs - qnorm(0.975) * Yu_z,
                    ymax = Y1_obs + qnorm(0.975) * Yu_z, col = "obs"
                ))
        }

        ZQplot <- ZQplot +
            scale_color_manual(values = c("sim" = "blue", "obs" = "red")) +
            labs(colour = NULL) +
            facet_wrap(~X1_obs, ncol = 2)
    }
    return(ZQplot)
}

plot_DIC <- function(
    dir_polynomial,
    DIC_criterion = "DIC3") {
    all_files <- list.files(
        file.path(
            dir_polynomial
        ),
        recursive = TRUE
    )

    DIC_path_Files <- all_files[grep(all_files, pattern = "Results_DIC.txt")]
    if (length(DIC_path_Files) == 0) stop("Any Results_DIC.txt file was found")
    DIC_results <- c()


    for (i in 1:length(DIC_path_Files)) {
        DIC_by_degree <- read.table(file.path(dir_polynomial, DIC_path_Files[i]), col.names = c("Criteria", "Value"))
        # Match the criterion chosen
        DIC_by_degree <- DIC_by_degree[which(DIC_by_degree[, 1] == DIC_criterion), ]
        # Assign the polynomial degree
        extraction <- strsplit(DIC_path_Files[i], "/")[[1]][1]
        # # Extraire le chiffre après le underscore
        # chiffre <- sub(".*_(\\d+)", "\\1", extraction)

        # # Convertir en numérique (optionnel)
        DIC_by_degree$case <- extraction

        DIC_results <- rbind(DIC_results, DIC_by_degree)
    }
    min_local <- DIC_results[which.min(DIC_results$Value), ]

    # Vérifier la présence des motifs
    if (!all(str_detect(DIC_results$case, "Kmin_n_"))) {
        stop("Erreur : certaines lignes n'ont pas 'Kmin_n_' !")
    }

    if (!all(str_detect(DIC_results$case, "Kflood_n_"))) {
        stop("Erreur : certaines lignes n'ont pas 'Kflood_n_' !")
    }

    DIC_plot <- ggplot(DIC_results, aes(x = factor(case), y = Value, col = factor(Criteria))) +
        geom_point(size = 3) +
        geom_point(data = min_local, aes(x = , factor(case), y = Value), col = "blue", size = 3) +
        annotate("segment",
            x = factor(min_local$Degree), y = min_local$Value * 1.005, xend = factor(min_local$Degree), yend = min_local$Value * 1.002,
            linewidth = 2, linejoin = "mitre",
            arrow = arrow(type = "closed", length = unit(0.01, "npc"))
        ) +
        theme_bw() +
        labs(
            x = "Calibration case",
            y = "Value",
            col = "Criterion",
            title = "DIC criterion :",
            subtitle = "Comparison of several polynomial degrees"
        ) +
        theme(
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5)
        )

    return(DIC_plot)
}

# Plot true values of friction (observations)
segment_layer_reference <- function(
    K_literature,
    mean_col,
    min_col = NULL,
    max_col = NULL) {
    # ---- CHECK REQUIRED BASE COLUMNS ----

    required_cols <- c("x_start", "x_end", "reaches_SU")

    missing_required <- setdiff(required_cols, names(K_literature))

    if (length(missing_required) > 0) {
        stop(
            paste(
                "Missing required columns:",
                paste(missing_required, collapse = ", ")
            )
        )
    }


    # ---- CHECK MEAN COLUMN ----

    if (is.null(mean_col)) {
        stop("mean_col must be provided")
    }

    if (!mean_col %in% names(K_literature)) {
        stop(
            paste(
                "mean_col not found in dataframe:",
                mean_col
            )
        )
    }


    # ---- CHECK MIN / MAX ----

    if (!is.null(min_col) && !min_col %in% names(K_literature)) {
        stop(
            paste(
                "min_col not found in dataframe:",
                min_col
            )
        )
    }


    if (!is.null(max_col) && !max_col %in% names(K_literature)) {
        stop(
            paste(
                "max_col not found in dataframe:",
                max_col
            )
        )
    }


    # Warn if only one provided
    if (xor(is.null(min_col), is.null(max_col))) {
        warning(
            "Only one of min_col or max_col provided. Usually both are needed for range."
        )
    }



    # ---- BUILD DATA ----

    K_segment_long <- data.frame()

    for (i in 1:nrow(K_literature)) {
        row <- K_literature[i, ]

        # Mean (main line)
        K_segment_long <- rbind(
            K_segment_long,
            data.frame(
                x_start = row$x_start,
                x_end = row$x_end,
                y_value = row[[mean_col]],
                reaches_SU = row$reaches_SU,
                line_type = "mean"
            )
        )


        # Optional Min
        if (!is.null(min_col) && min_col %in% names(K_literature)) {
            K_segment_long <- rbind(
                K_segment_long,
                data.frame(
                    x_start = row$x_start,
                    x_end = row$x_end,
                    y_value = row[[min_col]],
                    reaches_SU = row$reaches_SU,
                    line_type = "min"
                )
            )
        }


        # Optional Max
        if (!is.null(max_col) && max_col %in% names(K_literature)) {
            K_segment_long <- rbind(
                K_segment_long,
                data.frame(
                    x_start = row$x_start,
                    x_end = row$x_end,
                    y_value = row[[max_col]],
                    reaches_SU = row$reaches_SU,
                    line_type = "max"
                )
            )
        }
    }


    K_segment_layer <- geom_segment(
        data = K_segment_long,
        aes(
            x = x_start,
            y = y_value,
            xend = x_end,
            yend = y_value,
            linetype = line_type
        ),
        color = "gray"
    )

    return(K_segment_layer)
}
