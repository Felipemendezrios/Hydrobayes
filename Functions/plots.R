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
                    ymin = Y1_obs - 1.96 * Yu_z,
                    ymax = Y1_obs + 1.96 * Yu_z, col = "obs"
                ))
        }

        ZQplot <- ZQplot +
            scale_color_manual(values = c("sim" = "blue", "obs" = "red")) +
            labs(colour = NULL) +
            facet_wrap(~X1_obs, ncol = 2)
    }
    return(ZQplot)
}

# Plot spatialisation of k in main channel or floodplain
k_plot <- function(
    matrix_spatialisation,
    mcmc,
    n_degree_kmin,
    n_degree_kflood,
    covariate_discretization,
    MAP_param_vector,
    main_channel = TRUE,
    ks_literature) {
    if (main_channel) {
        indx <- 1:(n_degree_kmin + 1)
        label_title <- "Friction coefficient estimation in the main channel \nwith parametric uncertainty"
    } else {
        indx <- (n_degree_kmin + 1 + 1):(n_degree_kmin + 1 + n_degree_kflood + 1)
        label_title <- "Friction coefficient estimation in the floodplain \nwith parametric uncertainty"
    }
    k_estimated_MCMC <- as.data.frame(as.matrix(matrix_spatialisation) %*% as.matrix(t(mcmc[, indx])))

    k_estimated_MCMC$KP <- covariate_discretization

    # Convert to long format
    df_MCMC_sampling <- tidyr::pivot_longer(
        k_estimated_MCMC,
        cols = -KP,
        values_to = "Value"
    ) %>%
        select(KP, Value) %>%
        mutate(ID = "MCMC Sampling")

    k_estimated_MAP <- as.matrix(matrix_spatialisation) %*% MAP_param_vector[indx]

    df_MAP <- data.frame(
        KP = covariate_discretization,
        Value = k_estimated_MAP,
        ID = "MAP"
    )

    # Get 95% uncertainty for envelope curve : create ribbon data from MCMC
    df_envelope <- df_MCMC_sampling %>%
        filter(ID == "MCMC Sampling") %>%
        group_by(KP) %>%
        summarise(
            ymin = quantile(Value, probs = 0.025, na.rm = TRUE),
            ymax = quantile(Value, probs = 0.975, na.rm = TRUE),
            ID = "Parametric\nuncertainty", # so we can map to fill
            .groups = "drop"
        )

    k_plot <- ggplot() +
        geom_ribbon(
            data = df_envelope,
            aes(x = KP, ymin = ymin, ymax = ymax, fill = ID)
        ) +
        geom_line(
            data = df_MAP,
            aes(x = KP, y = Value, color = ID)
        ) +
        geom_hline(aes(yintercept = ks_literature$min, linetype = "Chow (1959)"), color = "gray") +
        geom_hline(aes(yintercept = ks_literature$max, linetype = "Chow (1959)"), color = "gray") +
        geom_hline(aes(yintercept = ks_literature$mean, linetype = "Chow (1959)"), color = "gray") +
        scale_fill_manual(values = c("Parametric\nuncertainty" = "pink")) +
        scale_color_manual(values = c("MAP" = "black")) +
        # Échelle de type de ligne
        scale_linetype_manual(
            name = "Reference\nvalues",
            values = c("Chow (1959)" = "dashed")
        ) +
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
        )


    return(list(df_MAP, k_plot, df_envelope))
}


CalData_plot <- function(
    data,
    scales_free = "free_y",
    y_label,
    title_label,
    col_label = NULL,
    plot_water_depth = TRUE,
    wrap = TRUE) {
    if (plot_water_depth) {
        plot_CalData <-
            ggplot(data, aes(
                x = x,
                color = ID,
                y = h_mean
            ))

        if (any(data$Yu != 0)) {
            plot_CalData <- plot_CalData +
                geom_errorbar(aes(
                    ymin = h_mean - 1.96 * Yu,
                    ymax = h_mean + 1.96 * Yu
                ))
        }
    } else {
        plot_CalData <-
            ggplot(data = data, aes(
                x = x,
                y = z_mean,
                col = ID
            )) +
            geom_line(
                aes(
                    y = z_riverbed,
                    col = "riverbed"
                )
            )

        if (any(data$Yu != 0)) {
            plot_CalData <- plot_CalData +
                geom_errorbar(aes(
                    ymin = z_mean - 1.96 * Yu,
                    ymax = z_mean + 1.96 * Yu
                ))
        }
    }
    plot_CalData <- plot_CalData +
        geom_point() +
        labs(
            x = "Streamwise position (meters)",
            y = y_label,
            title = title_label,
            col = col_label
        ) +
        theme_bw() +
        theme(
            plot.title = element_text(hjust = 0.5),
            legend.title = element_text(hjust = 0.5)
        )

    if (wrap) {
        plot_CalData <- plot_CalData +
            facet_wrap(~ID, scales = scales_free, ncol = 1)
    }
    return(plot_CalData)
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
        # Extraire le chiffre après le underscore
        chiffre <- sub(".*_(\\d+)", "\\1", extraction)

        # Convertir en numérique (optionnel)
        DIC_by_degree$Degree <- as.numeric(chiffre)

        DIC_results <- rbind(DIC_results, DIC_by_degree)
    }
    min_local <- DIC_results[which.min(DIC_results$Value), ]

    DIC_plot <- ggplot(DIC_results, aes(x = factor(Degree), y = Value, col = factor(Criteria))) +
        geom_point(size = 3) +
        geom_point(data = min_local, aes(x = , factor(Degree), y = Value), col = "blue", size = 3) +
        annotate("segment",
            x = factor(min_local$Degree), y = min_local$Value * 1.005, xend = factor(min_local$Degree), yend = min_local$Value * 1.002,
            linewidth = 2, linejoin = "mitre",
            arrow = arrow(type = "closed", length = unit(0.01, "npc"))
        ) +
        theme_bw() +
        labs(
            x = "Polynomial degree",
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
