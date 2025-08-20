# Plot the residual of ZQdX
ZQdX_residuals <- function(
    residuals_event_mm_m3_s,
    Q_input = NULL,
    Qplot = TRUE,
    title_label,
    xlabel = "Position (m)",
    ylabel) {
    if (Qplot == TRUE && is.null(Q_input)) {
        stop("Q_input must be provided if Qplot is TRUE")
    }

    # Replace -9999 and -1e9 with NA
    residuals_event_mm_m3_s[residuals_event_mm_m3_s == -9999] <- NA
    residuals_event_mm_m3_s[residuals_event_mm_m3_s == -1e9] <- NA

    ZQplot <- ggplot(residuals_event_mm_m3_s, aes(x = X3_obs)) +
        labs(title = title_label, x = xlabel, y = ylabel) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))
    if (Qplot == TRUE) {
        ZQplot <- ZQplot +
            geom_point(aes(y = Y2_sim, col = "sim")) +
            geom_hline(aes(yintercept = Q_input, col = "boundary\nconditions", linetype = "1")) +
            scale_color_manual(values = c("sim" = "blue", "boundary\nconditions" = "red")) +
            scale_linetype_manual(values = c("1" = "dashed", "2" = "solid", "3" = "dotted", "4" = "twodash")) +
            labs(colour = NULL, linetype = "boundary\ncondition\nevents") +
            facet_wrap(~X1_obs, ncol = 2)
    } else {
        ZQplot <- ZQplot +
            geom_point(aes(y = Y1_sim, col = "sim")) +
            geom_point(aes(y = Y1_obs, col = "obs"))

        if (any(residuals_event_mm_m3_s$Yu_z != 0)) {
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
    } else {
        indx <- (n_degree_kmin + 1 + 1):(n_degree_kmin + 1 + n_degree_kflood + 1)
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
            title = "Friction coefficient estimation in the floodplain \nwith parametric uncertainty",
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
    return(k_plot)
}
