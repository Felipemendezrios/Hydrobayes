cat("\014")
rm(list = ls())

# Plots:
library(ggplot2)

workspace <- "/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/Case_studies/HHLab"
setwd(workspace)

# Here, define calibration type
# 'Calibration_time_series'
# 'Calibration_water_surface_profiles'
path_results <- "Calibration_water_surface_profiles"
n_degree_max <- 4
do_prediction <- TRUE

# options  : "Prior","ParamU","Maxpost","TotalU"
all_prediction_case <- c("ParamU", "Maxpost", "TotalU")
# n_prediction <- length(all_prediction_case)


# Adjust size of plot : in pixel!
# Check : do not touch
check_cal_WS_profiles <- path_results == "Calibration_water_surface_profiles"


dimensions_plot_sim_vs_obs <- if (check_cal_WS_profiles) {
    c(
        1800, # Width
        1800 # Height
    )
} else {
    c(
        2000, # Width
        2000 # Height
    )
}
# Do not touche
conf_level <- 0.95 # Hard-coded
z_val <- qnorm(1 - (1 - conf_level) / 2)
# Desired level order
desired_order <- c("TotalU", "ParamU", "Maxpost", "Prior")
n_degree_seq <- seq(0, n_degree_max, 1)

for (n_degree in n_degree_seq) {
    all_uncertainty_data <- c()
    # Polynomial degree i
    path_polynomial <- file.path(path_results, paste0("n_", n_degree))
    path_post_traitement <- file.path(path_polynomial, "post_traitement")
    path_post_traitement_data <- file.path(path_post_traitement, "RData")
    # Read observation
    load(file.path(path_post_traitement, "RData", "output_var_comparison_data.RData"))

    observation_data <- sim_obs_output_variable_long %>%
        select(X1_obs, variable, type, value, value_unc) %>%
        filter(type == "obs") %>%
        rename(id_case = variable) %>%
        mutate(type = recode(type, "obs" = "Observations"))

    number_output_variables <- length(levels(factor(observation_data$id_case)))
    # Read names plot
    load(file.path(path_post_traitement, "RData", "names_plot_ordered.RData"))

    # Load data and model used during calibration
    load(file.path(path_polynomial, "BaM_object_calibration", "BaM_objects.RData"))
    names_Caldata <- colnames(data$data)[2:(number_output_variables + 1)]
    names(name_map) <- names_Caldata

    # Read grid
    grid <- read.table(file.path(path_polynomial, "X1.pred"), header = FALSE)$V1
    if (!dir.exists(path_polynomial)) stop(paste0("Polynomial degree ", n_degree, " is not performed yet. Please do the calibration"))
    # Check
    check_cal_WS_profiles <- path_results == "Calibration_water_surface_profiles"
    legend_name_col <- ifelse(check_cal_WS_profiles, "Calibration : water surface profiles", "Calibration : time series")

    for (prediction_case in all_prediction_case) {
        if (prediction_case == "Maxpost") {
            prediction_files <- list.files(path_polynomial,
                pattern = paste0("^", prediction_case, ".*\\.spag$"),
                full.names = TRUE
            )

            if (length(prediction_files) == 0) stop(paste0("The Prediction case: ", prediction_case, " is not performed yet or the envelop file is not created"))

            for (case in prediction_files) {
                id <- gsub(
                    x = basename(case),
                    pattern = paste0("^", basename(prediction_case), "|\\.spag$"),
                    replacement = ""
                )
                # Get uncertainty
                Uncertainty_env <- read.table(case, header = FALSE)
                Uncertainty_data <- Uncertainty_env %>%
                    rename(Maxpost = V1) %>%
                    mutate(
                        grid = grid,
                        id_case = recode(id, !!!name_map),
                        id_prediction = prediction_case,
                        q2.5 = NA,
                        q97.5 = NA,
                        Maxpost = Maxpost * 1000
                    ) %>%
                    select(id_case, id_prediction, grid, Maxpost, q2.5, q97.5)


                all_uncertainty_data <- rbind(
                    all_uncertainty_data,
                    Uncertainty_data
                )
            }
        } else {
            prediction_files <- list.files(path_polynomial,
                pattern = paste0("^", prediction_case, ".*\\.env$"),
                full.names = TRUE
            )
            if (length(prediction_files) == 0) stop(paste0("The Prediction case: ", prediction_case, " is not performed yet or the envelop file is not created"))

            for (case in prediction_files) {
                id <- gsub(
                    x = basename(case),
                    pattern = paste0("^", basename(prediction_case), "|\\.env$"),
                    replacement = ""
                )
                # Get uncertainty
                Uncertainty_env <- read.table(case, header = TRUE)
                Uncertainty_data <- Uncertainty_env %>%
                    select(q2.5, q97.5) %>%
                    mutate(
                        grid = grid,
                        id_case = recode(id, !!!name_map),
                        id_prediction = prediction_case,
                        Maxpost = NA,
                        q2.5 = q2.5 * 1000,
                        q97.5 = q97.5 * 1000
                    ) %>%
                    select(id_case, id_prediction, grid, Maxpost, q2.5, q97.5)
                all_uncertainty_data <- rbind(
                    all_uncertainty_data,
                    Uncertainty_data
                )
            }
        }
    }

    # Get actual levels present in the data
    present_levels <- intersect(
        desired_order,
        unique(all_uncertainty_data$id_prediction)
    )

    # Apply factor levels only for those that exist
    all_uncertainty_data$id_prediction <- factor(
        all_uncertainty_data$id_prediction,
        levels = present_levels
    )

    # Set plot only with uncertainty data
    sim_obs_plot <- ggplot(data = filter(all_uncertainty_data, !id_prediction %in% c("Maxpost", "Prior")), aes(x = grid)) +
        geom_ribbon(
            aes(
                ymin = q2.5,
                ymax = q97.5,
                fill = id_prediction
            ),
            alpha = 0.65
        )
    # Check if Maxpost is available
    if (any(all_prediction_case == "Maxpost")) {
        sim_obs_plot <- sim_obs_plot +
            geom_line(
                data = subset(all_uncertainty_data, !is.na(Maxpost)),
                aes(y = Maxpost, color = id_prediction),
                linewidth = 0.8
            )
    }
    # Plot simulation with uncertainty and observations
    sim_obs_plot <-
        sim_obs_plot +
        # Add observations
        geom_point(
            data = observation_data,
            aes(x = X1_obs, y = value, col = type)
        ) +
        geom_errorbar(
            data = observation_data,
            aes(
                x = X1_obs,
                ymin = value - z_val * value_unc,
                ymax = value + z_val * value_unc,
                col = type
            ),
            na.rm = TRUE
        ) +
        facet_wrap(~id_case, scales = "free_y") +
        theme_bw() +
        labs(
            x = "Lengthwise position (meters)",
            y = "Stage (mm)",
            title = paste0("Comparison of simulated and observed \nwater surface profiles \n", legend_name_col),
            fill = "Uncertainties",
            color = "Data"
        ) +
        scale_fill_manual(
            values = c(
                "ParamU" = "pink",
                "TotalU" = "red"
            )
        ) +
        scale_color_manual(
            values = c(
                "Observations" = "black",
                "Maxpost" = "blue"
            )
        ) +
        theme(
            strip.text = element_text(size = 12),
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5),
            legend.title = element_text(hjust = 0.5)
        )

    # Save plot and data
    ggsave(
        filename = file.path(path_post_traitement, "output_var_comparison_with_uncertainties.png"),
        plot = sim_obs_plot,
        dpi = 300, width = dimensions_plot_sim_vs_obs[1], height = dimensions_plot_sim_vs_obs[2],
        units = "px"
    )
    save(sim_obs_plot,
        file = file.path(path_post_traitement_data, "output_var_comparison_with_uncertainties_plot.RData")
    )
    if (any(all_prediction_case == "Prior")) {
        ### NOT TESTED YET!
        sim_obs_plot_prior <- ggplot(data = filter(all_uncertainty_data, id_prediction == "Prior"), aes(x = grid)) +
            geom_ribbon(
                aes(
                    ymin = q2.5,
                    ymax = q97.5,
                    fill = id_prediction
                ),
                alpha = 0.65
            )
        # Plot prior simulation and observations
        sim_obs_plot_prior <-
            sim_obs_plot_prior +
            # Add observations
            geom_point(
                data = observation_data,
                aes(x = X1_obs, y = value, col = type)
            ) +
            geom_errorbar(
                data = observation_data,
                aes(
                    x = X1_obs,
                    ymin = value - z_val * value_unc,
                    ymax = value + z_val * value_unc,
                    col = type
                ),
                na.rm = TRUE
            ) +
            facet_wrap(~id_case, scales = "free_y") +
            theme_bw() +
            labs(
                x = "Lengthwise position (meters)",
                y = "Stage (mm)",
                title = paste0("Comparison of prior simulation and observed \nwater surface profiles \n", legend_name_col),
                fill = NULL,
                color = "Data"
            ) +
            scale_fill_manual(
                values = c(
                    "Prior" = "skyblue"
                )
            ) +
            scale_color_manual(
                values = c(
                    "Observations" = "black"
                )
            ) +
            theme(
                strip.text = element_text(size = 12),
                axis.text.x = element_text(angle = 45, hjust = 1),
                plot.title = element_text(hjust = 0.5),
                legend.title = element_text(hjust = 0.5)
            )

        # Save plot and data
        ggsave(
            filename = file.path(path_post_traitement, "prior_output_var_comparison_with_uncertainties.png"),
            plot = sim_obs_plot_prior,
            dpi = 300, width = dimensions_plot_sim_vs_obs[1], height = dimensions_plot_sim_vs_obs[2],
            units = "px"
        )
        save(sim_obs_plot_prior,
            file = file.path(path_post_traitement_data, "prior_output_var_comparison_with_uncertainties_plot.RData")
        )
    }

    # Residual plot
    observation_data_quantiles <- observation_data %>%
        mutate(
            obs_u95_inf = value - z_val * value_unc,
            obs_u95_sup = value + z_val * value_unc
        )

    # Join datasets by id_case and position
    res_data <- observation_data_quantiles %>%
        left_join(
            all_uncertainty_data %>%
                select(id_case, id_prediction, grid, Maxpost) %>%
                filter(!is.na(Maxpost)), # keep only rows with Maxpost
            by = c("id_case" = "id_case", "X1_obs" = "grid")
        ) %>%
        filter(!is.na(Maxpost)) %>% # Ensure Maxpost exists
        mutate(
            residual = value - Maxpost,
            residual_low = obs_u95_inf - Maxpost,
            residual_high = obs_u95_sup - Maxpost
        ) %>%
        select(id_case, X1_obs, id_prediction, residual, residual_low, residual_high)

    res_uncertainty <-
        observation_data_quantiles %>%
        crossing(id_prediction = c("ParamU", "TotalU")) %>% # créer une ligne pour chaque id_prediction
        left_join(
            all_uncertainty_data %>%
                select(id_case, id_prediction, grid, q2.5, q97.5),
            by = c("id_case", "id_prediction", "X1_obs" = "grid")
        ) %>%
        mutate(
            residual_low  = value - q2.5,
            residual_high = value - q97.5
            # type_residual = "Uncertainty"
        ) %>%
        select(id_case, X1_obs, id_prediction, residual_low, residual_high)

    # Puis on rassemble tout (avec les mêmes colonnes clés)
    residuals_all <- bind_rows(
        res_data,
        res_uncertainty
    )

    # Get actual levels present in the data
    present_levels <- intersect(
        desired_order,
        unique(residuals_all$id_prediction)
    )

    # Apply factor levels only for those that exist
    residuals_all$id_prediction <- factor(
        residuals_all$id_prediction,
        levels = present_levels
    )

    # Set residual plot
    res_sim_obs_plot <- ggplot(
        data = filter(residuals_all, !id_prediction %in% c("Maxpost", "Prior")),
        aes(
            x = X1_obs,
            y = residual,
            ymin = residual_low,
            ymax = residual_high,
            fill = id_prediction
        )
    ) +
        # Residual uncertainty bands as ribbon
        geom_ribbon(
            alpha = 0.65
        )

    # Check if Maxpost is available
    if (any(all_prediction_case == "Maxpost")) {
        res_sim_obs_plot <-
            res_sim_obs_plot +
            geom_point(
                data = residuals_all %>%
                    filter(id_prediction == "Maxpost"),
                aes(y = residual, color = "Maxpost")
            ) +
            geom_errorbar(
                data = residuals_all %>%
                    filter(id_prediction == "Maxpost"),
                aes(
                    ymin = residual_low,
                    ymax = residual_high,
                    color = "Maxpost"
                )
            )
    }

    res_sim_obs_plot <- res_sim_obs_plot +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
        facet_wrap(~id_case) +
        theme_bw() +
        scale_fill_manual(
            values = c(
                "ParamU" = "pink",
                "TotalU" = "red"
            )
        ) +
        scale_color_manual(
            values = c(
                "Maxpost" = "black"
            )
        ) +
        labs(
            title = paste0("Residuals water surface profiles with uncertainties \n", legend_name_col),
            x = "Lengthwise position (meters)",
            y = "Residuals (mm)",
            fill = "Uncertainty",
            col = "Residuals \n(obs-sim)"
        ) +
        theme(
            strip.text = element_text(size = 12),
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5),
            legend.title = element_text(hjust = 0.5)
        )


    # Save plot and data
    ggsave(
        filename = file.path(path_post_traitement, "residual_output_var_comparison_with_uncertainties.png"),
        plot = res_sim_obs_plot,
        dpi = 300, width = dimensions_plot_sim_vs_obs[1], height = dimensions_plot_sim_vs_obs[2],
        units = "px"
    )
    save(res_sim_obs_plot,
        file = file.path(path_post_traitement_data, "residual_output_var_comparison_with_uncertainties_plot.RData")
    )
}
