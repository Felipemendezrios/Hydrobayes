rm(list = ls())
graphics.off()

function_list <- list.files("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/Functions/", full.names = TRUE)
for (i in function_list) {
    source(i)
}


# Libraries
library(RBaM)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(stringr)

# Set directory
dir_workspace <- here::here()

do_prediction <- TRUE
do_plots <- TRUE

prediction_file <- c("Prior", "ParamU", "Maxpost", "TotalU")
n_prediction <- length(prediction_file)
nsim_prior <- 500 # number of simulation for propagation (def : 100 for priori prediction)

n_degree_max <- 3
Experiment_id <- "1_WSE_main_channel_real_uncertainty"
path_experiment <- file.path("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/Case_studies/HHLab/Rectangular_channel_1_event/Calibration_experiments", Experiment_id)


# Define the calculation grid: Extract from net folder and ST File

# Events and number of reaches in total
nb_events <- c(1)
nb_reaches <- 1

# Define xmin and xmax for each reach
# xmin_values <- c(0, 10) # xmin for reach 1 and reach 2 if nb_reaches = 2
xmin_values <- c(0.06) # Minimum value must higher or equal to the first Calibration data or second (?) cross section in the model
xmax_values <- c(15)


if (nb_reaches != length(xmin_values) || nb_reaches != length(xmax_values)) {
    stop("nb_reaches must be equal to the number of element in xmin_values and xmax_values")
}
grid_user <- data.frame(
    event = rep(nb_events, times = nb_reaches),
    reach = rep(1:nb_reaches, each = length(nb_events)),
    xmin = rep(xmin_values, each = length(nb_events)),
    xmax = rep(xmax_values, each = length(nb_events))
)

# Order by event
grid_user <- grid_user[order(grid_user$event), ]

# Number of discretization points
n_points <- 100

# Generate discretization for each row and expand into a long-format data frame
X_temp <- do.call(rbind, lapply(1:nrow(grid_user), function(i) {
    data.frame(
        event = grid_user$event[i],
        reach = grid_user$reach[i],
        x = seq(from = grid_user$xmin[i], to = grid_user$xmax[i], length.out = n_points)
        # x = seq(0.06,16)
    )
}))

X_temp$t <- 3420 # any number to complet the input variables to keep the format

# Do not touch
n_degree_seq <- seq(0, n_degree_max, 1)
dir_exe_BaM <- "/home/famendezrios/Documents/Git/BaM/makefile/"
nY <- 5
###################################################

for (n_degree in n_degree_seq) {
    # Intitialization
    pred_confi_file_name <- c() # vector with setting prediction file name (Config_Pred_"type of prediction")
    doStructural_logical <- c() # logical vector with setting in function of type of prediction
    doParametric_logical <- c() # logical vector with setting in function of type of prediction
    priorNsim_int <- c() # list of number of simulation for
    pred_list <- c() # list with information about all predictions
    pred_var_name <- c() # list with names of each variable
    mod_list <- list() # list to write several models for running in parallel

    # Polynomial degree i
    path_polynomial <- file.path(path_experiment, paste0("n_", n_degree), "Calibration")
    if (!dir.exists(path_polynomial)) stop(paste0("Polynomial degree ", n_degree, " is not performed yet. Please do the calibration"))

    if (!file.exists(file.path(
        path_polynomial,
        "Results_Cooking.txt"
    ))) {
        stop("Calibration is not performed yet")
    }

    # Load data and model used during calibration
    path_BaM_object <- file.path(path_polynomial, "post_traitement", "RData")
    load(file.path(path_BaM_object, "BaM_objects.RData"))

    Caldata <- data$data
    X <- full_join(X_temp, Caldata[, 1:4], by = c("event", "reach", "x", "t")) %>%
        arrange(event, reach, x)

    names_file_prediction <- colnames(Caldata[, !grepl(colnames(Caldata), pattern = "^Yu_")])[c(-1:-4)]
    for (i in 1:nY) {
        # Initialize temporary storage for prediction file names
        pred_var_name_temps <- c()

        for (j in 1:n_prediction) {
            # File naming
            pred_var_name_temps <- cbind(
                pred_var_name_temps,
                paste0(
                    prediction_file[j], "_",
                    names_file_prediction[i], ".spag"
                )
            )

            pred_confi_file_name[j] <- c(paste0(
                "Config_Pred_",
                prediction_file[j],
                ".txt"
            ))
            #  Logical flags for structural and parametric predictions
            doStructural_logical[[j]] <- switch(prediction_file[j],
                "Prior" = c(rep(FALSE, nY)),
                "ParamU" = c(rep(FALSE, nY)),
                "Maxpost" = c(rep(FALSE, nY)),
                "TotalU" = c(rep(TRUE, nY)),
            )
            doParametric_logical[j] <- switch(prediction_file[j],
                "Prior" = c(TRUE),
                "ParamU" = c(TRUE),
                "Maxpost" = c(FALSE),
                "TotalU" = c(TRUE),
            )
            # Simulation count for Prior
            priorNsim_int[j] <- switch(prediction_file[j],
                "Prior" = nsim_prior,
                "ParamU" = c(-1),
                "Maxpost" = c(-1),
                "TotalU" = c(-1),
            )
            # To run in parallel, it is necessary to have several mage models to avoid write in the same file
            mod_list[[j]] <- mod
            # Move model_mage folder to the result folder
            dir_origin <- mod_list[[j]]$xtra$object$mageDir
            dir_copied_all <- file.path(normalizePath(dirname(dir_origin)), prediction_file[j], basename(dir_origin))
            for (source_dir in dir_origin) {
                source_dir <- normalizePath(source_dir)
                dir_name <- basename(source_dir)
                dir_copied <- file.path(dirname(source_dir), prediction_file[j], dir_name)

                if (dir.exists(dir_copied)) {
                    unlink(dir_copied, recursive = TRUE)
                }

                dir.create(dir_copied, recursive = TRUE)

                files_to_copy <- list.files(
                    source_dir,
                    full.names = TRUE,
                    recursive = TRUE,
                    include.dirs = TRUE
                )

                # Copie chaque fichier/sous-dossier vers dir_copied
                for (file in files_to_copy) {
                    # Chemin relatif du fichier/sous-dossier par rapport à source_dir
                    relative_path <- sub(source_dir, "", file)
                    # Chemin de destination final
                    dest_path <- file.path(dir_copied, relative_path)

                    # Crée les sous-dossiers nécessaires dans la destination
                    if (dir.exists(file)) {
                        if (!dir.exists(dest_path)) {
                            dir.create(dest_path, recursive = TRUE)
                        }
                    } else {
                        if (!dir.exists(dirname(dest_path))) {
                            dir.create(dirname(dest_path), recursive = TRUE)
                        }
                        file.copy(file, dest_path)
                    }
                }
            }
            mod_list[[j]]$xtra$object$mageDir <- paste0(dir_copied_all, "/")
            mod_list[[j]]$xtra$fname <- sub("\\.txt$", paste0("_", prediction_file[j], ".txt"), mod$xtra$fname)
        }
        pred_var_name <- rbind(pred_var_name, pred_var_name_temps)
    }
    # create a pred_list only for spagfile and another for fname

    for (i in 1:n_prediction) {
        pred_list[[i]] <- prediction(
            X = X,
            spagFiles = pred_var_name[, i],
            data.dir = workspace_user,
            fname = pred_confi_file_name[i],
            priorNsim = priorNsim_int[i],
            doParametric = doParametric_logical[i],
            doStructural = doStructural_logical[[i]]
        )
        # Run BaM, but in prediction mode
        BaM(
            mod = mod_list[[i]],
            data = data,
            remnant = remant_error_list,
            mcmc = mcmcOptions_user,
            cook = mcmcCooking_user,
            summary = mcmcSummary_user,
            residuals = residualOptions(),
            pred = pred_list[[i]],
            doCalib = FALSE,
            doPred = TRUE,
            na.value = -9999,
            run = FALSE,
            preClean = FALSE,
            workspace = path_polynomial,
            dir.exe = dir_exe_BaM,
            name.exe = "BaM",
            predMaster_fname = paste0("Config_Pred_Master_", prediction_file[i], ".txt")
        )

        # Rename Config_BaM.txt file to the result folder
        cf_file <- file.path(path_polynomial, paste0("Config_BaM_", prediction_file[i], ".txt"))
        if (file.exists(file.path(path_polynomial, "Config_BaM.txt"))) {
            invisible(file.rename(
                from = file.path(path_polynomial, "Config_BaM.txt"),
                to = cf_file
            ))
        }

        if (do_prediction) {
            system2(
                command = file.path(dir_exe_BaM, "BaM"),
                args = c("-cf", cf_file),
                wait = FALSE
            )
            Sys.sleep(0.1)
        }
    }
}


##########################################
### Plots
##########################################
Kmin_literature <- data.frame(
    x_start = c(0),
    x_end = c(16),
    min_value = c(1 / 0.013),
    max_value = c(1 / 0.009),
    mean_value = c(1 / 0.010)
)

Q_observed <- c(120)
final_calibraion <- TRUE
counter <- 1

all_events <- c(
    "Q_120")
factor_scaled <- 1000 # WSE in m all the time

# Do not touche
conf_level <- 0.95 # Hard-coded
z_val <- qnorm(1 - (1 - conf_level) / 2)
# Desired level order
desired_order <- c("Total", "Parametric", "Maxpost", "Observations")

all_prediction_case <- c("ParamU", "Maxpost", "TotalU")

if (do_plots) {
    #############
    # DIC Estimation
    #############
    Plot_dic_rectangular <- plot_DIC(path_experiment) +
        theme(legend.position = "none")

    ggsave(
        plot = Plot_dic_rectangular, filename = file.path(path_experiment, "DIC.png"),
        width = 9, height = 7, units = "cm",
        dpi = 600
    )



    for (n_degree in n_degree_seq) {
        # n_degree <- 0

        all_uncertainty_data <- c()
        # Polynomial degree i
        path_polynomial <- file.path(path_experiment, paste0("n_", n_degree))
        path_temp_plots <- file.path(path_polynomial, "Calibration")
        if (!dir.exists(path_polynomial)) stop(paste0("Polynomial degree ", n_degree, " is not performed yet. Please do the calibration"))
        path_post_traitement <- file.path(path_temp_plots, "post_traitement")
        path_post_traitement_data <- file.path(path_post_traitement, "RData")
        # # Read observation
        # load(file.path(path_post_traitement, "RData", "Data_Z_sim_vs_obs.RData"))

        load(file.path(path_post_traitement_data, "BaM_objects.RData"))

        if (!dir.exists(path_post_traitement)) {
            dir.create(path_post_traitement)
        }

        if (!dir.exists(path_post_traitement_data)) {
            dir.create(path_post_traitement_data)
        }

        path_model_mage <- c(paste0(path_polynomial, "/model_mage/", all_events, "/"))

        mcmc_not_cooked <- readMCMC(file.path(path_temp_plots, "Results_MCMC.txt"))
        plots <- tracePlot(mcmc_not_cooked)

        mcmcplot <- wrap_plots(plots, ncol = 3)
        # Save plot and data
        ggsave(
            file.path(
                path_post_traitement,
                "MCMC_not_cooked.png"
            ),
            mcmcplot,
            width = 20,
            height = 20,
            units = "cm"
        )
        # Density plot for each parameter
        plots <- densityPlot(mcmc_not_cooked)
        pdf_plot <- wrap_plots(plots, ncol = 3)

        ggsave(
            file.path(
                path_post_traitement,
                "densityplot_not_cooked.png"
            ),
            pdf_plot,
            width = 20,
            height = 20,
            units = "cm"
        )

        if (final_calibraion) {
            # zmatrix
            Z_MatrixKmin <- mod$xtra$object$zKmin
            Z_MatrixKmoy <- mod$xtra$object$zKmoy
            # Covariate discretization
            covariate_grid <- read.table(
                file.path(path_temp_plots, "legendre_covariate_non_normalized.txt"),
                header = T
            )
            covariate_grid <- covariate_grid$Covariate

            if (!file.exists(file.path(path_temp_plots, "Results_Cooking.txt"))) stop("MCMC is still running or calculation is not going to the end. Verify if calibration is already finished or verify that calibration has not error messages. Please put final_results = FALSE as input")

            mcmc <- readMCMC(file.path(path_temp_plots, "Results_Cooking.txt"))
            plots <- tracePlot(mcmc)

            mcmcplot <- wrap_plots(plots, ncol = 3)


            ggsave(
                file.path(
                    path_post_traitement,
                    "MCMC_Cooked.png"
                ),
                mcmcplot,
                width = 20,
                height = 20,
                units = "cm"
            )

            # Density plot for each parameter
            plots <- densityPlot(mcmc)
            pdf_plot <- wrap_plots(plots, ncol = 3)

            ggsave(
                file.path(
                    path_post_traitement,
                    "densityplot_Cooked.png"
                ),
                pdf_plot,
                width = 20,
                height = 20,
                units = "cm"
            )


            png(
                file.path(path_post_traitement, "corelation_cooked.png"),
                width = 800,
                height = 800,
                res = 120
            )
            pairs(mcmc)
            dev.off()

            getSummary <- read.table(
                file = file.path(path_temp_plots, "Results_Summary.txt"),
                header = TRUE,
                stringsAsFactors = FALSE
            )

            # Values of error model in meter for WSE, in m3/s for discharge and m/s for velocity.
            knitr::kable(getSummary,
                align = "c"
            )


            # Zoom into the MAP and standard deviation of the error model
            getSummary_zoom <- getSummary[c(11, 16), ]
            getSummary_zoom[, c("Y1_intercept")] <- getSummary_zoom[, c("Y1_intercept")] * factor_scaled # WSE in mm

            # Get MAP simulation
            MAP_param_matrix <- as.numeric(getSummary_zoom[2, c(1:(n_degree + 1))])
        } else {
            stop("Prediction is not conducted yet")
        }

        Kmin_calculations <- k_plot(
            matrix_spatialisation = Z_MatrixKmin,
            mcmc = mcmc,
            n_degree_kmin = n_degree,
            n_degree_kflood = 0,
            covariate_discretization = covariate_grid,
            MAP_param_vector = MAP_param_matrix,
            main_channel = TRUE,
            ks_literature = Kmin_literature
        )
        df_MAP_Kmin <- Kmin_calculations[[1]]
        Kmin_plot <- Kmin_calculations[[2]]
        df_envelope_Kmin <- Kmin_calculations[[3]]

        ls_spatial_friction_Kmin <- list(
            df_envelope = df_envelope_Kmin,
            df_MAP = df_MAP_Kmin
        )
        save(ls_spatial_friction_Kmin,
            file = file.path(path_post_traitement_data, "Data_friction_estimation_ls_spatial_friction_Kmin.RData")
        )
        ggsave(
            file.path(
                path_post_traitement,
                "Kmin.png"
            ),
            Kmin_plot,
            width = 13,
            height = 10,
            units = "cm",
            dpi = 600
        )
        # Kflood_calculations <- K_plot(
        #         matrix_spatialisation = Z_MatrixKflood,
        #         mcmc = mcmc,
        #         n_param_Kmin = length(Kmin_prior),
        #         n_param_Kflood = length(Kflood_prior),
        #         covariate_discretization = covariate_grid,
        #         MAP_param_vector = MAP_param_matrix,
        #         main_channel = FALSE
        #     )
        #     df_MAP_Kflood <- Kflood_calculations[[1]]
        #     Kflood_plot <- Kflood_calculations[[2]]
        #     df_envelope_Kflood <- Kflood_calculations[[3]]

        #     ls_spatial_friction_Kflood <- list(
        #         df_envelope = df_envelope_Kflood,
        #         df_MAP = df_MAP_Kflood
        #     )
        #     save(ls_spatial_friction_Kflood,
        #         file = file.path(path_post_traitement_data, "Data_friction_estimation_ls_spatial_friction_Kflood.RData")
        #     )

        ####################################
        # Residuals
        ###################################
        CalData <- read.table(file.path(path_temp_plots, "CalibrationData.txt"), header = TRUE)
        # Read residuals
        if (final_calibraion) {
            residuals <- read.table(
                file = file.path(path_temp_plots, "Results_Residuals.txt"),
                header = TRUE,
                stringsAsFactors = FALSE
            )

            # Convert residuals to mm and m3/s
            residuals_mm_m3_s <- data.frame(
                residuals[, 1:4],
                residuals[, c(9, 19, 20, 24, 29)] * factor_scaled
            )
            sim_event_mm_m3_s <- data.frame(residuals_mm_m3_s, Yu_z = CalData$Yu_z * factor_scaled)
        } else {
            stop("Calibration is not conducted yet")
        }

        Q_sim_vs_obs <- ZQdX_sim(
            sim_event_mm_m3_s = sim_event_mm_m3_s,
            Q_input = Q_observed,
            Qplot = TRUE,
            title_label = "MAP simulations vs boundary condition \n(Discharge upstream)",
            ylabel = "Q (L/s)"
        )
        ggsave(
            file.path(
                path_post_traitement,
                "sim_vs_obs_Q_MAP.png"
            ),
            Q_sim_vs_obs,
            width = 12,
            height = 10,
            units = "cm"
        )
        save(Q_sim_vs_obs,
            file = file.path(path_post_traitement_data, "Plot_sim_vs_obs_Q_MAP.RData")
        )
        Z_sim_vs_obs <- ZQdX_sim(
            sim_event_mm_m3_s = sim_event_mm_m3_s,
            Q_input = NULL,
            Qplot = FALSE,
            title_label = "MAP simulations vs Observations (WSE)",
            ylabel = "Water surface elevation (mm)"
        )
        ggsave(
            file.path(
                path_post_traitement,
                "sim_vs_obs_WSE_MAP.png"
            ),
            Z_sim_vs_obs,
            width = 12,
            height = 12,
            units = "cm"
        )
        save(Z_sim_vs_obs,
            file = file.path(path_post_traitement_data, "Plot_sim_vs_obs_WSE_MAP.RData")
        )
        save(sim_event_mm_m3_s,
            file = file.path(path_post_traitement_data, "Data_Z_sim_vs_obs.RData")
        )


        #############
        # WSE predictions
        #############

        X3 <- read.table(file.path(path_temp_plots, "X3.pred"))
        X1 <- read.table(file.path(path_temp_plots, "X1.pred"))
        # X2 = read.table(file.path(path_temp_plots,'X2.pred'))

        X_axis <- cbind(X1, X3)
        colnames(X_axis) <- c("Event", "x")
        WSE_TotalU <- read.table(file.path(path_temp_plots, "TotalU_WSE.env"), header = TRUE)
        WSE_TotalU <- WSE_TotalU[, c("q2.5", "q97.5")]

        WSE_ParamU <- read.table(file.path(path_temp_plots, "ParamU_WSE.env"), header = TRUE)
        WSE_ParamU <- WSE_ParamU[, c("q2.5", "q97.5")]

        WSE_MAP <- read.table(file.path(path_temp_plots, "Maxpost_WSE.spag"))
        colnames(WSE_MAP) <- "Maxpost"

        all_uncertainty_data <- rbind(
            data.frame(
                X_axis,
                WSE_TotalU * factor_scaled,
                id = "Total",
                Maxpost = NA
            ),
            data.frame(
                X_axis,
                WSE_ParamU * factor_scaled,
                id = "Parametric",
                Maxpost = NA
            ),
            data.frame(
                X_axis,
                `q2.5` = NA,
                `q97.5` = NA,
                id = "Maxpost",
                Maxpost = WSE_MAP * factor_scaled
            )
        )

        observation_data <- CalData %>%
            select(event, x, Yu_z, WSE) %>%
            mutate(
                "q2.5" = (WSE - z_val * Yu_z) * factor_scaled,
                "q97.5" = (WSE + z_val * Yu_z) * factor_scaled,
                id = "Observations",
                Value = WSE * factor_scaled
            ) %>%
            rename(
                Event = event
            ) %>%
            select(Event, x, q2.5, q97.5, id, Value)

        # Get actual levels present in the data
        present_levels <- intersect(
            desired_order,
            unique(all_uncertainty_data$id)
        )

        # Apply factor levels only for those that exist
        all_uncertainty_data$id <- factor(
            all_uncertainty_data$id,
            levels = present_levels
        )



        # Custom names for facet labels
        facet_names <- c(
            "1" = "Q = 120",
            "2" = "Q = 60",
            "3" = "Q = 30",
            "4" = "Q = 14"
        )
        # Plot simulation with uncertainty and observations
        sim_obs_plot <-
            ggplot(
                data = filter(all_uncertainty_data, !id %in% c("Maxpost", "Prior")),
                aes(
                    x = x,
                    ymin = q2.5,
                    ymax = q97.5
                ), alpha = 0.65
            ) +
            # Residual uncertainty bands as ribbon
            geom_ribbon(
                aes(fill = id),
                alpha = 0.65
            )

        # Check if Maxpost is available
        if (any(all_prediction_case == "Maxpost")) {
            sim_obs_plot <- sim_obs_plot +
                geom_line(
                    data = subset(all_uncertainty_data, !is.na(Maxpost)),
                    aes(y = Maxpost, color = id),
                    linewidth = 0.8
                )
        }
        # Plot simulation with uncertainty and observations
        sim_obs_plot <- sim_obs_plot +
            # Add observations
            geom_point(
                data = observation_data,
                aes(x = x, y = Value, col = id)
            ) +
            geom_errorbar(
                data = observation_data,
                aes(
                    x = x,
                    ymin = q2.5,
                    ymax = q97.5,
                    col = id
                ),
                na.rm = TRUE
            ) +
            facet_wrap(~Event,
                scales = "free_y",
                labeller = as_labeller(facet_names)
            ) +
            theme_bw() +
            labs(
                x = "Streamwise position (meters)",
                y = "Water surface elevation (mm)",
                title = paste0("Comparison of simulated and observed \nwater surface profiles \n"),
                fill = "Uncertainties",
                color = "Data"
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
            ) +
            theme(
                strip.text = element_text(size = 12),
                axis.text.x = element_text(angle = 45, hjust = 1),
                plot.title = element_text(hjust = 0.5),
                legend.title = element_text(hjust = 0.5)
            )

        ggsave(
            file.path(
                path_post_traitement,
                "sim_vs_obs_WSE_with_uncertainties.png"
            ),
            sim_obs_plot,
            width = 15,
            height = 12,
            units = "cm",
            dpi = 600
        )

        save(sim_obs_plot,
            file = file.path(path_post_traitement_data, "output_var_comparison_with_uncertainties_plot.RData")
        )

        ############################
        # Residuals with uncertainty
        ############################
        res_data <-
            observation_data %>%
            left_join(
                all_uncertainty_data %>%
                    filter(!is.na(Maxpost)) %>%
                    select(Event, x, Maxpost), # keep join columns!
                by = c("Event", "x")
            ) %>%
            mutate(
                id_prediction = "Maxpost"
            ) %>%
            mutate(
                residual = Value - Maxpost,
                residual_low = q2.5 - Maxpost,
                residual_high = q97.5 - Maxpost
            ) %>%
            select(Event, x, id_prediction, residual, residual_low, residual_high)


        res_uncertainty <-
            observation_data %>%
            crossing(id_prediction = c("Parametric", "Total")) %>% # créer une ligne pour chaque id_prediction
            left_join(
                all_uncertainty_data %>%
                    select(Event, id, x, q2.5, q97.5) %>%
                    rename(
                        q2.5_sim = q2.5,
                        q97.5_sim = q97.5,
                    ),
                by = c("Event", "id_prediction" = "id", "x")
            ) %>%
            mutate(
                residual_low  = Value - q2.5_sim,
                residual_high = Value - q97.5_sim
            ) %>%
            select(Event, x, id_prediction, residual_low, residual_high)

        # # all_uncertainty_data() %>%
        #     filter(id %in% c("Observations", "Maxpost")) %>%
        #     pivot_wider(
        #         names_from = id,
        #         values_from = c(Value, q2.5, q97.5),
        #         names_sep = "_"
        #     ) %>%
        #     # Compute residuals
        #     mutate(
        #         residual = Value_Observations - Value_Maxpost,
        #         residual_low = q2.5_Observations - Value_Maxpost,
        #         residual_high = q97.5_Observations - Value_Maxpost,
        #         id = "Maxpost"
        #     ) %>%
        #     select(Event, x, id, residual, residual_low, residual_high)

        # res_data <- res_data[!is.na(res_data$residual), ]

        # residual_total <- all_uncertainty_data %>%
        #     filter(id %in% c("Parametric", "Total", "Maxpost")) %>%
        #     pivot_wider(
        #         names_from = id,
        #         values_from = c(Value, q2.5, q97.5),
        #         names_sep = "_"
        #     ) %>%
        #     mutate(
        #         residual_low = Value_Maxpost - q2.5_Total,
        #         residual_high = Value_Maxpost - q97.5_Total,
        #         id = "Total"
        #     ) %>%
        #     select(Event, x, id, residual_low, residual_high)

        # residual_parametric <- all_uncertainty_data %>%
        #     filter(id %in% c("Parametric", "Total", "Maxpost")) %>%
        #     pivot_wider(
        #         names_from = id,
        #         values_from = c(Value, q2.5, q97.5),
        #         names_sep = "_"
        #     ) %>%
        #     mutate(
        #         residual_low = Value_Maxpost - q2.5_Parametric,
        #         residual_high = Value_Maxpost - q97.5_Parametric,
        #         id = "Parametric"
        #     ) %>%
        #     select(Event, x, id, residual_low, residual_high)

        # residuals_all <- bind_rows(
        #     res_data,
        #     residual_total,
        #     residual_parametric
        # )

        residuals_all <- bind_rows(
            res_data,
            res_uncertainty
        )

        # Get actual levels present in the data
        present_levels <- intersect(
            desired_order,
            unique(residuals_all$id)
        )
        # Apply factor levels only for those that exist
        residuals_all$id <- factor(
            residuals_all$id,
            levels = present_levels
        )
        # Set residual plot
        res_sim_obs_plot <-
            ggplot(
                data = filter(residuals_all, !id %in% c("Maxpost", "Prior")),
                aes(
                    x = x,
                    y = residual,
                    ymin = residual_low,
                    ymax = residual_high,
                    fill = id
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
                        filter(id == "Maxpost"),
                    aes(y = residual, color = "Maxpost")
                ) +
                geom_errorbar(
                    data = residuals_all %>%
                        filter(id == "Maxpost"),
                    aes(
                        ymin = residual_low,
                        ymax = residual_high,
                        color = "Maxpost"
                    )
                )
        }

        res_sim_obs_plot <-
            res_sim_obs_plot +
            geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
            facet_wrap(~Event, labeller = as_labeller(facet_names)) +
            theme_bw() +
            scale_fill_manual(
                values = c(
                    "Parametric" = "pink",
                    "Total" = "red"
                )
            ) +
            scale_color_manual(
                values = c(
                    "Maxpost" = "black"
                )
            ) +
            labs(
                title = "Residuals water surface profiles with uncertainties",
                x = "Streamwise position (meters)",
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
            dpi = 600,
            width = 15,
            height = 12,
            units = "cm"
        )
        save(res_sim_obs_plot,
            file = file.path(path_post_traitement_data, "residual_output_var_comparison_with_uncertainties_plot.RData")
        )
    }
}
