rm(list = ls())
graphics.off()

# function_list <- list.files("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/Functions/", full.names = TRUE)
# for (i in function_list) {
#     source(i)
# }


# Libraries
library(RBaM)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(stringr)

# do_prediction <- TRUE
do_prediction <- FALSE
do_plots <- TRUE

prediction_file <- c("Prior", "ParamU", "Maxpost", "TotalU")
n_prediction <- length(prediction_file)
nsim_prior <- 500 # number of simulation for propagation (def : 100 for priori prediction)

Experiment_id <- "Syn_case_without_ks_measures_all_obs_flat_prior"
# "Syn_case_with_manual_Ks_obs_limites_all_obs_flat_prior"
# "Syn_case_rec_5_3_with_Ks_obs_limites_all_obs_flat_prior"
# "Syn_case_without_ks_measures_all_obs_flat_prior"

path_experiment <- file.path("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/Case_studies/Synthetic_case/Simplified_ks_rectangular_MC/Calibration_experiments", Experiment_id)

# Experiments input data to be used during calibration setting

all_cal_case <- c(
    "Kmin_n_0_Kflood_n_0_flat.r",
    "Kmin_n_1_Kflood_n_0_flat.r",
    "Kmin_n_2_Kflood_n_0_flat.r",
    "Kmin_n_3_Kflood_n_0_flat.r",
    "Kmin_n_4_Kflood_n_0_flat.r"
    # "Kmin_n_0_Kflood_n_0.r",
    # "Kmin_n_1_Kflood_n_0.r",
    # "Kmin_n_2_Kflood_n_0.r",
    # "Kmin_n_3_Kflood_n_0.r",
    # "Kmin_n_4_Kflood_n_0.r"
)

# Define the calculation grid: Extract from net folder and ST File

# Events and number of reaches in total
info_events_reaches <- list(
    event_1 = list(
        data.frame(
            event = c(1),
            reach = c(1, 3),
            xmin = c(0, 18),
            xmax = c(18, 25)
        )
    ),
    event_2 = list(
        data.frame(
            event = c(2),
            reach = c(2),
            xmin = (0),
            xmax = c(20)
        )
    )
)

for (ev in names(info_events_reaches)) {
    df <- info_events_reaches[[ev]][[1]]

    nb_reaches <- nrow(df)

    if (nb_reaches != length(df$xmin) || nb_reaches != length(df$xmax)) {
        stop(paste(
            "Erreur dans", ev,
            ": nb_reaches doit être égal au nombre de xmin et xmax"
        ))
    }
}

grid_user <- do.call(
    rbind,
    lapply(info_events_reaches, function(x) x[[1]])
)
rownames(grid_user) <- NULL

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

X_temp$t <- 10800 # any number to complet the input variables to keep the format

# Do not touch
dir_exe_BaM <- "/home/famendezrios/Documents/Git/BaM/makefile/"
nY <- 5
###################################################

for (id_cal_case in 1:length(all_cal_case)) {
    # Intitialization
    pred_confi_file_name <- c() # vector with setting prediction file name (Config_Pred_"type of prediction")
    doStructural_logical <- c() # logical vector with setting in function of type of prediction
    doParametric_logical <- c() # logical vector with setting in function of type of prediction
    priorNsim_int <- c() # list of number of simulation for
    pred_list <- c() # list with information about all predictions
    pred_var_name <- c() # list with names of each variable
    mod_list <- list() # list to write several models for running in parallel


    path_polynomial <- file.path(
        path_experiment,
        sub("\\.r$", "", all_cal_case[id_cal_case])
    )
    if (!dir.exists(path_polynomial)) stop(paste0("Case ", all_cal_case[id_cal_case], " is not performed yet. Please do the calibration"))

    if (!file.exists(file.path(
        path_polynomial, "BaM",
        "Results_Cooking.txt"
    ))) {
        stop("Calibration is not performed yet")
    }

    # Path to save results
    workspace_user <- file.path(path_polynomial, "BaM")
    path_post_traitement <- file.path(workspace_user, "post_traitement")
    path_post_traitement_data <- file.path(path_post_traitement, "RData")

    # Load data and model used during calibration
    load(file.path(path_post_traitement_data, "BaM_objects.RData"))

    Caldata <- data$data
    X <- full_join(X_temp, Caldata[, 1:4], by = c("event", "reach", "x", "t")) %>%
        arrange(event, reach, x)

    names_file_prediction <- colnames(Caldata[, !grepl(colnames(Caldata), pattern = "^Yu_")])[c(-1:-4)]
    if (length(names_file_prediction) > 5) {
        names_file_prediction <- names_file_prediction[1:5]
    }

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
            workspace = workspace_user,
            dir.exe = dir_exe_BaM,
            name.exe = "BaM",
            predMaster_fname = paste0("Config_Pred_Master_", prediction_file[i], ".txt")
        )

        # Rename Config_BaM.txt file to the result folder
        cf_file <- file.path(workspace_user, paste0("Config_BaM_", prediction_file[i], ".txt"))
        if (file.exists(file.path(workspace_user, "Config_BaM.txt"))) {
            invisible(file.rename(
                from = file.path(workspace_user, "Config_BaM.txt"),
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

final_calibraion <- TRUE
counter <- 1


all_events <- c(
    "5_1",
    "3_2"
)

factor_scaled <- 1 # WSE in m all the time

# Do not touche
conf_level <- 0.95 # Hard-coded
z_val <- qnorm(1 - (1 - conf_level) / 2)
# Desired level order
desired_order <- c("Total", "Parametric", "Maxpost", "Observations")

all_prediction_case <- c("ParamU", "Maxpost", "TotalU")

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

## Functions:
# Plot spatialisation of k in main channel or floodplain
K_plot <- function(
    matrix_spatialisation,
    mcmc,
    n_param_Kmin,
    n_param_Kflood,
    SR_reaches,
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

    k_estimated_MCMC$KP <- SR_reaches$KP
    k_estimated_MCMC$reaches_nb <- SR_reaches$reach
    k_estimated_MCMC$reaches_sr <- SR_reaches$id_river

    # Convert to long format
    df_MCMC_sampling <- tidyr::pivot_longer(
        k_estimated_MCMC,
        cols = -c(reaches_nb, reaches_sr, KP),
        values_to = "Value"
    ) %>%
        select(reaches_nb, reaches_sr, KP, Value) %>%
        mutate(ID = "MCMC Sampling")

    k_estimated_MAP <- as.matrix(matrix_spatialisation) %*% MAP_param_vector[indx]

    df_MAP <- data.frame(
        reaches_nb = SR_reaches$reach,
        reaches_sr = SR_reaches$id_river,
        KP = SR_reaches$KP,
        Value = k_estimated_MAP,
        ID = "MAP"
    )

    # Get 95% uncertainty for envelope curve : create ribbon data from MCMC
    df_envelope <- df_MCMC_sampling %>%
        filter(ID == "MCMC Sampling") %>%
        group_by(KP, reaches_sr, reaches_nb) %>%
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
        facet_wrap(~reaches_sr, scales = "free", ncol = 1)


    return(list(df_MAP, K_plot, df_envelope))
}


get_init_prior <- function(parameter, FIX_dist = FALSE) {
    # Identify if parameter is remnantErrorModel class
    logical_test <- class(parameter[[1]]) == "remnantErrorModel"
    if (logical_test) { # if remnantErrorModel, a list is needed
        init_priors <- list()
    } else { # if not a vector is needed
        init_priors <- numeric(0)
    }

    counter_gamma <- 1
    for (i in parameter) {
        # Handle if parameters is remnantErrorModel
        if (logical_test) {
            param <- i$par
            number_var_error_model <- seq_along(param)

            for (local_counter in number_var_error_model) {
                if (FIX_dist) {
                    init_priors[[counter_gamma]] <- param[[local_counter]]$init
                    counter_gamma <- counter_gamma + 1
                } else {
                    if (param[[local_counter]]$prior$dist != "FIX") {
                        init_priors[[counter_gamma]] <- param[[local_counter]]$init
                        counter_gamma <- counter_gamma + 1
                    }
                }
            }
        } else {
            # Handle if parameters comes from theta
            param <- i
            if (FIX_dist) {
                init_priors <- c(init_priors, param$init)
            } else {
                if (param$prior$dist != "FIX") {
                    init_priors <- c(init_priors, param$init)
                }
            }
        }
    }
    return(init_priors)
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
            facet_wrap(~X1_obs, ncol = 1)

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
            facet_wrap(~X1_obs, ncol = 1)
    }
    return(ZQplot)
}


get_init_prior(mod$par)
data_n_degree <- data.frame(all_cal_case,
    n_degree = c(0, 1, 2, 3, 4)
)

if (do_plots) {
    for (id_cal_case in 1:length(all_cal_case)) {
        all_uncertainty_data <- c()


        path_polynomial <- file.path(
            path_experiment,
            sub("\\.r$", "", all_cal_case[id_cal_case])
        )
        if (!dir.exists(path_polynomial)) stop(paste0("Polynomial degree ", n_degree, " is not performed yet. Please do the calibration"))

        if (!file.exists(file.path(
            path_polynomial, "BaM",
            "Results_Cooking.txt"
        ))) {
            stop("Calibration is not performed yet")
        }

        # Path to save results
        workspace_user <- file.path(path_polynomial, "BaM")
        path_post_traitement <- file.path(workspace_user, "post_traitement")
        path_post_traitement_data <- file.path(path_post_traitement, "RData")

        # # Read observation
        load(file.path(path_post_traitement_data, "BaM_objects.RData"))

        if (!dir.exists(path_post_traitement)) {
            dir.create(path_post_traitement)
        }

        if (!dir.exists(path_post_traitement_data)) {
            dir.create(path_post_traitement_data)
        }

        path_model_mage <- c(paste0(path_polynomial, "/model_mage/", all_events, "/"))


        if (final_calibraion) {
            # zmatrix
            Z_MatrixKmin <- mod$xtra$object$zKmin
            Z_MatrixKmoy <- mod$xtra$object$zKmoy
            # Covariate discretization
            covariate_grid <- read.table(
                file.path(workspace_user, "grid_covariate_non_normalized.txt"),
                header = T
            )
            covariate_grid <- covariate_grid$Covariate

            if (!file.exists(file.path(workspace_user, "Results_Cooking.txt"))) stop("MCMC is still running or calculation is not going to the end. Verify if calibration is already finished or verify that calibration has not error messages. Please put final_results = FALSE as input")

            getSummary <- read.table(
                file = file.path(workspace_user, "Results_Summary.txt"),
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
            MAP_param_matrix <- as.numeric(getSummary_zoom[2, c(1:(data_n_degree$n_degree[id_cal_case] + 1))])
        } else {
            stop("Prediction is supported only when final results are already obtained")
        }
        ####################################
        # Residuals
        ###################################
        CalData <- read.table(file.path(workspace_user, "CalibrationData.txt"), header = TRUE)
        # Read residuals
        if (final_calibraion) {
            residuals <- read.table(
                file = file.path(workspace_user, "Results_Residuals.txt"),
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
            stop("Calibration step are not conducted yet or it is in process")
        }
        #############
        # WSE predictions
        #############

        X3 <- read.table(file.path(workspace_user, "X3.pred"))
        X1 <- read.table(file.path(workspace_user, "X1.pred"))
        # X2 = read.table(file.path(workspace_user,'X2.pred'))

        X_axis <- cbind(X1, X3)
        colnames(X_axis) <- c("Event", "x")
        WSE_TotalU <- read.table(file.path(workspace_user, "TotalU_WSE.env"), header = TRUE)
        WSE_TotalU <- WSE_TotalU[, c("q2.5", "q97.5")]

        WSE_ParamU <- read.table(file.path(workspace_user, "ParamU_WSE.env"), header = TRUE)
        WSE_ParamU <- WSE_ParamU[, c("q2.5", "q97.5")]

        WSE_MAP <- read.table(file.path(workspace_user, "Maxpost_WSE.spag"))
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
            "1" = "Q(SU 1) = 5 m3/s",
            "2" = "Q(SU 2) = 2 m3/s",
            "Upstream" = "Upstream",
            "Downstream" = "Downstream",
            "Tributary" = "Tributary"
        )
        all_uncertainty_data$id_mage <- NA

        idx <- which(all_uncertainty_data$x >= 18 & all_uncertainty_data$Event == 1)
        all_uncertainty_data$id_mage[idx] <- "Downstream"
        idx <- which(all_uncertainty_data$x <= 18 & all_uncertainty_data$Event == 1)
        all_uncertainty_data$id_mage[idx] <- "Upstream"
        idx <- which(all_uncertainty_data$Event == 2)
        all_uncertainty_data$id_mage[idx] <- "Tributary"

        all_uncertainty_data$Event <- factor(
            all_uncertainty_data$Event,
            levels = c("1", "2")
        )
        all_uncertainty_data$id_mage <- factor(
            all_uncertainty_data$id_mage,
            levels = c("Upstream", "Downstream", "Tributary")
        )
        observation_data$id_mage <- NA
        idx <- which(observation_data$x >= 18 & observation_data$Event == 1)
        observation_data$id_mage[idx] <- "Downstream"
        idx <- which(observation_data$x <= 18 & observation_data$Event == 1)
        observation_data$id_mage[idx] <- "Upstream"
        idx <- which(observation_data$Event == 2)
        observation_data$id_mage[idx] <- "Tributary"

        observation_data$Event <- factor(
            observation_data$Event,
            levels = c("1", "2")
        )
        observation_data$id_mage <- factor(
            observation_data$id_mage,
            levels = c("Upstream", "Downstream", "Tributary")
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
                    aes(y = Maxpost, color = id)
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
                scales = "free",
                labeller = as_labeller(facet_names),
                ncol = 1
            ) +
            theme_bw() +
            labs(
                x = "Streamwise position (meters)",
                y = "Water surface elevation (meters)",
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
            width = 17,
            height = 22,
            units = "cm",
            dpi = 300
        )

        save(sim_obs_plot,
            file = file.path(path_post_traitement_data, "output_var_comparison_with_uncertainties_plot.RData")
        )

        sim_obs_plot_by_reaches <- sim_obs_plot +
            facet_wrap(
                ~ id_mage + Event,
                scales = "free",
                labeller = labeller(Event = as_labeller(facet_names)),
                ncol = 1
            )


        ggsave(
            file.path(
                path_post_traitement,
                "sim_vs_obs_WSE_with_uncertainties_by_reaches.png"
            ),
            sim_obs_plot_by_reaches,
            width = 17,
            height = 26,
            units = "cm",
            dpi = 300
        )

        save(sim_obs_plot_by_reaches,
            file = file.path(path_post_traitement_data, "output_var_comparison_with_uncertainties_plot_by_reaches.RData")
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
            facet_wrap(~Event,
                labeller = as_labeller(facet_names), scales = "free_x",
                ncol = 1
            ) +
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
                y = "Residuals (meters)",
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
            dpi = 300,
            width = 15,
            height = 12,
            units = "cm"
        )
        save(res_sim_obs_plot,
            file = file.path(path_post_traitement_data, "residual_output_var_comparison_with_uncertainties_plot.RData")
        )
    }
}
