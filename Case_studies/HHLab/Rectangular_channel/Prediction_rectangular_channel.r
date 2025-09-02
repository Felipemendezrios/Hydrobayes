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

# Set directory
dir_workspace <- here::here()

do_prediction <- FALSE
do_plots <- TRUE

prediction_file <- c("Prior", "ParamU", "Maxpost", "TotalU")
n_prediction <- length(prediction_file)
nsim_prior <- 500 # number of simulation for propagation (def : 100 for priori prediction)

n_degree_max <- 3
Experiment_id <- "4_WSE_main_channel_real_uncertainty"
path_experiment <- file.path("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/Case_studies/HHLab/Rectangular_channel/Calibration_experiments", Experiment_id)
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
    X <- data$data[, data$col.X]
    grid <- data.frame(grid = Caldata[, c("x")])

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



### Plots (not finished yet!)

# Do not touche
conf_level <- 0.95 # Hard-coded
z_val <- qnorm(1 - (1 - conf_level) / 2)
# Desired level order
desired_order <- c("TotalU", "ParamU", "Maxpost", "Prior")

if (do_plots) {
    # for (n_degree in n_degree_seq) {}

    n_degree <- 0

    all_uncertainty_data <- c()
    # Polynomial degree i
    path_polynomial <- file.path(path_experiment, paste0("n_", n_degree), "Calibration")
    if (!dir.exists(path_polynomial)) stop(paste0("Polynomial degree ", n_degree, " is not performed yet. Please do the calibration"))
    path_post_traitement <- file.path(path_polynomial, "post_traitement")
    path_post_traitement_data <- file.path(path_post_traitement, "RData")
    # Read observation
    load(file.path(path_post_traitement, "RData", "Data_Z_sim_vs_obs.RData"))
}
