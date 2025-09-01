rm(list = ls())
graphics.off()

function_list <- list.files("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/Functions/", full.names = TRUE)
for (i in function_list) {
    source(i)
}


# Libraries
library(RBaM)
library(ggplot2)

# Set directory
dir_workspace <- here::here()

do_prediction <- TRUE

prediction_file <- c("Prior", "ParamU", "Maxpost", "TotalU")
n_prediction <- length(prediction_file)
nsim_prior <- 500 # number of simulation for propagation (def : 100 for priori prediction)

n_degree_max <- 3
Experiment_id <- "4_WSE_main_channel_real_uncertainty"
path_experiment <- file.path("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/Case_studies/HHLab/Rectangular_channel/Calibration_experiments", Experiment_id)
# Do not touch
n_degree_seq <- seq(0, n_degree_max, 1)
dir_exe_BaM <- "/home/famendezrios/Documents/Git/BaM/makefile/"
###################################################

# Intitialization
pred_confi_file_name <- c() # vector with setting prediction file name (Config_Pred_"type of prediction")
doStructural_logical <- c() # logical vector with setting in function of type of prediction
doParametric_logical <- c() # logical vector with setting in function of type of prediction
priorNsim_int <- c() # list of number of simulation for
pred_list <- c() # list with information about all predictions
pred_var_name <- c() # list with names of each variable
mod_list <- list() # list to write several models for running in parallel


for (n_degree in n_degree_seq) {
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
    path_BaM_object <- file.path(path_polynomial, "BaM_object_calibration")
    load(file.path(path_BaM_object, "BaM_objects.RData"))

    Caldata <- data$data

    grid <- data.frame(grid = Caldata[, 1])

    names_file_prediction <- colnames(Caldata[, !grepl(colnames(Caldata), pattern = "^Yu_")])[-1]
    nY <- length(names_file_prediction)
    for (i in 1:nY) {
        # name of prediction files (.spag)
        pred_var_name_temps <- c()
        for (j in 1:n_prediction) {
            pred_var_name_temps <- cbind(
                pred_var_name_temps,
                paste0(
                    prediction_file[j],
                    names_file_prediction[i], ".spag"
                )
            )

            pred_confi_file_name[j] <- c(paste0(
                "Config_Pred_",
                prediction_file[j],
                ".txt"
            ))
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
            dir_copied <- file.path(dirname(dir_origin), paste0(basename(dir_origin), "_", prediction_file[j]))

            if (!dir.exists(dir_copied)) {
                dir.create(dir_copied)
            } else {
                file.remove(list.files(dir_copied, full.names = TRUE, recursive = TRUE))
            }

            files_to_copy <- list.files(dir_origin, full.names = TRUE, recursive = TRUE)
            # Relative paths (to preserve directory structure)
            relative_paths <- list.files(dir_origin, full.names = FALSE, recursive = TRUE)
            # Ensure destination directories exist
            subfolder <- dirname(relative_paths)[dirname(relative_paths) != "."]
            subfolderDir <- file.path(dir_copied, subfolder)
            if (!dir.exists(subfolderDir)) {
                dir.create(subfolderDir, recursive = TRUE)
            } else {
                file.remove(list.files(subfolderDir, full.names = TRUE))
            }
            file.copy(
                from = files_to_copy,
                to = file.path(dir_copied, relative_paths),
                overwrite = TRUE
            )
            mod_list[[j]]$xtra$object$mageDir <- paste0(dir_copied, "/")
            mod_list[[j]]$xtra$fname <- sub("\\.txt$", paste0("_", prediction_file[j], ".txt"), mod$xtra$fname)
        }
        pred_var_name <- rbind(pred_var_name, pred_var_name_temps)
    }
    # create a pred_list only for spagfile and another for fname
    for (i in 1:n_prediction) {
        pred_list[[i]] <- prediction(
            X = grid,
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
        # Move Config_BaM.txt file to the result folder
        file.copy(from = dir_cf, to = file.path(
            path_polynomial,
            paste0("Config_BaM_", prediction_file[i], ".txt")
        ), overwrite = TRUE)

        if (do_prediction) {
            system2(
                command = file.path(dir_exe_BaM, "BaM"),
                args = c("-cf", file.path(
                    workspace_user,
                    paste0("Config_BaM_", prediction_file[i], ".txt")
                )),
                wait = FALSE
            )
            Sys.sleep(0.1)
        }
    }
}
