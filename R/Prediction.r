grid_user <- function(
    info_events_reaches) {
    grid_user <- c()
    # Loop through each event
    for (event in info_events_reaches) {
        # Initialize an empty list to store all data frames
        all_dfs <- list()

        # Loop through each element in the event (skip the 'type' field)
        for (i in 2:length(event)) {
            all_dfs <- c(all_dfs, list(event[[i]]))
        }
        # Checks
        check_type_pred(event$type)

        # Combine all SU into one event
        df_event <- do.call(rbind, all_dfs)
        rownames(df_event) <- NULL

        if (event$type == "ZdX") {
            check_dX(df_event)

            # Generate discretization for each row and expand into a long-format data frame
            grid_temp <- do.call(rbind, lapply(1:nrow(df_event), function(i) {
                data.frame(
                    event = df_event$event[i],
                    reach = df_event$reach[i],
                    x = seq(from = df_event$xmin[i], to = df_event$xmax[i], length.out = df_event$nb_discretization[i]),
                    t = df_event$tmin[i]
                )
            }))
        } else if (event$type == "QdXT") {
            check_dXT(df_event)

            # Generate discretization for each row and expand into a long-format data frame
            grid_temp <- do.call(rbind, lapply(1:nrow(df_event), function(i) {
                data.frame(
                    event = df_event$event[i],
                    reach = df_event$reach[i],
                    x = df_event$xmin[i],
                    t = df_event$tmin[i]
                )
            }))
        } else if (event$type == "ZdT") {
            check_dT(df_event)

            # Generate discretization for each row and expand into a long-format data frame
            grid_temp <- do.call(rbind, lapply(1:nrow(df_event), function(i) {
                data.frame(
                    event = df_event$event[i],
                    reach = df_event$reach[i],
                    x = df_event$xmin[i],
                    t = seq(from = df_event$tmin[i], to = df_event$tmin[i], length.out = df_event$nb_discretization[i]),
                )
            }))
        } else {
            stop("Never arrive here, if not, a check must be created before")
        }
        grid_user <- rbind(grid_user, grid_temp)
    }
    return(grid_user = grid_user)
}

# Function to initialize variables for each calibration case
initialize_variables <- function() {
    # Initialization
    doStructural_logical <- c()
    doParametric_logical <- c()
    priorNsim_int <- c()
    pred_var_name <- c()
    mod_list <- list()
    list(
        pred_confi_file_name = character(), # character vector with the prediction file name
        doStructural_logical = list(), # logical vector with setting in function of type of prediction
        doParametric_logical = list(), # logical vector with setting in function of type of prediction
        priorNsim_int = integer(), # vector of number of simulation for prior simulation
        pred_list = list(), # list with information about all predictions
        pred_var_name = character(), # character data frame with names of each variable
        mod_list = list() # list to write several models for running in parallel
    )
}

# Function to load data and prepare X
Add_calData_grid_user <- function(
    X_pred,
    Caldata,
    nX,
    nY) {
    X <- full_join(X_pred, Caldata[, 1:nX], by = c("event", "reach", "x", "t")) %>% arrange(event, reach, x)

    names_file_prediction <- colnames(Caldata[, (nX + 1):(nX + nY)])

    list(X = X, names_file_prediction = names_file_prediction)
}


# Function to set prediction file names and logical flags
set_prediction_configs <- function(prediction_file, nY, nsim_prior = 500) {
    pred_confi_file_name <- sapply(prediction_file, function(j) paste0("Config_Pred_", j, ".txt"))
    doStructural_logical <- lapply(prediction_file, function(j) {
        switch(j,
            "Prior" = rep(FALSE, nY),
            "ParamU" = rep(FALSE, nY),
            "Maxpost" = rep(FALSE, nY),
            "TotalU" = rep(TRUE, nY)
        )
    })
    doParametric_logical <- sapply(prediction_file, function(j) {
        switch(j,
            "Prior" = TRUE,
            "ParamU" = TRUE,
            "Maxpost" = FALSE,
            "TotalU" = TRUE
        )
    })
    priorNsim_int <- sapply(prediction_file, function(j) {
        switch(j,
            "Prior" = nsim_prior,
            "ParamU" = -1,
            "Maxpost" = -1,
            "TotalU" = -1
        )
    })
    list(
        pred_confi_file_name = pred_confi_file_name,
        doStructural_logical = doStructural_logical,
        doParametric_logical = doParametric_logical,
        priorNsim_int = priorNsim_int
    )
}

# Function to copy model directory for parallel runs
copy_model_directory <- function(mod, prediction_file, idx) {
    dir_origin <- mod$xtra$object$mageDir
    dir_copied_all <- file.path(normalizePath(dirname(dir_origin)), prediction_file[idx], basename(dir_origin))
    for (source_dir in dir_origin) {
        source_dir <- normalizePath(source_dir)
        dir_name <- basename(source_dir)
        dir_copied <- file.path(dirname(source_dir), prediction_file[idx], dir_name)

        if (dir.exists(dir_copied)) {
            unlink(dir_copied, recursive = TRUE)
        }
        dir.create(dir_copied, recursive = TRUE)

        files_to_copy <- list.files(source_dir, full.names = TRUE, recursive = TRUE, include.dirs = TRUE)
        for (file in files_to_copy) {
            relative_path <- sub(paste0("^", source_dir, "/?"), "", file)
            dest_path <- file.path(dir_copied, relative_path)
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
    mod$xtra$object$mageDir <- paste0(dir_copied_all, "/")
    mod$xtra$fname <- sub("\\.txt$", paste0("_", prediction_file[idx], ".txt"), mod$xtra$fname)
    mod
}

prediction_MAGE <- function(
    cal_case,
    paths,
    prediction_file,
    data,
    do_prediction,
    X_pred,
    mod,
    remant_error_list,
    mcmcOptions = RBaM::mcmcOptions(),
    mcmcCooking = RBaM::mcmcCooking(),
    mcmcSummary = RBaM::mcmcSummary(xtendedMCMC.fname = "Results_xtendedMCMC.txt"),
    nsim_prior = 500) {
    # Check calibration
    check_calibration_case(paths$path_polynomial)

    # Initialize variables
    vars <- initialize_variables()
    n_prediction <- length(prediction_file)

    nX <- mod$nX
    nY <- mod$nY

    # Add calibration data to grid
    data_to_add_grid <- Add_calData_grid_user(
        X_pred = X_pred,
        Caldata = data$data,
        nX = nX,
        nY = nY
    )
    X <- data_to_add_grid$X
    names_file_prediction <- data_to_add_grid$names_file_prediction

    # Set prediction configs
    configs <- set_prediction_configs(
        prediction_file = prediction_file,
        nY = nY,
        nsim_prior = nsim_prior
    )

    # Prepare prediction file names and copy model directories
    for (i in 1:nY) {
        pred_var_name_temps <- matrix(nrow = 1, ncol = n_prediction)
        for (j in 1:n_prediction) {
            pred_var_name_temps[1, j] <- paste0(prediction_file[j], "_", names_file_prediction[i], ".spag")
            vars$mod_list[[j]] <- copy_model_directory(
                mod = mod,
                prediction_file = prediction_file,
                idx = j
            )
        }
        vars$pred_var_name <- rbind(vars$pred_var_name, pred_var_name_temps)
    }

    cf_file <- c()
    for (i in 1:n_prediction) {
        pred_list <- prediction(
            X = X,
            spagFiles = vars$pred_var_name[, i],
            data.dir = paths$path_BaM_folder,
            fname = configs$pred_confi_file_name[i],
            priorNsim = configs$priorNsim_int[i],
            doParametric = configs$doParametric_logical[i],
            doStructural = configs$doStructural_logical[[i]],
        )

        BaM(
            mod = vars$mod_list[[i]],
            data = data,
            remnant = remant_error_list,
            mcmc = mcmcOptions,
            cook = mcmcCooking,
            summary = mcmcSummary,
            residuals = RBaM::residualOptions(),
            pred = pred_list,
            doCalib = FALSE,
            doPred = TRUE,
            na.value = -9999,
            run = FALSE,
            preClean = FALSE,
            workspace = paths$path_BaM_folder,
            predMaster_fname = paste0("Config_Pred_Master_", prediction_file[i], ".txt")
        )
        cf_file[i] <- file.path(
            paths$path_BaM_folder,
            paste0(
                "Config_BaM_",
                prediction_file[i],
                ".txt"
            )
        )
        template_file <- file.path(
            paths$path_BaM_folder,
            "Config_BaM.txt"
        )

        if (file.exists(template_file)) {
            invisible(file.copy(
                from = template_file,
                to = cf_file[i],
                overwrite = TRUE
            ))
        }
    }
    return(
        list(
            X_pred_grid = X,
            cf_file = cf_file
        )
    )
}
