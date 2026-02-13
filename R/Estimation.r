#### pre-processing to obtain the Calibration data

assign_calibration_and_validation_data <- function(
    Input_Typology,
    Input_Model_Reach,
    CalData,
    All_observations) {
    # 1. Calibration set
    CalData_event <- c()
    for (i in 1:length(Input_Typology)) {
        # Remove the WSE at the nodes
        nodes_extraction <- Input_Model_Reach %>%
            filter(reach %in%
                Input_Typology[[i]]) %>%
            select(-reach)

        nodes_to_remove_cal <- unique(c(nodes_extraction$KP_start, nodes_extraction$KP_end))

        if (nrow(CalData %>%
            filter(id_reach_CAL %in% Input_Typology[[i]])) > 0) {
            if (i == 1) {
                CalData_event <- CalData %>%
                    mutate(Reach_groupped_Cal = case_when(
                        id_reach_CAL %in% Input_Typology[[i]] ~ names(Input_Typology)[i],
                        TRUE ~ Reach_groupped_Cal
                    )) %>%
                    group_by(id_case, Reach_groupped_Cal, id_reach_CAL) %>%
                    # If KP duplicated, take it off
                    filter(!KP %in% nodes_to_remove_cal) %>%
                    ungroup()
            } else {
                CalData_event_temp <- CalData %>%
                    mutate(Reach_groupped_Cal = case_when(
                        id_reach_CAL %in% Input_Typology[[i]] ~ names(Input_Typology)[i],
                        TRUE ~ Reach_groupped_Cal
                    )) %>%
                    group_by(id_case, Reach_groupped_Cal, id_reach_CAL) %>%
                    # If KP duplicated, take it off
                    filter(!KP %in% nodes_to_remove_cal) %>%
                    ungroup()
                CalData_event <- rbind(CalData_event, CalData_event_temp)
            }
        }
    }

    # 2. Validation set
    ValData_event <- All_observations %>%
        anti_join(CalData_event, by = colnames(All_observations)) %>%
        mutate(set = "validation")

    # 3. Combine into a single data frame
    CalValData <- bind_rows(CalData_event, ValData_event) %>% arrange(KP)

    return(CalValData)
}


constructor_CalData <- function(
    observed_data,
    do_manual_uncertainty = FALSE,
    sd_WSE_fixed = 0.05, # in meters
    sd_Q_fixed = 8, # in %
    sd_V_fixed = 3, # in %
    sd_Kmin_fixed = 5, # in m1/3/s
    sd_Kflood_fixed = 5 # in m1/3/s
    ) {
    size_df <- nrow(observed_data)
    # ---------------- WSE ----------------
    if ("WSE" %in% colnames(observed_data)) {
        Y <- data.frame(WSE = observed_data$WSE)

        if (do_manual_uncertainty) {
            Yu <- data.frame(Yu_WSE = rep(sd_WSE_fixed, size_df))
        } else {
            Yu <- data.frame(Yu_WSE = observed_data$Yu_WSE)
        }
    } else {
        Y <- data.frame(WSE = rep(-9999, size_df))
        Yu <- data.frame(Yu_WSE = rep(-9999, size_df))
    }

    # ---------------- Q ----------------
    if ("Q" %in% colnames(observed_data)) {
        Y$Q <- observed_data$Q

        if (do_manual_uncertainty) {
            Yu$Yu_Q <- rep(sd_Q_fixed, size_df)
        } else {
            Yu$Yu_Q <- observed_data$Yu_Q
        }
    } else {
        Y$Q <- rep(-9999, size_df)
        Yu$Yu_Q <- rep(-9999, size_df)
    }

    # ---------------- V ----------------
    if ("V" %in% colnames(observed_data)) {
        Y$V <- observed_data$V

        if (do_manual_uncertainty) {
            Yu$Yu_V <- rep(sd_V_fixed, size_df)
        } else {
            Yu$Yu_V <- observed_data$Yu_V
        }
    } else {
        Y$V <- rep(-9999, size_df)
        Yu$Yu_V <- rep(-9999, size_df)
    }

    # ---------------- Kmin ----------------
    if ("Kmin" %in% colnames(observed_data)) {
        Y$Kmin <- observed_data$Kmin

        if (do_manual_uncertainty) {
            Yu$Yu_Kmin <- rep(sd_Kmin_fixed, size_df)
        } else {
            Yu$Yu_Kmin <- observed_data$Yu_Kmin
        }
    } else {
        Y$Kmin <- rep(-9999, size_df)
        Yu$Yu_Kmin <- rep(-9999, size_df)
    }

    # ---------------- Kflood ----------------
    if ("Kflood" %in% colnames(observed_data)) {
        Y$Kflood <- observed_data$Kflood

        if (do_manual_uncertainty) {
            Yu$Yu_Kflood <- rep(sd_Kflood_fixed, size_df)
        } else {
            Yu$Yu_Kflood <- observed_data$Yu_Kflood
        }
    } else {
        Y$Kflood <- rep(-9999, size_df)
        Yu$Yu_Kflood <- rep(-9999, size_df)
    }

    return(list(
        Y = Y,
        Yu = Yu
    ))
}


get_Key_Info_Typology_Model_Reach <- function(
    Input_Typology,
    Input_Model_Reach,
    total_points_discretization = 100) {
    # At this level, all information is already available to link BaM and HM
    # Objective: common information to share in BaM and HM environment
    # This step is so important to ensure the transmission of information

    # Define the grid of interpolation by Typology
    Key_Info_Typology_Model_Reach <- vector(mode = "list", length = length(Input_Typology))
    names(Key_Info_Typology_Model_Reach) <- names(Input_Typology)
    for (i in seq_along(Key_Info_Typology_Model_Reach)) {
        # Get information from hydraulic model by Typology
        mask_Model_Reach_by_Typology <- Input_Model_Reach[match(Input_Typology[[i]], Input_Model_Reach$reach), ]

        if (nrow(mask_Model_Reach_by_Typology) != 1) {
            if (all(diff(mask_Model_Reach_by_Typology$KP_start) < 0)) {
                Logical_decreasing <- TRUE

                mask_Model_Reach_by_Typology <- mask_Model_Reach_by_Typology %>%
                    rename(
                        KP_start = KP_end,
                        KP_end = KP_start
                    )
            } else if (all(diff(mask_Model_Reach_by_Typology$KP_start) > 0)) {
                Logical_decreasing <- FALSE
            } else {
                stop("All KP_start values must be increased or decreased")
            }
        } else if (nrow(mask_Model_Reach_by_Typology) == 0) {
            stop("Any matched found between MR and XR, please verify")
        } else {
            if ((mask_Model_Reach_by_Typology$KP_start - mask_Model_Reach_by_Typology$KP_end) < 0) {
                Logical_decreasing <- FALSE
            } else {
                Logical_decreasing <- TRUE

                mask_Model_Reach_by_Typology <- mask_Model_Reach_by_Typology %>%
                    rename(
                        KP_start = KP_end,
                        KP_end = KP_start
                    )
            }
        }

        # Get nodes of Model_Reach
        Model_Reach_nodes <- unique(c(
            mask_Model_Reach_by_Typology$KP_start,
            mask_Model_Reach_by_Typology$KP_end
        ))
        Model_Reach_nodes <- sort(Model_Reach_nodes, decreasing = Logical_decreasing)

        Key_Info_Typology_Model_Reach[[i]]$Model_Reach_nodes <- Model_Reach_nodes

        # Grid from interpolation passing by Model_Reach_nodes
        # KP_Model_Reach_nodes is defined by XR
        KP_Model_Reach_nodes <- interpolation_specific_points(
            total_points = total_points_discretization,
            specific_nodes = Model_Reach_nodes
        )

        KP_Model_Reach_nodes <- sort(KP_Model_Reach_nodes, decreasing = Logical_decreasing)

        # Assign reach at the KP_Model_Reach_nodes
        KP_reach_Model_Reach_nodes <- assign_reach_from_a_grid(
            reach_KP_boundaries = mask_Model_Reach_by_Typology,
            grid = KP_Model_Reach_nodes,
            Logical_decreasing = Logical_decreasing
        )

        # These information allow build RUGFile
        if (Logical_decreasing) {
            RUGFile <- data.frame(
                id_reach = KP_reach_Model_Reach_nodes$reaches[-nrow(KP_reach_Model_Reach_nodes)],
                KP_start = KP_reach_Model_Reach_nodes$grid[-1],
                KP_end = KP_reach_Model_Reach_nodes$grid[-nrow(KP_reach_Model_Reach_nodes)],
                Kmin = 0,
                Kflood = 0
            )
        } else {
            RUGFile <- data.frame(
                id_reach = KP_reach_Model_Reach_nodes$reaches[-nrow(KP_reach_Model_Reach_nodes)],
                KP_start = KP_reach_Model_Reach_nodes$grid[-nrow(KP_reach_Model_Reach_nodes)],
                KP_end = KP_reach_Model_Reach_nodes$grid[-1],
                Kmin = 0,
                Kflood = 0
            )
        }
        Key_Info_Typology_Model_Reach[[i]]$RUGFile <- RUGFile
        # Now, the middle value of each interval is taken to get the KP_grid, which will be used for all the calculations

        # Calculate the middle value for each row
        Key_Info_Typology_Model_Reach[[i]]$KP_grid <- (RUGFile$KP_start + RUGFile$KP_end) / 2
        Key_Info_Typology_Model_Reach[[i]]$reach <- RUGFile$id_reach
        Key_Info_Typology_Model_Reach[[i]]$Logical_decreasing <- Logical_decreasing
    }
    return(Key_Info_Typology_Model_Reach)
}



constructor_RUGFile_covariate_grid <- function(
    Key_Info_Typology_Model_Reach) {
    ########################
    # Get discretization of the covariate, same in Kmin or Kflood
    ########################
    covariate_grid_lists <- lapply(Key_Info_Typology_Model_Reach, function(SR) {
        SR$KP_grid
    })
    covariate_grid <- data.frame(covariate = unlist(covariate_grid_lists))
    rownames(covariate_grid) <- NULL

    ########################
    # RUGFile
    ########################
    # Extract all RUGFile data frames from the main list
    all_RUGFiles <- lapply(Key_Info_Typology_Model_Reach, function(channel) {
        channel$RUGFile
    })

    # Row-bind all RUGFile data frames into a single data frame
    RUGFile_Mage_Final <- do.call(rbind, all_RUGFiles) %>% arrange(id_reach)
    # Remove row names
    rownames(RUGFile_Mage_Final) <- NULL
    #######################################
    # End RUGFile
    #######################################

    return(
        list(
            RUGFile_Mage_Final = RUGFile_Mage_Final,
            covariate_grid = covariate_grid
        )
    )
}


constructor_spatialization_matrix <- function(
    K_SR) {
    # Get Z file for all XR in Kmin
    Z_all_K_SR <- lapply(K_SR, function(channel) {
        lapply(channel, function(SR) {
            SR$Z
        })
    })

    # Get the reaches of the Z file for all XR in Kmin
    reaches_Z_all_K_SR <- lapply(K_SR, function(channel) {
        lapply(channel, function(SR) {
            SR$reach
        })
    })

    # Flatten the nested list of matrices
    flat_K_SR <- unlist(Z_all_K_SR, recursive = FALSE)
    reaches_flat_K_SR <- cbind(unlist(unlist(reaches_Z_all_K_SR, recursive = FALSE)))
    rownames(reaches_flat_K_SR) <- NULL

    # Pass the flattened list to block_diagonal_matrix
    Z_MatrixK_not_ordered <- do.call(block_diagonal_matrix, flat_K_SR)

    Z_MatrixK <- as.matrix(
        data.frame(
            cbind(
                reaches_flat_K_SR,
                Z_MatrixK_not_ordered
            )
        ) %>%
            arrange(X1) %>%
            select(-X1)
    )

    return(Z_MatrixK)
}

Estimation_Mage <- function(
    Key_Info_Typology_Model_Reach,
    Input_Typology,
    path_experiment,
    file_main_path,
    all_cal_case,
    do_calibration) {
    mod_polynomials <- list()
    counter_model <- 1

    for (id_cal_case in 1:length(all_cal_case)) {
        # Source the input data for the experiments
        path_Experiment_Input_Data <- file.path(file_main_path, "Experiments_Input_Data", all_cal_case[id_cal_case])
        if (!file.exists(path_Experiment_Input_Data)) {
            stop(paste0(
                "The experiment input data named : '",
                all_cal_case[id_cal_case],
                "' does not exist"
            ))
        }
        source(path_Experiment_Input_Data)
        path_polynomial <- file.path(
            path_experiment,
            sub("\\.r$", "", all_cal_case[id_cal_case])
        )
        # Path to save results
        workspace_user <- file.path(path_polynomial, "BaM")
        path_post_traitement <- file.path(workspace_user, "post_traitement")
        path_post_traitement_data <- file.path(path_post_traitement, "RData")
        if (!dir.exists(path_polynomial)) {
            dir.create(path_polynomial)
        } else {
            if (do_calibration) {
                file.remove(list.files(path_polynomial, full.names = TRUE, recursive = TRUE))
            }
        }
        MAGE_polynomial_subfolder <- file.path(path_polynomial, "model_mage")

        if (!dir.exists(MAGE_polynomial_subfolder)) {
            dir.create(MAGE_polynomial_subfolder)
        }
        # Copy Main MAGE folder to local folder for calibration
        copy_folder(MAGE_main_folder, MAGE_polynomial_subfolder)

        if (!dir.exists(workspace_user)) {
            dir.create(workspace_user)
        } else {
            if (do_calibration) {
                file.remove(list.files(workspace_user, full.names = TRUE, recursive = TRUE))
            }
        }
        if (!dir.exists(path_post_traitement)) {
            dir.create(path_post_traitement)
        }

        if (!dir.exists(path_post_traitement_data)) {
            dir.create(path_post_traitement_data)
        }

        constructor_results <- constructor_RUGFile_covariate_grid(Key_Info_Typology_Model_Reach = Key_Info_Typology_Model_Reach)

        RUGFile_Mage_Final <- constructor_results$RUGFile_Mage_Final
        covariate_grid <- constructor_results$covariate_grid
        ############################################
        # Kmin environment (encapsulated in Kmin_SR)
        ############################################

        # Check if size is respected between ID and Kmin
        if (length(Input_Typology) != length(Input_Kmin_Key_SR_MR)) stop("Size must be equal between Input_Kmin_Key_SR_MR and Input_Typology")

        # Assign properties of each SR in XR structure
        Kmin_SR <- SR_constructor(
            SR_key_HM = Input_Kmin_Key_SR_MR,
            Key_Info_Typology_Model_Reach = Key_Info_Typology_Model_Reach
        )

        Z_MatrixKmin <- constructor_spatialization_matrix(K_SR = Kmin_SR)

        # Check size between RUGFile_Mage_Final and Z_file
        if (nrow(RUGFile_Mage_Final) != nrow(Z_MatrixKmin)) stop("RUGFile_Mage_Final must have the same size as Z file (spatialisation)")

        Kmin_prior <- extract_priors(Kmin_SR)
        ############################################
        # End Kmin environment
        ############################################

        ############################################
        # Kflood environment (encapsulated in Kflood_SR)
        ############################################

        # Check if size is respected between ID and Kflood
        if (length(Input_Typology) != length(Input_Kflood_Key_SR_MR)) stop("Size must be equal between Input_Kflood_Key_SR_MR and Input_Typology")

        # Assign properties of each SR in XR structure
        Kflood_SR <- SR_constructor(
            SR_key_HM = Input_Kflood_Key_SR_MR,
            Key_Info_Typology_Model_Reach = Key_Info_Typology_Model_Reach
        )

        Z_MatrixKflood <- constructor_spatialization_matrix(K_SR = Kflood_SR)

        # Check size between RUGFile_Mage_Final and Z_file
        if (nrow(RUGFile_Mage_Final) != nrow(Z_MatrixKflood)) stop("RUGFile_Mage_Final must have the same size as Z file (spatialisation)")

        Kflood_prior <- extract_priors(Kflood_SR)

        ############################################
        # End Kflood environment
        ############################################

        if (!(nrow(covariate_grid) == nrow(Z_MatrixKmin) && nrow(covariate_grid) == nrow(Z_MatrixKflood))) {
            stop("Number of rows of covariate_grid must be equal to both Z_MatrixKmin and Z_MatrixKflood")
        }

        # Write covariate discretization, available for Kmin and Kflood:
        write.table(
            covariate_grid,
            file = file.path(workspace_user, "grid_covariate_non_normalized.txt"), row.names = F
        )
        mageDir <- c(paste0(path_polynomial, "/model_mage/", all_events, "/"))

        RUGFile_paths <- paste0(mageDir, mage_projet_name, ".RUG")

        # j is a index for multi-events
        for (id_multi_event in seq_along(RUGFile_paths)) {
            write_RUGFile(
                RUG_path = RUGFile_paths[id_multi_event],
                RUGFile_data = RUGFile_Mage_Final,
                RUG_format = "%1s%3d      %10.3f%10.3f%10.2f%10.2f"
            )
        }

        zFileKmin <- file.path(workspace_user, "Zfile_Kmin.txt")
        zFileKmoy <- file.path(workspace_user, "Zfile_Kflood.txt")

        xtra <- xtraModelInfo(
            fname = "Config_setup.txt",
            object = list(
                exeFile = paste0(MAGE_executable, " ", command_line_MAGE),
                # commandLine = command_line_MAGE,
                version = "8",
                mageDir = mageDir,
                repFile = paste0(mage_projet_name, ".REP"),
                zKmin = Z_MatrixKmin,
                zFileKmin = zFileKmin,
                doExpKmin = FALSE,
                zKmoy = Z_MatrixKflood,
                zFileKmoy = zFileKmoy,
                doExpKmoy = FALSE
            )
        )

        theta_param <- c(Kmin_prior, Kflood_prior)

        mod_polynomials[[counter_model]] <- model(
            ID = "MAGE_ZQV",
            nX = 4,
            nY = 5,
            par = theta_param,
            xtra = xtra
        )
        mod <- mod_polynomials[[counter_model]]

        prior_theta_param <- get_init_prior(theta_param)

        data <- dataset(X = X, Y = Y, Yu = Yu, data.dir = file.path(workspace_user))

        jump_MCMC_theta_param <- ifelse(prior_theta_param != 0,
            prior_theta_param[(prior_theta_param != 0)] * 0.1,
            jump_MCMC_theta_param_user
        )

        jump_MCMC_error_model <- list()
        for (ind in 1:length(prior_error_model)) {
            jump_MCMC_error_model[[ind]] <- ifelse(prior_error_model[[ind]] > threshold_jump_MCMC_error_model,
                prior_error_model[[ind]][(prior_error_model[[ind]] != 0)] * 0.1,
                jump_MCMC_error_model_user
            )
        }

        mcmcOptions_user <- mcmcOptions(
            nCycles = 100,
            nAdapt = 100,
            manualMode = TRUE,
            thetaStd = jump_MCMC_theta_param,
            gammaStd = jump_MCMC_error_model
        )
        mcmcCooking_user <- mcmcCooking(
            burn = 0.5,
            nSlim = 10
        )
        mcmcSummary_user <- mcmcSummary(xtendedMCMC.fname = "Results_xtendedMCMC.txt")

        # Save all data used during calibration for prediction
        save(mod, data, remant_error_list,
            mcmcOptions_user, mcmcCooking_user,
            mcmcSummary_user, workspace_user,
            file = file.path(path_post_traitement_data, "BaM_objects.RData")
        )
        BaM(
            mod = mod,
            data = data,
            remnant = remant_error_list,
            mcmc = mcmcOptions_user,
            cook = mcmcCooking_user,
            summary = mcmcSummary_user,
            residuals = residualOptions(),
            pred = NULL,
            doCalib = TRUE,
            doPred = FALSE,
            na.value = -9999,
            run = FALSE,
            preClean = FALSE,
            workspace = workspace_user,
            # dir.exe = file.path(find.package("RBaM"), "bin"),
            # name.exe = "BaM",
            predMaster_fname = "Config_Pred_Master.txt"
        )
        counter_model <- counter_model + 1

        dir_cf <- file.path(workspace_user, "Config_BaM.txt")

        if (do_calibration) {
            system2(
                command = file.path(dir_exe_BaM, "BaM"),
                args = c("-cf", file.path(workspace_user, "Config_BaM.txt")),
                wait = FALSE
            )
        }
    }
    # return(list(
    #     Z_MatrixKmin = Z_MatrixKmin,
    #     Z_MatrixKflood = Z_MatrixKflood,
    #     mcmc = mcmc,
    #     Kmin_prior = Kmin_prior,
    #     Kflood_prior = Kflood_prior
    # ))
}
