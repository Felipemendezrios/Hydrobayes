last <- function(data) {
    utils::tail(data, n = 1)
}

get_ST_fortran <- function(
    path_ST, skip = 1) {
    bathy <- readLines(file.path(path_ST))
    bathy <- bathy[-c(1:skip)]

    # Remove headers and get xyz values

    fields <- strsplit(trimws(bathy), "\\s+")
    n_fields <- lengths(fields)

    is_header <- n_fields >= 5

    is_end <- sapply(fields, function(x) {
        length(x) >= 3 && all(x[1:3] == "999.9990")
    })

    is_xyz <- (n_fields %in% c(3, 4)) & !is_end

    block_id <- cumsum(is_header) # each header starts a new block
    block_id[is_end] <- NA # remove end marker lines

    # keep only headers
    header_lines <- fields[is_header]

    # get first and before-last values
    header_info <- t(sapply(header_lines, function(x) {
        first_val <- x[1]
        before_last <- x[length(x) - 1]
        c(profile_id = first_val, KP = before_last)
    }))

    header_info <- as.data.frame(header_info, stringsAsFactors = FALSE)

    xyz_list <- lapply(which(is_xyz), function(i) {
        nums <- as.numeric(fields[[i]][1:3])
        lbl <- if (length(fields[[i]]) >= 4) fields[[i]][4] else NA
        blk <- block_id[i] # which header/block this line belongs to
        c(nums, label = lbl, block = blk)
    })

    xyz_df <- as.data.frame(do.call(rbind, xyz_list), stringsAsFactors = FALSE)
    colnames(xyz_df) <- c("X", "Y", "Z", "label", "block")
    xyz_df$profile_id <- header_info$profile_id[as.numeric(xyz_df$block)]
    xyz_df$KP <- header_info$KP[as.numeric(xyz_df$block)]

    # optional: remove temporary block column
    xyz_df$block <- NULL

    # Convert numeric columns
    xyz_df$X <- as.numeric(xyz_df$X)
    xyz_df$Y <- as.numeric(xyz_df$Y)
    xyz_df$Z <- as.numeric(xyz_df$Z)
    xyz_df$profile_id <- as.numeric(xyz_df$profile_id)
    xyz_df$KP <- as.numeric(xyz_df$KP)
    return(list(
        xyz_df = xyz_df,
        header_info = header_info
    ))
}


Cross_sections_interpolation_and_export_rectangular_channel <- function(
    So,
    Pk,
    positions_banks,
    borders_heigh,
    path_export) {
    PK <- sort(Pk)
    ST <- data.frame(PK = PK)
    ST$z <- NA

    j <- 1
    for (i in PK) {
        ST$z[j] <- So * (PK[length(PK)] - i)
        j <- j + 1
    }

    con <- file(path_export, open = "w")

    format_str_entete <- "%6d%6d%6d%6d%13.8f  %-20s"
    format_str_xyz <- "%10.8f %10.8f %10.8f %-3s"

    # Use a for loop to write lines
    for (i in 1:nrow(ST)) {
        # Appliquer ligne par ligne
        entete_line <- sprintf(
            format_str_entete,
            as.integer(i),
            as.integer(0),
            as.integer(0),
            as.integer(6),
            as.numeric(ST$PK[i]),
            as.character("")
        )
        # Écrire dans un fichier
        writeLines(entete_line, con)

        XYZ_df <- data.frame(
            x = c(rep(ST$PK[i], 6), 999.999),
            y = c(rep(positions_banks[1], 3), rep(positions_banks[2], 3), 999.999),
            z = c(
                rep(borders_heigh + ST$z[i], 2),
                rep(ST$z[i], 2),
                rep(borders_heigh + ST$z[i], 2),
                999.999
            ),
            char = c("", "rg", rep("", 2), "rd", "", "")
        )

        # Générer toutes les lignes formatées
        XYZ <- apply(XYZ_df, 1, function(row) {
            sprintf(
                format_str_xyz,
                as.numeric(row["x"]),
                as.numeric(row["y"]),
                as.numeric(row["z"]),
                as.character(row["char"])
            )
        })
        # Écrire dans un fichier
        writeLines(XYZ, con)
    }

    # Close the connection
    close(con)
}


# Assuming Kflood_SU is your nested list
extract_priors <- function(nested_list) {
    priors <- list()

    # Recursive function to extract 'prior' elements
    extract <- function(x) {
        if (is.list(x)) {
            if ("prior" %in% names(x)) {
                priors <<- c(priors, list(x$prior$Config_Model))
            } else {
                lapply(x, extract)
            }
        }
    }

    # Apply the recursive function
    extract(nested_list)

    # Flatten the list of priors
    flat_priors <- unlist(priors, recursive = FALSE)

    return(flat_priors)
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

get_all_init_prior_theta <- function(parameter) {
    init_priors <- numeric(0)
    for (i in parameter) {
        param <- i
        init_priors <- c(init_priors, param$init)
    }
    return(init_priors)
}

get_prior_info_plot_SU <- function(K_SU) {
    data_prior_info_to_plot <-
        lapply(
            names(K_SU),
            function(typology) {
                su_list <- lapply(
                    names(K_SU[[typology]]),
                    function(su) {
                        plot_info <- K_SU[[typology]][[su]]$prior$plot

                        prior_summary <- plot_info$prior_summary
                        spatial <- plot_info$spatial_positions_prior
                        correlation <- plot_info$prior_correlation_SU
                        prior_spat <- plot_info$prior_spatialization_matrix
                        variances <- plot_info$variances

                        all_info <- data.frame(
                            Typology = prior_summary$Typology,
                            SU = prior_summary$SU,
                            param_name = prior_summary$param_name,
                            prior.dist = prior_summary$prior.dist,
                            par_1 = prior_summary$par_1,
                            par_2 = prior_summary$par_2,
                            covariate = spatial$covariate,
                            scaled = spatial$scaled
                        )

                        list(
                            prior_summary = prior_summary,
                            spatial_positions_prior = spatial,
                            prior_correlation_SU = correlation,
                            prior_spat = prior_spat,
                            variances = variances,
                            all_info = all_info
                        )
                    }
                )

                names(su_list) <- names(K_SU[[typology]])

                su_list
            }
        )
    names(data_prior_info_to_plot) <- names(K_SU)

    return(data_prior_info_to_plot)
}

traitement_prior <- function(data_prior_info_K) {
    prior_envelope_all <- c()
    for (Typology in names(data_prior_info_K)) {
        for (SU in seq_along(data_prior_info_K[[Typology]])) {
            data_prior_local <- data_prior_info_K[[Typology]][[SU]]$all_info

            data_matrix_spatial_prior_local <- data_prior_info_K[[Typology]][[SU]]$prior_spat

            data_variances_prior_local <- data_prior_info_K[[Typology]][[SU]]$variances

            if (nrow(data_prior_local) == 1 & all(data_prior_local$prior.dist == "FIX")) {
                prior_envelope_all <- rbind(
                    prior_envelope_all,
                    data.frame(
                        x = NA,
                        ymin = NA,
                        ymax = NA,
                        ID = "Prior",
                        typology = Typology,
                        id_reach_SU = SU
                    )
                )
            } else {
                # Simulate T replicates
                simT <- t(mvtnorm::rmvnorm(
                    1000,
                    data_prior_local$par_1,
                    data_variances_prior_local
                ))
                # Transform into K replicates
                simK <- data_matrix_spatial_prior_local %*% simT
                # Prior realization in K space
                prior_realization <- data.frame(
                    x = data_prior_local$scaled,
                    Value = simK,
                    ID = "Prior"
                ) %>% tidyr::pivot_longer(
                    cols = -c(x, ID),
                    names_to = "Iteration",
                    values_to = "Value"
                )

                # Quantify uncertainty at 95%
                prior_envelope <- prior_realization %>%
                    group_by(x) %>%
                    summarise(
                        ymin = quantile(Value, probs = 0.025, na.rm = TRUE),
                        ymax = quantile(Value, probs = 0.975, na.rm = TRUE),
                        ID = "Prior", # so we can map to fill
                        .groups = "drop"
                    ) %>%
                    mutate(
                        typology = Typology,
                        id_reach_SU = SU
                    )
                # If n=0, add a point at the begging for plotting
                if (nrow(prior_envelope) == 1) {
                    prior_envelope <-
                        rbind(
                            prior_envelope %>% mutate(x = -1),
                            prior_envelope
                        )
                }
                prior_envelope_all <- rbind(prior_envelope_all, prior_envelope)
            }
        }
    }
    return(prior_envelope_all)
}

write_RUGFile <- function(RUG_path,
                          RUGFile_data,
                          RUG_format) {
    # Open a .RUG file for writing
    fileConn <- file(RUG_path, "w")
    # Write the first line as a comment
    writeLines("* This file is generated by PAMHYR, please don't modify", fileConn)
    formatted_lines <- sapply(1:nrow(RUGFile_data), function(i) {
        sprintf(
            RUG_format,
            "K",
            RUGFile_data$id_reach[i],
            RUGFile_data$KP_start[i],
            RUGFile_data$KP_end[i],
            RUGFile_data$Kmin[i],
            RUGFile_data$Kflood[i]
        )
    })

    # Write the formatted lines to the file
    writeLines(formatted_lines, fileConn)


    # Close the file
    close(fileConn)
}

copy_folder <- function(source, destination) {
    # Create the destination directory if it doesn't exist
    if (!dir.exists(destination)) {
        dir.create(destination, recursive = TRUE)
    } else {
        # Clean the destination folder: delete all files and subdirectories
        files_to_delete <- list.files(destination, full.names = TRUE, recursive = TRUE)
        if (length(files_to_delete) > 0) {
            file.remove(files_to_delete)
        }
        # Remove empty subdirectories (if any)
        dirs_to_delete <- list.dirs(destination, full.names = TRUE, recursive = FALSE)
        if (length(dirs_to_delete) > 0) {
            unlink(dirs_to_delete, recursive = TRUE)
        }
    }

    # List all files and subdirectories in the source
    items <- list.files(source, full.names = TRUE, recursive = FALSE)

    for (item in items) {
        item_relative <- basename(item)
        dest_item <- file.path(destination, item_relative)

        if (dir.exists(item)) {
            # If the item is a subdirectory, recursively copy it
            copy_folder(item, dest_item)
        } else {
            # If the item is a file, copy it
            file.copy(item, dest_item)
        }
    }
}

convert_9999_to_NA <- function(values) {
    values[values == -9999] <- NA
    return(values)
}

read_fortran_data <- function(file_path, col_widths_RUGFile, skip = 0) {
    # Read the file with the fixed-width format
    data <- utils::read.fwf(file_path, widths = col_widths_RUGFile, header = FALSE, skip = skip)
    data <- data[, -3]
    colnames(data) <-
        c(
            "",
            "id_reach",
            "KP_start",
            "KP_end",
            "Kmin",
            "Kflood"
        )
    return(data)
}

get_prior_distribution <- function(K_prior) {
    sapply(K_prior, function(x) x$prior$dist)
}

get_param_vector_MAP_values <- function(SU_Kmin, SU_Kflood, MAP) {
    param_MAP_values <- vector("list", sum(lengths(SU_Kmin)) + sum(lengths(SU_Kflood)))
    names_all <- c(
        unlist(lapply(names(SU_Kmin), function(l1) {
            paste0("Kmin_", l1, "_", names(SU_Kmin[[l1]]))
        })),
        unlist(lapply(names(SU_Kflood), function(l1) {
            paste0("Kflood_", l1, "_", names(SU_Kflood[[l1]]))
        }))
    )
    names(param_MAP_values) <- names_all

    counter <- 1
    MAP_idx <- 1
    SU_Kmin_Kflood <- c(SU_Kmin, SU_Kflood)

    for (id_reach in SU_Kmin_Kflood) { # Browse the id reach (tributary, main channel, etc) in the Kmin and Kflood
        for (id_SU in id_reach) { # Browse the spatial unit (SU) of the id_reach
            # id_SU <- SU_Kmin[[3]][[1]]
            # idea: keep the order of the config model !

            # Get initial values of all parameters
            param_temp <- get_init_prior(extract_priors(id_SU), FIX_dist = TRUE)
            # If prior distribution is FIX, then keep the initial value.
            # Otherwise, modify by the MAP value to calculate residuals
            idx <- which(get_prior_distribution(id_SU$prior) != "FIX")

            # Get the number of distribution different to FIX distribution to move indicator in the MAP variable
            n <- length(idx)

            if (n != 0) {
                param_temp[idx] <- MAP[MAP_idx:(MAP_idx + n - 1)]
            }

            MAP_idx <- MAP_idx + n

            param_MAP_values[[counter]] <- param_temp

            counter <- counter + 1
        }
    }
    return(unlist(param_MAP_values))
}

remnantErrorModel_default <- function(name) {
    RBaM::remnantErrorModel(
        fname = name,
        funk = "Constant",
        par = list(parameter(
            name = "intercept",
            init = 0.01,
            prior.dist = "Exponential",
            prior.par = c(0, 0.01) # 0 is the threshold and 0.01 is the scale
        ))
    )
}


# Set paths

load_experiment <- function(file_main_path, cal_case, path_experiment, all_events) {
    path_input <- file.path(file_main_path, "Experiments_Input_Data", cal_case)

    if (!file.exists(path_input)) {
        stop("Experiment input file does not exist: ", cal_case)
    }

    source(path_input)

    path_polynomial <- file.path(
        path_experiment,
        sub("\\.r$", "", cal_case)
    )

    path_BaM_folder <- file.path(path_polynomial, "BaM")

    path_plot_folder <- file.path(path_BaM_folder, "post_traitement")

    path_RData <- file.path(path_plot_folder, "RData")

    dir.create(path_plot_folder, showWarnings = FALSE)
    dir.create(path_RData, showWarnings = FALSE)

    path_model_HM <- file.path(
        path_polynomial,
        "model_mage"
    )

    path_model_HM_events <- paste0(
        file.path(
            path_model_HM,
            all_events
        ),
        "/"
    )


    return(list(
        path_experiment = path_input,
        path_polynomial = path_polynomial,
        path_BaM_folder = path_BaM_folder,
        path_plot_folder = path_plot_folder,
        path_RData = path_RData,
        path_model_HM_events = path_model_HM_events,
        path_model_HM = path_model_HM
    ))
}

RUGFile_post_estimation <- function(
    RUGFile_structure,
    Z_MatrixKmin,
    Z_MatrixKflood,
    Kmin_SU,
    Kflood_SU,
    MAP_Kmin_Kflood,
    dist_prior_Kmin,
    dist_prior_Kflood) {
    # Check RUGFile structure
    if (any(colnames(RUGFile_structure) != c("id_reach", "KP_start", "KP_end", "Kmin", "Kflood"))) stop('RUGFILE_structure must have this information: "id_reach" "KP_start" "KP_end"   "Kmin"     "Kflood"  ')

    # Order RUGFile_structure following id_reach
    RUGFile_structure_ordered <- RUGFile_structure[order(RUGFile_structure$id_reach), ]

    # Check MAP is a vector
    if (!is.vector(MAP_Kmin_Kflood)) stop("MAP_Kmin_Kflood must be vectors")

    # Kmin
    # Get index of the estimated parameters
    id_param_estimated_Kmin <- which(dist_prior_Kmin != "FIX")
    Fix_dist_positions_Kmin <- which(dist_prior_Kmin == "FIX")

    # Kflood
    # Get index of the estimated parameters
    id_param_estimated_Kflood <- which(dist_prior_Kflood != "FIX")
    Fix_dist_positions_Kflood <- which(dist_prior_Kflood == "FIX")


    # Check size of MAP and Z_Matrix
    if ((length(id_param_estimated_Kmin) + length(id_param_estimated_Kflood)) != length(MAP_Kmin_Kflood)) stop(" Summing number of columns between id_param_estimated_Kmin and id_param_estimated_Kflood, the result must have the same number of columns than length of MAP_Kmin_Kflood")

    param_MAP_values <- get_param_vector_MAP_values(
        SU_Kmin = Kmin_SU,
        SU_Kflood = Kflood_SU,
        MAP = MAP_Kmin_Kflood
    )

    RUGFile_structure_ordered$Kmin <- Z_MatrixKmin %*% param_MAP_values[1:ncol(Z_MatrixKmin)]
    RUGFile_structure_ordered$Kflood <- Z_MatrixKflood %*% param_MAP_values[(ncol(Z_MatrixKmin) + 1):(ncol(Z_MatrixKmin) + ncol(Z_MatrixKflood))]

    return(RUGFile_structure_ordered)
}

prior_distributions <- function(distribution,
                                param1,
                                param2) {
    if (distribution != "Gaussian" & distribution != "LogNormal" &
        distribution != "Uniform" & distribution != "FIX") {
        stop("Prior distribution is not supported. \nPlease ensure that distribution is Gaussian, \n LogNormal, Uniform or FIX")
    }

    if (distribution == "Gaussian") {
        prior_realization <- stats::rnorm(1000,
            mean = as.numeric(param1),
            sd = as.numeric(param2)
        )
    } else if (distribution == "LogNormal") {
        prior_realization <- stats::rlnorm(1000,
            meanlog = as.numeric(param1),
            sdlog = as.numeric(param2)
        )
    } else if (distribution == "Uniform") {
        if (param1 <= param2) stop("To use Uniform distribution, \nfirst value introduce must be lower than second value in the vector")
        prior_realization <- stats::runif(1000,
            min = as.numeric(param1),
            max = as.numeric(param2)
        )
    }
    return(prior_realization)
}

# Generate prior realizations
get_prior_density <- function(prior_list) {
    prior_list <- Filter(function(x) x$prior$dist != "FIX", prior_list)

    priors_realization <- lapply(prior_list, function(x) {
        prior_distributions(
            distribution = x$prior$dist,
            param1 = x$prior$par[1],
            param2 = x$prior$par[2]
        )
    })

    names(priors_realization) <- vapply(prior_list, `[[`, "", "name")

    as.data.frame(priors_realization)
}

# Combine prior and posterior
combine_prior_posterior_MAP <- function(prior_density, mcmc, MAP) {
    if (length(prior_density) == 0) {
        return(NULL)
    }
    mcmc_extracted <- mcmc[, names(prior_density), drop = FALSE]

    if (any(length(MAP) != ncol(mcmc_extracted) | ncol(mcmc_extracted) != length(prior_density))) stop("Inconsistency of the number of estimated parameters")

    names(MAP) <- colnames(mcmc_extracted)

    if (length(mcmc_extracted) != 0) {
        do.call(rbind, lapply(names(prior_density), function(nm) {
            rbind(
                data.frame(
                    value = prior_density[[nm]],
                    Distributions = "Prior",
                    id = nm
                ),
                data.frame(
                    value = mcmc_extracted[[nm]],
                    Distributions = "Posterior",
                    id = nm
                ),
                data.frame(
                    value = MAP[[nm]],
                    Distributions = "MAP",
                    id = nm
                )
            )
        }))
    } else {
        NULL
    }
}
extract_param_values <- function(Input_Key_SU_MR) {
    param_values_df <- data.frame()

    for (typology in names(Input_Key_SU_MR)) {
        for (SU in names(Input_Key_SU_MR[[typology]])) {
            if (Input_Key_SU_MR[[typology]][[SU]]$prior$config$distribution == "FIX") {
                param_values_df <- rbind(
                    param_values_df,
                    data.frame(
                        typology = typology,
                        SU = SU,
                        mu = NA,
                        sigma = NA
                    )
                )
            } else {
                param_values <- Input_Key_SU_MR[[typology]][[SU]]$prior$config$param_values

                # Valeurs NULL -> NA
                mu <- if (!is.null(param_values$mu)) param_values$mu else NA_real_
                sigma <- if (!is.null(param_values$sigma)) param_values$sigma else NA_real_

                param_values_df <- rbind(
                    param_values_df,
                    data.frame(
                        typology = typology,
                        SU = SU,
                        mu = mu,
                        sigma = sigma
                    )
                )
            }
        }
    }
    return(param_values_df)
}
