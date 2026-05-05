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
                priors <<- c(priors, list(x$prior))
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
    MAP_Kmin_Kflood,
    init_guess_prior_Kmin,
    dist_prior_Kmin,
    init_guess_prior_Kflood,
    dist_prior_Kflood) {
    # Check RUGFile structure
    if (any(colnames(RUGFile_structure) != c("id_reach", "KP_start", "KP_end", "Kmin", "Kflood"))) stop('RUGFILE_structure must have this information: "id_reach" "KP_start" "KP_end"   "Kmin"     "Kflood"  ')
    # Check MAP is a vector
    if (!is.vector(MAP_Kmin_Kflood)) stop("MAP_Kmin_Kflood must be vectors")

    # Check size of MAP and Z_Matrix
    ncol_Z_Kmin <- switch(as.character(ncol(Z_MatrixKmin)),
        "1" = 0,
        ncol(Z_MatrixKmin)
    )

    ncol_Z_Kflood <- switch(as.character(ncol(Z_MatrixKflood)),
        "1" = 0,
        ncol(Z_MatrixKflood)
    )

    if ((ncol_Z_Kmin + ncol_Z_Kflood) != length(MAP_Kmin_Kflood)) stop(" Summing number of columns between Z_MatrixKmin and Z_MatrixKflood, the result must have the same number of columns than length of MAP_Kmin_Kflood")

    # Kmin

    # Get index of the estimated parameters
    id_param_estimated_Kmin <- which(dist_prior_Kmin != "FIX")
    Fix_dist_positions_Kmin <- which(dist_prior_Kmin == "FIX")

    # Case 1: all parameters are fixed
    if (length(id_param_estimated_Kmin) == 0) {
        MAP_all_param_Kmin <- init_guess_prior_Kmin
    } else if (
        # Case 2: No one is fixed
        length(Fix_dist_positions_Kmin) == 0) {
        MAP_all_param_Kmin <- MAP[1:ncol(Z_MatrixKmin)]
    } else { # Case 3: mixted
        # Extract only MCMC of the parameters theta (estimated)
        MAP_estimated <- cbind(value = MAP[1:ncol(Z_MatrixKmin)], id_param = id_param_estimated_Kmin)
        # Get values constant value of the fix distributions
        vals <- init_guess_prior_Kmin[Fix_dist_positions_Kmin]

        for (i in seq_along(Fix_dist_positions_Kmin)) {
            MAP_estimated <- rbind(
                MAP_estimated,
                data.frame(
                    value = vals[i],
                    id_param = Fix_dist_positions_Kmin[i]
                )
            )
        }
        MAP_estimated$row_id <- ave(MAP_estimated$id_param,
            MAP_estimated$id_param,
            FUN = seq_along
        )
        MAP_all_param_Kmin <- as.vector(unlist(MAP_estimated %>%
            pivot_wider(
                id_cols = row_id,
                names_from = id_param,
                values_from = value
            ) %>%
            select(-row_id)))
    }

    # Kflood
    # Get index of the estimated parameters
    id_param_estimated_Kflood <- which(dist_prior_Kflood != "FIX")
    Fix_dist_positions_Kflood <- which(dist_prior_Kflood == "FIX")

    # Case 1: all parameters are fixed
    if (length(id_param_estimated_Kflood) == 0) {
        MAP_all_param_Kflood <- init_guess_prior_Kflood
    } else if (
        # Case 2: No one is fixed
        length(Fix_dist_positions_Kflood) == 0) {
        MAP_all_param_Kflood <- MAP[1:ncol(Z_MatrixKflood)]
    } else { # Case 3: mixted
        # Extract only MCMC of the parameters theta (estimated)
        MAP_estimated <- cbind(value = MAP[1:ncol(Z_MatrixKflood)], id_param = id_param_estimated_Kflood)
        # Get values constant value of the fix distributions
        vals <- init_guess_prior_Kflood[Fix_dist_positions_Kflood]

        for (i in seq_along(Fix_dist_positions_Kflood)) {
            MAP_estimated <- rbind(
                MAP_estimated,
                data.frame(
                    value = vals[i],
                    id_param = Fix_dist_positions_Kflood[i]
                )
            )
        }
        MAP_estimated$row_id <- ave(MAP_estimated$id_param,
            MAP_estimated$id_param,
            FUN = seq_along
        )
        MAP_all_param_Kflood <- as.vector(unlist(MAP_estimated %>%
            pivot_wider(
                id_cols = row_id,
                names_from = id_param,
                values_from = value
            ) %>%
            select(-row_id)))
    }

    RUGFile_structure$Kmin <- Z_MatrixKmin %*% MAP_all_param_Kmin
    RUGFile_structure$Kflood <- Z_MatrixKflood %*% MAP_all_param_Kflood
    return(RUGFile_structure)
}
