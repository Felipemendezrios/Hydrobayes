RUGFile_constructor_SU_correction <- function(Key_Info_XR_MR) {
    # Calculate the length of RUGFile using Key_Info_XR_MR
    length_RUGFile <- sum(
        sapply(Key_Info_XR_MR, function(id) {
            length(id$reach)
        })
    )
    # Initialize RUGFile_first_part as a data frame with 5 columns
    RUGFile_first_part <- data.frame(
        matrix(
            data = 0,
            nrow = length_RUGFile,
            ncol = 5 # RUGFILE columns : reach, kp_start, kp_end, kmin,Kflood
        )
    )
    colnames(RUGFile_first_part) <- c("reach", "KP_start", "KP_end", "Kmin", "Kflood")
    i1 <- 1
    for (id_Key in Key_Info_XR_MR) {
        i2 <- i1 + length(id_Key$KP_grid)
        RUGFile_first_part[i1:i2, c(1, 2, 3)] <- data.frame(
            id_reach = id_Key$reach[-length(id_Key$reach)],
            KP_start = id_Key$KP_grid[-length(id_Key$KP_grid)],
            KP_end = id_Key$KP_grid[-1]
        )
        i1 <- i2 + 1
    }
    return(RUGFile = RUGFile_first_part)
}


SU_constructor <- function(
    SU_key_HM,
    Key_Info_Typology_Model_Reach) {
    # Initialize structure by fixing layer 1 : XR
    # Layer 2 : Kmin or Kflood is handle outside this function
    SU <- vector(mode = "list", length = length(SU_key_HM))
    names(SU) <- names(SU_key_HM)

    # Extract only KP boundaries points keeping the same structure of XR
    SU_KP_boundaries_structure <- lapply(SU_key_HM, function(element) {
        lapply(element, function(SU) {
            SU$KP_boundaries_points
        })
    })

    # Loop through Typology
    for (id_typology in seq_along(SU_key_HM)) {
        # Get KP boundaries of all SU of the each typology
        SU_KP_boundaries_list <- SU_KP_boundaries_structure[[id_typology]]
        # Get HM KP values of the typology
        data_id_reach_HM <- Key_Info_Typology_Model_Reach[[id_typology]]$Model_Reach_nodes

        if (min(data_id_reach_HM) != min(unlist(SU_KP_boundaries_list))) stop(paste0("Min value of KP in the model (", min(data_id_reach_HM), ") is not presented in the boundaries of the spatialisation of the typology: ", names(SU_key_HM)[id_typology], ". Please verify the input Kmin or Kflood Key SU_MR"))

        if (max(data_id_reach_HM) != max(unlist(SU_KP_boundaries_list))) stop("Min value of KP in the model (", max(data_id_reach_HM), ") is not presented in the boundaries of the spatialisation of the typology: ", names(SU_key_HM)[id_typology], ". Please verify the input Kmin or Kflood Key SU_MR")

        RUG_min_boundary_KP <- min(Key_Info_Typology_Model_Reach[[id_typology]]$RUGFile$KP_start, Key_Info_Typology_Model_Reach[[id_typology]]$RUGFile$KP_end)
        RUG_max_boundary_KP <- max(Key_Info_Typology_Model_Reach[[id_typology]]$RUGFile$KP_start, Key_Info_Typology_Model_Reach[[id_typology]]$RUGFile$KP_end)

        # Check if KP boundaries of each SU respect the KP boundaries of XR
        if (min(unlist(SU_KP_boundaries_list)) != RUG_min_boundary_KP ||
            max(unlist(SU_KP_boundaries_list)) != RUG_max_boundary_KP) {
            return(stop(sprintf(
                "KP boundary mismatch: data min/max (%s, %s) vs HM min/max (%s, %s)",
                RUG_min_boundary_KP,
                RUG_max_boundary_KP,
                min(unlist(SU_KP_boundaries_list)),
                max(unlist(SU_KP_boundaries_list))
            )))
        }

        data_KP <- Key_Info_Typology_Model_Reach[[id_typology]]$KP_grid
        data_reaches <- Key_Info_Typology_Model_Reach[[id_typology]]$reach

        # Adding layer 3 : declare all SUs into the structure
        SU[[id_typology]] <- vector("list", length(SU_KP_boundaries_list))
        names(SU[[id_typology]]) <- names(SU_KP_boundaries_list)
        # Loop through each SU
        for (id_SU in seq_along(SU_KP_boundaries_list)) {
            # Adding layer 4: assign properties of each SU
            # Starting with KP and reach
            boundaries <- SU_KP_boundaries_list[[id_SU]]

            if (Key_Info_Typology_Model_Reach[[id_typology]]$Logical_decreasing) {
                # If logical_decreasing is TRUE
                # Include first value and exclude last if it is the last SU of the typology
                if (id_SU != length(SU_KP_boundaries_list)) {
                    position_match <- which(
                        data_KP <= boundaries[1] &
                            data_KP > boundaries[2]
                    )
                } else {
                    position_match <- which(
                        data_KP <= boundaries[1] &
                            data_KP >= boundaries[2]
                    )
                }
            } else {
                # Include start (left) and exclude end (right)
                if (id_SU != length(SU_KP_boundaries_list)) {
                    position_match <- which(
                        data_KP >= boundaries[1] &
                            data_KP < boundaries[2]
                    )
                } else {
                    position_match <- which(
                        data_KP >= boundaries[1] &
                            data_KP <= boundaries[2]
                    )
                }
            }

            SU[[id_typology]][[id_SU]]$KP <- data_KP[position_match]
            SU[[id_typology]][[id_SU]]$reach <- data_reaches[position_match]
            SU[[id_typology]][[id_SU]]$id_reach_SU_boundaries <- SU_KP_boundaries_list[[id_SU]]

            # Add Z
            if (identical(
                SU_key_HM[[id_typology]][[id_SU]]$function_SU,
                getCovariate_piecewise
            )) {
                if (!("shiftPoints" %in% names(SU_key_HM[[id_typology]][[id_SU]]))) {
                    return(stop(paste0(
                        "To apply getCovariate_piecewise function, shiftPoints must be passed as argument. Please check the XR = ", names(SU_key_HM)[[id_typology]], ", SU = ", names(SU_key_HM[[id_typology]])[[id_SU]]
                    )))
                }
                SU[[id_typology]][[id_SU]]$Z <- getCovariate_piecewise(
                    shiftPoints = SU_key_HM[[id_typology]][[id_SU]]$shiftPoints,
                    KP_grid = SU[[id_typology]][[id_SU]]$KP
                )

                df <- data.frame(
                    KP = SU[[id_typology]][[id_SU]]$KP,
                    SU[[id_typology]][[id_SU]]$Z
                )

                df$pattern <- apply(df[, names(df) != "KP"], 1, paste0, collapse = "")

                change_idx <- which(df$pattern != dplyr::lag(df$pattern))

                if (length(change_idx) != 0) {
                    KP_change <- df$KP[change_idx - 1]

                    SU[[id_typology]][[id_SU]]$id_reach_SU_boundaries <- c(
                        SU[[id_typology]][[id_SU]]$id_reach_SU_boundaries,
                        KP_change
                    )
                }
            } else if (identical(
                SU_key_HM[[id_typology]][[id_SU]]$function_SU,
                getCovariate_Legendre
            )) {
                if (!("max_polynomial_degree" %in% names(SU_key_HM[[id_typology]][[id_SU]]))) {
                    return(stop(paste0(
                        "To apply getCovariate_Legendre function, max_polynomial_degree must be passed as argument. Please check the XR = ", names(SU_key_HM)[[id_typology]], ", SU = ", names(SU_key_HM[[id_typology]])[[id_SU]]
                    )))
                }
                SU[[id_typology]][[id_SU]]$Z <- getCovariate_Legendre(
                    max_polynomial_degree = SU_key_HM[[id_typology]][[id_SU]]$max_polynomial_degree,
                    covariate_discretization = SU[[id_typology]][[id_SU]]$KP
                )
            } else {
                return(stop("function_SU given in input is not supported. Please select either getCovariate_Legendre or getCovariate_piecewise"))
            }

            # Add all boundaries of reach_SU
            SU[[id_typology]][[id_SU]]$id_reach_SU_boundaries <- sort(
                SU[[id_typology]][[id_SU]]$id_reach_SU_boundaries,
                decreasing = Key_Info_Typology_Model_Reach[[id_typology]]$Logical_decreasing
            )

            all_boundaries <- SU[[id_typology]][[id_SU]]$id_reach_SU_boundaries

            # Add id_reach_SU
            if (identical(
                SU_key_HM[[id_typology]][[id_SU]]$function_SU,
                getCovariate_piecewise
            )) {
                SU[[id_typology]][[id_SU]]$id_reach_SU <- sapply(
                    SU[[id_typology]][[id_SU]]$KP,
                    function(x) {
                        bounds <- all_boundaries
                        n <- length(bounds) - 1

                        idx <- sapply(seq_len(n), function(i) {
                            left <- bounds[i]
                            right <- bounds[i + 1]
                            # Decreasing
                            if (Key_Info_Typology_Model_Reach[[id_typology]]$Logical_decreasing) {
                                # Case 1: first interval → closed both sides
                                if (i == 1) {
                                    cond <- x <= left & x >= right

                                    # Case 2: open left, closed right
                                } else {
                                    cond <- x < left & x >= right
                                }
                            } else {
                                # Case 1: first interval → closed both sides
                                if (i == 1) {
                                    cond <- x >= left & x <= right

                                    # Case 2: open left, closed right
                                } else {
                                    cond <- x > left & x <= right
                                }
                            }

                            cond
                        })

                        idx <- which(idx)
                        if (length(idx) == 0) {
                            print(list(x = x, bounds = bounds))
                            stop("No interval match for KP = ", x)
                        }
                        # Safety checks
                        if (length(idx) == 0) {
                            stop("No interval match for KP = ", x)
                        } else if (length(idx) > 1) {
                            stop("KP belongs to multiple intervals: ", x)
                        } else {
                            idx
                        }
                    }
                )
            } else if (identical(
                SU_key_HM[[id_typology]][[id_SU]]$function_SU,
                getCovariate_Legendre
            )) {
                SU[[id_typology]][[id_SU]]$id_reach_SU <- rep(id_SU, length(SU[[id_typology]][[id_SU]]$KP))
            } else {
                return(stop("function_SU given in input is not supported. Please select either getCovariate_Legendre or getCovariate_piecewise"))
            }

            # Add priors
            SU[[id_typology]][[id_SU]]$prior <- SU_key_HM[[id_typology]][[id_SU]]$prior

            if (dim(SU[[id_typology]][[id_SU]]$Z)[2] != length(SU[[id_typology]][[id_SU]]$prior)) {
                stop(paste0(
                    "Error identified in XR = ", names(SU_key_HM)[[id_typology]], ", SU = ", names(SU_key_HM[[id_typology]])[[id_SU]], ". Number of columns of Z (", dim(SU[[id_typology]][[id_SU]]$Z)[2], ") must be equal to the length(prior) = ", length(SU[[id_typology]][[id_SU]]$prior)
                ))
            }
        }
    }

    # Summary SU
    summary_SU <- do.call(
        rbind,
        lapply(names(SU), function(section_name) {
            section <- SU[[section_name]]

            do.call(
                rbind,
                lapply(names(section), function(su_name) {
                    SU_local <- section[[su_name]]
                    id_vec <- SU_local$id_reach_SU_boundaries

                    data.frame(
                        typology = section_name,
                        id_SU = su_name,
                        id_reach_SU = seq_len(length(id_vec) - 1),
                        KP_max_SU = max(id_vec),
                        KP_min_SU = min(id_vec)
                    )
                })
            )
        })
    )

    return(list(
        SU = SU,
        summary_SU = summary_SU
    ))
}
