RUGFile_constructor_SR_correction <- function(Key_Info_XR_MR) {
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


SR_constructor <- function(SR_key_HM,
                           Key_Info_Typology_Model_Reach) {
    # Initialize structure by fixing layer 1 : XR
    # Layer 2 : Kmin or Kmoy is handle outside this function
    SR <- vector(mode = "list", length = length(SR_key_HM))
    names(SR) <- names(SR_key_HM)

    # Extract only KP boundaries points keeping the same structure of XR
    SR_KP_boundaries_structure <- lapply(SR_key_HM, function(element) {
        lapply(element, function(SR) {
            SR$KP_boundaries_points
        })
    })

    # Loop through XR
    for (id_XR in seq_along(SR_key_HM)) {
        # Get KP boundaries of all SR at a XR
        SR_KP_boundaries_list <- SR_KP_boundaries_structure[[id_XR]]
        RUG_min_boundary_KP <- min(Key_Info_Typology_Model_Reach[[id_XR]]$RUGFile$KP_start)
        RUG_max_boundary_KP <- max(Key_Info_Typology_Model_Reach[[id_XR]]$RUGFile$KP_end)

        # Check if KP boundaries of each SR respect the KP boundaries of XR
        if (min(unlist(SR_KP_boundaries_list)) != RUG_min_boundary_KP ||
            max(unlist(SR_KP_boundaries_list)) != RUG_max_boundary_KP) {
            return(stop(sprintf(
                "KP boundary mismatch: data min/max (%s, %s) vs HM min/max (%s, %s)",
                RUG_min_boundary_KP,
                RUG_max_boundary_KP,
                min(unlist(SR_KP_boundaries_list)),
                max(unlist(SR_KP_boundaries_list))
            )))
        }
        data_KP <- Key_Info_Typology_Model_Reach[[id_XR]]$KP_grid
        data_reaches <- Key_Info_Typology_Model_Reach[[id_XR]]$reach

        # Adding layer 3 : declare all SRs into the structure
        SR[[id_XR]] <- vector("list", length(SR_KP_boundaries_list))
        names(SR[[id_XR]]) <- names(SR_KP_boundaries_list)
        # Loop through each SR
        for (id_SR in seq_along(SR_KP_boundaries_list)) {
            # Adding layer 4: assign properties of each SR
            # Starting with KP and reach
            boundaries <- SR_KP_boundaries_list[[id_SR]]

            if (Key_Info_Typology_Model_Reach[[id_XR]]$Logical_decreasing) {
                # If logical_decreasing is TRUE, boundaries must be reordered to the calculation
                boundaries <- sort(boundaries)
                # Include end (right) and exclude start (left)
                if (id_SR != length(SR_KP_boundaries_list)) {
                    position_match <- which(
                        data_KP > boundaries[1] &
                            data_KP <= boundaries[2]
                    )
                } else {
                    position_match <- which(
                        data_KP >= boundaries[1] &
                            data_KP <= boundaries[2]
                    )
                }
            } else {
                # Include start (left) and exclude end (right)
                if (id_SR != length(SR_KP_boundaries_list)) {
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

            SR[[id_XR]][[id_SR]]$KP <- data_KP[position_match]
            SR[[id_XR]][[id_SR]]$reach <- data_reaches[position_match]

            # Add Z and priors
            if (identical(
                SR_key_HM[[id_XR]][[id_SR]]$function_SR,
                getCovariate_piecewise
            )) {
                if (!("shiftPoints" %in% names(SR_key_HM[[id_XR]][[id_SR]]))) {
                    return(stop(paste0(
                        "To apply getCovariate_piecewise function, shiftPoints must be passed as argument. Please check the XR = ", names(SR_key_HM)[[id_XR]], ", SR = ", names(SR_key_HM[[id_XR]])[[id_SR]]
                    )))
                }
                SR[[id_XR]][[id_SR]]$Z <- getCovariate_piecewise(
                    shiftPoints = SR_key_HM[[id_XR]][[id_SR]]$shiftPoints,
                    KP_grid = SR[[id_XR]][[id_SR]]$KP
                )
            } else if (identical(
                SR_key_HM[[id_XR]][[id_SR]]$function_SR,
                getCovariate_Legendre
            )) {
                if (!("max_polynomial_degree" %in% names(SR_key_HM[[id_XR]][[id_SR]]))) {
                    return(stop(paste0(
                        "To apply getCovariate_Legendre function, max_polynomial_degree must be passed as argument. Please check the XR = ", names(SR_key_HM)[[id_XR]], ", SR = ", names(SR_key_HM[[id_XR]])[[id_SR]]
                    )))
                }
                SR[[id_XR]][[id_SR]]$Z <- getCovariate_Legendre(
                    max_polynomial_degree = SR_key_HM[[id_XR]][[id_SR]]$max_polynomial_degree,
                    covariate_discretization = SR[[id_XR]][[id_SR]]$KP
                )
            } else {
                return(stop("function_SR given in input is not supported. Please select either getCovariate_Legendre or getCovariate_piecewise"))
            }

            SR[[id_XR]][[id_SR]]$prior <- SR_key_HM[[id_XR]][[id_SR]]$prior

            if (dim(SR[[id_XR]][[id_SR]]$Z)[2] != length(SR[[id_XR]][[id_SR]]$prior)) {
                stop(paste0(
                    "Error identified in XR = ", names(SR_key_HM)[[id_XR]], ", SR = ", names(SR_key_HM[[id_XR]])[[id_SR]], ". Number of columns of Z (", dim(SR[[id_XR]][[id_SR]]$Z)[2], ") must be equal to the length(prior) = ", length(SR[[id_XR]][[id_SR]]$prior)
                ))
            }
        }
    }
    return(SR)
}
