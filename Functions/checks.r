check_kp_first_reach <- function(MAGE_specific_points, Kmin_specific_ponits, Kflood_specific_ponits) {
    # Get the first reach's KP_start and KP_end for the lists
    kp_start_MAGE <- MAGE_specific_points$KP_start[1]
    kp_end_MAGE <- last(MAGE_specific_points$KP_end)

    kp_start_Kmin <- Kmin_specific_ponits$KP_start[1]
    kp_end_Kmin <- last(Kmin_specific_ponits$KP_end)

    kp_start_Kflood <- Kflood_specific_ponits$KP_start[1]
    kp_end_Kflood <- last(Kflood_specific_ponits$KP_end)

    # Check if any KP_start or KP_end values are not equal
    if (!all(kp_start_Kmin == kp_start_Kflood) ||
        !all(kp_start_Kmin == kp_start_MAGE) ||
        !all(kp_end_Kmin == kp_end_Kflood) ||
        !all(kp_end_Kmin == kp_end_MAGE)) {
        return(stop("KP_start or KP_end values in Kmin, Kflood, and MAGE_specific_points must be the same."))
    }
}


check_metadata <- function(metadata) {
    if (length(metadata$MAGE$spatial_specific_points$reach) != length(unique(metadata$MAGE$spatial_specific_points$reach))) {
        return(stop(paste0("Duplicated values detected in MAGE. Here the values ", metadata$MAGE$spatial_specific_points$reach)))
    }

    if (length(metadata$Kmin$spatial_specific_points$reach) != length(unique(metadata$Kmin$spatial_specific_points$reach))) {
        return(stop(paste0("Duplicated values detected in Kmin. Here the values ", metadata$Kmin$spatial_specific_points$reach)))
    }

    if (length(metadata$Kflood$spatial_specific_points$reach) != length(unique(metadata$Kflood$spatial_specific_points$reach))) {
        return(stop(paste0("Duplicated values detected in Kflood Here the values ", metadata$Kflood$spatial_specific_points$reach)))
    }

    if (nrow(metadata$MAGE$spatial_specific_points) > 1) {
        for (id in 1:(nrow(metadata$MAGE$spatial_specific_points) - 1)) {
            if (metadata$MAGE$spatial_specific_points$KP_end[id] != metadata$MAGE$spatial_specific_points$KP_start[id + 1]) {
                return(stop(paste0("The specified points must respect that KP_start of reach i must be equal to the KP_end of the reach i + 1. Please check the valmues in MAGE specifications. Here the values", metadata$MAGE$spatial_specific_points)))
            }
        }
    }
    if (nrow(metadata$Kmin$spatial_specific_points) > 1) {
        for (id in 1:(nrow(metadata$Kmin$spatial_specific_points) - 1)) {
            if (metadata$Kmin$spatial_specific_points$KP_end[id] != metadata$Kmin$spatial_specific_points$KP_start[id + 1]) {
                return(stop(paste0("The specified points must respect that KP_start of reach i must be equal to the KP_end of the reach i + 1. Please check the in Kmin specifications. Here the values", metadata$Kmin$spatial_specific_points)))
            }
        }
    }
    if (nrow(metadata$Kflood$spatial_specific_points) > 1) {
        for (id in 1:(nrow(metadata$Kflood$spatial_specific_points) - 1)) {
            if (metadata$Kflood$spatial_specific_points$KP_end[id] != metadata$Kflood$spatial_specific_points$KP_start[id + 1]) {
                return(stop(paste0("The specified points must respect that KP_start of reach i must be equal to the KP_end of the reach i + 1. Please check the in Kflood specifications. Here the values", metadata$Kflood$spatial_specific_points)))
            }
        }
    }
}

# check_prior(metadata){
#     if(metadata$spatialisation_function != "constant_piecewise_function()"){
#     if((metadata$n_degree+1) != length(metadata$prior) )stop('')

#     }
# }
