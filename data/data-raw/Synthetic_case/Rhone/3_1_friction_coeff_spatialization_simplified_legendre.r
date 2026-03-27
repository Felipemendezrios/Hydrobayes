rm(list = ls())
graphics.off()

function_list <- list.files("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/R", full.names = TRUE)
for (i in function_list) {
    source(i)
}
#####################

Input_MR <- data.frame(
    reach = c(1, 2, 3),
    KP_start = c(55900, 34500, 19330),
    KP_end = c(34500, 26750, 202)
)

Input_XR <- list(
    MR = c(1, 2),
    TR = c(3)
)

library(ggplot2)
library(dplyr)

Key_Info_XR_MR <- vector(mode = "list", length = length(Input_XR))
names(Key_Info_XR_MR) <- names(Input_XR)
for (i in seq_along(Key_Info_XR_MR)) {
    # Get information from MH by XR
    mask_MR_by_XR <- Input_MR[match(Input_XR[[i]], Input_MR$reach), ]
    # Input_MR[which(Input_MR$reach %in% Input_XR[[i]]), ]

    if (nrow(mask_MR_by_XR) != 1) {
        if (all(diff(mask_MR_by_XR$KP_start) < 0)) {
            Logical_decreasing <- TRUE

            mask_MR_by_XR <- mask_MR_by_XR %>%
                rename(
                    KP_start = KP_end,
                    KP_end = KP_start
                )
        } else if (all(diff(mask_MR_by_XR$KP_start) > 0)) {
            Logical_decreasing <- FALSE
        } else {
            stop("All KP_start values must be increased or decreased")
        }
    } else if (nrow(mask_MR_by_XR) == 0) {
        stop("Any matched found between MR and XR, please verify")
    } else {
        if ((mask_MR_by_XR$KP_start - mask_MR_by_XR$KP_end) < 0) {
            Logical_decreasing <- FALSE
        } else {
            Logical_decreasing <- TRUE

            mask_MR_by_XR <- mask_MR_by_XR %>%
                rename(
                    KP_start = KP_end,
                    KP_end = KP_start
                )
        }
    }


    # Get nodes of MR
    MR_nodes <- unique(c(
        mask_MR_by_XR$KP_start,
        mask_MR_by_XR$KP_end
    ))
    MR_nodes <- sort(MR_nodes, decreasing = Logical_decreasing)

    Key_Info_XR_MR[[i]]$MR_nodes <- MR_nodes

    # Grid from interpolation passing by MR_nodes
    # KP_MR_nodes is defined by XR
    KP_MR_nodes <- interpolation_specific_points(
        total_points = 15,
        specific_nodes = MR_nodes
    )

    KP_MR_nodes <- sort(KP_MR_nodes, decreasing = Logical_decreasing)

    # Assign reach at the KP_MR_nodes
    KP_reach_MR_nodes <- assign_reach_from_a_grid(
        reach_KP_boundaries = mask_MR_by_XR,
        grid = KP_MR_nodes,
        Logical_decreasing = Logical_decreasing
    )

    # These information allow build RUGFile
    if (Logical_decreasing) {
        RUGFile <- data.frame(
            id_reach = KP_reach_MR_nodes$reaches[-nrow(KP_reach_MR_nodes)],
            KP_start = KP_reach_MR_nodes$grid[-1],
            KP_end = KP_reach_MR_nodes$grid[-nrow(KP_reach_MR_nodes)],
            Kmin = 0,
            Kflood = 0
        )
    } else {
        RUGFile <- data.frame(
            id_reach = KP_reach_MR_nodes$reaches[-nrow(KP_reach_MR_nodes)],
            KP_start = KP_reach_MR_nodes$grid[-nrow(KP_reach_MR_nodes)],
            KP_end = KP_reach_MR_nodes$grid[-1],
            Kmin = 0,
            Kflood = 0
        )
    }
    Key_Info_XR_MR[[i]]$RUGFile <- RUGFile
    # Now, the middle value of each interval is taken to get the KP_grid, which will be used for all the calculations

    # Calculate the middle value for each row
    Key_Info_XR_MR[[i]]$KP_grid <- (RUGFile$KP_start + RUGFile$KP_end) / 2
    Key_Info_XR_MR[[i]]$reach <- RUGFile$id_reach
    Key_Info_XR_MR[[i]]$Logical_decreasing <- Logical_decreasing
}


##########
kmin_spat <- list(
    MR = c(35, -0.5, 4, 3),
    TR = c(30, -1.1, -3, 3.5)
)

kflood_spat <- list(
    MR = c(15),
    TR = c(15)
)


# toto <- list()
# toto[[1]] <- getCovariate_piecewise(
#     KP_grid = Key_Info_XR_MR[[1]]$KP_grid,
#     shiftPoints = c(47000, 34500)
# )

# toto[[2]] <- getCovariate_piecewise(
#     KP_grid = Key_Info_XR_MR[[2]]$KP_grid,
#     shiftPoints = c(12500)
# )

HM_Ks_results <- spatialisation_results <- list()
for (i in 1:length(Key_Info_XR_MR)) {
    Z <- getCovariate_Legendre(
        max_polynomial_degree = 3,
        covariate_discretization = Key_Info_XR_MR[[i]]$KP_grid
    )

    Zflood <- getCovariate_Legendre(
        max_polynomial_degree = 0,
        covariate_discretization = Key_Info_XR_MR[[i]]$KP_grid
    )

    kmin_spat_res <- Z %*% kmin_spat[[i]]
    kflood_spat_res <- Zflood %*% kflood_spat[[i]]


    plot(kmin_spat_res)
    plot(kflood_spat_res)
    spatialisation_results[[i]] <- data.frame(
        id_reach = i,
        KP_start = Key_Info_XR_MR[[i]]$RUGFile$KP_start,
        KP_end = Key_Info_XR_MR[[i]]$RUGFile$KP_end,
        Kmin = kmin_spat_res,
        Kflood = kflood_spat_res,
        K_covariate = Key_Info_XR_MR[[i]]$KP_grid
    )
    # # SU:
    # RUG_path <- paste0("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/data-raw/Synthetic_case/Rhone/PamHyr_Legendre/Calibration/RUGFiles/RUGFile_spatialized_SU_", i, ".RUG")

    # write_RUGFile(
    #     RUG_path = RUG_path,
    #     RUGFile_data = spatialisation_results[[i]],
    #     RUG_format = "%1s%3d      %10.3f%10.3f%10.2f%10.2f"
    # )

    # MAGE:
    RUG_path_HM <- paste0("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/data-raw/Synthetic_case/Rhone/PamHyr_Legendre/Calibration/RUGFiles/RUGFile_spatialized_SU_by_HM_reach_", i, ".RUG")

    HM_Ks_results[[i]] <- data.frame(
        id_reach = Key_Info_XR_MR[[i]]$RUGFile$id_reach,
        KP_start = Key_Info_XR_MR[[i]]$RUGFile$KP_start,
        KP_end = Key_Info_XR_MR[[i]]$RUGFile$KP_end,
        Kmin = kmin_spat_res,
        Kflood = kflood_spat_res,
        K_covariate = Key_Info_XR_MR[[i]]$KP_grid
    )

    write_RUGFile(
        RUG_path = RUG_path_HM,
        RUGFile_data = HM_Ks_results[[i]],
        RUG_format = "%1s%3d      %10.3f%10.3f%10.2f%10.2f"
    )
}

RUGFile_data <- do.call(rbind, spatialisation_results)

RUG_path <- "/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/data-raw/Synthetic_case/Rhone/PamHyr_Legendre/Calibration/RUGFiles/RUGFile_spatialized_SU.RUG"


write_RUGFile(
    RUG_path = RUG_path,
    RUGFile_data = RUGFile_data,
    RUG_format = "%1s%3d      %10.3f%10.3f%10.2f%10.2f"
)


RUGFile_data <- do.call(rbind, HM_Ks_results)

RUG_path <- "/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/data-raw/Synthetic_case/Rhone/PamHyr_Legendre/Calibration/RUGFiles/RUGFile_spatialized_HM.RUG"


write_RUGFile(
    RUG_path = RUG_path,
    RUGFile_data = RUGFile_data,
    RUG_format = "%1s%3d      %10.3f%10.3f%10.2f%10.2f"
)
