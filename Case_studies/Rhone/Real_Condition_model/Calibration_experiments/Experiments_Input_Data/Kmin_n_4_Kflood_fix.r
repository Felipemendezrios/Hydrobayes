############################################
# Module 5 : BaM environment
############################################

# First layer: XR
############################################
# All reaches in the model (MR) defined within the same coordinate reference

# Must be careful with the order, the order of the reaches must be given upstream to downstream
# Input_XR <- list(
#     Rhone = c(1, 2, 3),
#     Ain = c(8, 4, 5),
#     Cassier = c(6, 7)
# )

# At this level, all information is already available to link BaM and HM
# Objective: common information to share in BaM and HM environment
# This step is so important to ensure the transmission of information

# Define the grid of interpolation by XR
Key_Info_XR_MR <- vector(mode = "list", length = length(Input_XR))
names(Key_Info_XR_MR) <- names(Input_XR)
for (i in seq_along(Key_Info_XR_MR)) {
    # Get information from MH by XR
    mask_MR_by_XR <- Input_MR[match(Input_XR[[i]], Input_MR$reach), ]
    # Input_MR[which(Input_MR$reach %in% Input_XR[[i]]), ]

    # Get nodes of MR
    MR_nodes <- unique(c(
        mask_MR_by_XR$KP_start,
        mask_MR_by_XR$KP_end
    ))

    Key_Info_XR_MR[[i]]$MR_nodes <- MR_nodes

    # Grid from interpolation passing by MR_nodes
    # KP_MR_nodes is defined by XR
    KP_MR_nodes <- interpolation_specific_points(
        total_points = 100,
        specific_nodes = MR_nodes
    )

    # Assign reach at the KP_MR_nodes
    KP_reach_MR_nodes <- assign_reach_from_a_grid(
        reach_KP_boundaries = mask_MR_by_XR,
        grid = KP_MR_nodes
    )

    # These information allow build RUGFile

    RUGFile <- data.frame(
        id_reach = KP_reach_MR_nodes$reaches[-nrow(KP_reach_MR_nodes)],
        KP_start = KP_reach_MR_nodes$grid[-nrow(KP_reach_MR_nodes)],
        KP_end = KP_reach_MR_nodes$grid[-1],
        Kmin = 0,
        Kflood = 0
    )

    Key_Info_XR_MR[[i]]$RUGFile <- RUGFile
    # Now, the middle value of each interval is taken to get the KP_grid, which will be used for all the calculations

    # Calculate the middle value for each row
    Key_Info_XR_MR[[i]]$KP_grid <- (RUGFile$KP_start + RUGFile$KP_end) / 2
    Key_Info_XR_MR[[i]]$reach <- RUGFile$id_reach
}

########################
# Get discretization of the covariate, same in Kmin or Kflood
########################
covariate_grid_lists <- lapply(Key_Info_XR_MR, function(SR) {
    SR$KP_grid
})
covariate_grid <- data.frame(covariate = unlist(covariate_grid_lists))
rownames(covariate_grid) <- NULL

########################
# RUGFile
########################
# Extract all RUGFile data frames from the main list
all_RUGFiles <- lapply(Key_Info_XR_MR, function(channel) {
    channel$RUGFile
})

# Row-bind all RUGFile data frames into a single data frame
RUGFile_Mage_Final <- do.call(rbind, all_RUGFiles) %>% arrange(id_reach)
# Remove row names
rownames(RUGFile_Mage_Final) <- NULL
#######################################
# End RUGFile
#######################################

# Second layer: Parameters for vertical estimation of friction
####################################################################

# In MAGE context:
# Constant piecewise function : Kmin and Kflood

# In Dassflow context:
# K = alpha * (H) ^ beta
# Exponential function : so, alpha and beta are the parameters

############################################
# Kmin environment (encapsulated in Kmin_SR)
############################################

# Third layer: Set all SRs to each XR by Kmin or Kflood
# Fourth layer: Set properties of each SR in terms of KP and spatialisation function
####################################################################
# Key to relate all SR by XR
Input_Kmin_Key_SR_MR <- list(
    # XR information
    Rhone =
        list(
            SR1 = list(
                # KP boundary points
                KP_boundaries_points = c(55900, 26750),
                # Function to apply at this SR
                function_SR = getCovariate_Legendre,
                # Arguments of this SR
                max_polynomial_degree = 4,
                prior = list(
                    RBaM::parameter(
                        name = "Kmin_Rh_SR1_a0",
                        init = 15,
                        prior.dist = "FlatPrior+"
                    ),
                    RBaM::parameter(
                        name = "Kmin_Rh_SR1_a1",
                        init = 0
                    ),
                    RBaM::parameter(
                        name = "Kmin_Rh_SR1_a2",
                        init = 0
                    ),
                    RBaM::parameter(
                        name = "Kmin_Rh_SR1_a3",
                        init = 0
                    ),
                    RBaM::parameter(
                        name = "Kmin_Rh_SR1_a4",
                        init = 0
                    )
                )
            )
        ),
    # XR information
    Ain =
        list(
            SR1 = list(
                # KP boundary points
                KP_boundaries_points = c(22333, 41461),
                # Function to apply at this SR
                function_SR = getCovariate_Legendre,
                # Arguments of this SR
                max_polynomial_degree = 4,
                prior = list(
                    RBaM::parameter(
                        name = "Kmin_Ai_SR1_a0",
                        init = 15,
                        prior.dist = "FlatPrior+"
                    ),
                    RBaM::parameter(
                        name = "Kmin_Ai_SR1_a1",
                        init = 0
                    ),
                    RBaM::parameter(
                        name = "Kmin_Ai_SR1_a2",
                        init = 0
                    ),
                    RBaM::parameter(
                        name = "Kmin_Ai_SR1_a3",
                        init = 0
                    ),
                    RBaM::parameter(
                        name = "Kmin_Ai_SR1_a4",
                        init = 0
                    )
                )
            )
        ),
    Cassier =
        list(
            SR1 = list(
                # KP boundary points
                KP_boundaries_points = c(0, 40),
                # Function to apply at this SR
                function_SR = getCovariate_Legendre,
                # Arguments of this SR
                max_polynomial_degree = 0,
                prior = list(
                    RBaM::parameter(
                        name = "Kmin_Cassier",
                        init = 35,
                        prior.dist = "FIX"
                    )
                )
            )
        )
)
# Check if size is respected between ID and Kmin
if (length(Input_XR) != length(Input_Kmin_Key_SR_MR)) stop("Size must be equal between Input_Kmin_Key_SR_MR and Input_XR")

# Assign properties of each SR in XR structure
Kmin_SR <- SR_constructor(
    SR_key_HM = Input_Kmin_Key_SR_MR,
    Key_Info_XR_MR = Key_Info_XR_MR
)

# Get Z file for all XR in Kmin
Z_all_Kmin_SR <- lapply(Kmin_SR, function(channel) {
    lapply(channel, function(SR) {
        SR$Z
    })
})

# Get the reaches of the Z file for all XR in Kmin
reaches_Z_all_Kmin_SR <- lapply(Kmin_SR, function(channel) {
    lapply(channel, function(SR) {
        SR$reach
    })
})

# Flatten the nested list of matrices
flat_Kmin_SR <- unlist(Z_all_Kmin_SR, recursive = FALSE)
reaches_flat_Kmin_SR <- cbind(unlist(unlist(reaches_Z_all_Kmin_SR, recursive = FALSE)))
rownames(reaches_flat_Kmin_SR) <- NULL

# Pass the flattened list to block_diagonal_matrix
Z_MatrixKmin_not_ordered <- do.call(block_diagonal_matrix, flat_Kmin_SR)

Z_MatrixKmin <- as.matrix(
    data.frame(
        cbind(
            reaches_flat_Kmin_SR,
            Z_MatrixKmin_not_ordered
        )
    ) %>%
        arrange(X1) %>%
        select(-X1)
)

# Check size between RUGFile_Mage_Final and Z_file
if (nrow(RUGFile_Mage_Final) != nrow(Z_MatrixKmin)) stop("RUGFile_Mage_Final must have the same size as Z file (spatialisation)")

Kmin_prior <- extract_priors(Kmin_SR)
############################################
# End Kmin environment
############################################

############################################
# Kflood environment (encapsulated in Kflood_SR)
############################################

# Third layer: Set all SRs to each XR by Kflood
# Fourth layer: Set properties of each SR in terms of KP and spatialisation function
####################################################################
# Key to relate all SR by XR
Input_Kflood_Key_SR_MR <- list(
    # XR information
    Rhone =
        list(
            SR1 = list(
                # KP boundary points
                KP_boundaries_points = c(55900, 26750),
                # Function to apply at this SR
                function_SR = getCovariate_Legendre,
                # Arguments of this SR
                max_polynomial_degree = 0,
                prior = list(
                    RBaM::parameter(
                        name = "Kflood_Rh_SR1_a0",
                        init = 15,
                        prior.dist = "FIX"
                    )
                )
            )
        ),
    Ain =
        list(
            SR1 = list(
                # KP boundary points
                KP_boundaries_points = c(22333, 41461),
                # Function to apply at this SR
                function_SR = getCovariate_Legendre,
                # Arguments of this SR
                max_polynomial_degree = 0,
                prior = list(
                    RBaM::parameter(
                        name = "Kflood_Ai_SR1_a0",
                        init = 15,
                        prior.dist = "FIX"
                    )
                )
            )
        ),
    Cassier =
        list(
            SR1 = list(
                # KP boundary points
                KP_boundaries_points = c(0, 40),
                # Function to apply at this SR
                function_SR = getCovariate_Legendre,
                # Arguments of this SR
                max_polynomial_degree = 0,
                prior = list(
                    RBaM::parameter(
                        name = "Kflood_Cassier",
                        init = 15,
                        prior.dist = "FIX"
                    )
                )
            )
        )
)

# Check if size is respected between ID and Kmin
if (length(Input_XR) != length(Input_Kflood_Key_SR_MR)) stop("Size must be equal between Input_Kflood_Key_SR_MR and Input_XR")

# Assign properties of each SR in XR structure
Kflood_SR <- SR_constructor(
    SR_key_HM = Input_Kflood_Key_SR_MR,
    Key_Info_XR_MR = Key_Info_XR_MR
)

# Get Z file for all XR in Kflood
Z_all_Kflood_SR <- lapply(Kflood_SR, function(channel) {
    lapply(channel, function(SR) {
        SR$Z
    })
})

# Get the reaches of the Z file for all XR in Kmin
reaches_Z_all_Kflood_SR <- lapply(Kflood_SR, function(channel) {
    lapply(channel, function(SR) {
        SR$reach
    })
})

# Flatten the nested list of matrices
flat_Kflood_SR <- unlist(Z_all_Kflood_SR, recursive = FALSE)
reaches_flat_Kflood_SR <- cbind(unlist(unlist(reaches_Z_all_Kflood_SR, recursive = FALSE)))
rownames(reaches_flat_Kflood_SR) <- NULL

# Pass the flattened list to block_diagonal_matrix
Z_MatrixKflood_not_ordered <- do.call(block_diagonal_matrix, flat_Kflood_SR)


Z_MatrixKflood <- as.matrix(
    data.frame(
        cbind(
            reaches_flat_Kflood_SR,
            Z_MatrixKflood_not_ordered
        )
    ) %>%
        arrange(X1) %>%
        select(-X1)
)



############


# Check size between RUGFile_Mage_Final and Z_file
if (nrow(RUGFile_Mage_Final) != nrow(Z_MatrixKflood)) stop("RUGFile_Mage_Final must have the same size as Z file (spatialisation)")

# Get discretization of the covariate
covariate_Kflood_lists <- lapply(Kflood_SR, function(channel) {
    lapply(channel, function(SR) {
        SR$KP
    })
})

Kflood_prior <- extract_priors(Kflood_SR)

############################################
# End Kflood environment
############################################

if (!(nrow(covariate_grid) == nrow(Z_MatrixKmin) && nrow(covariate_grid) == nrow(Z_MatrixKflood))) {
    stop("Number of rows of covariate_grid must be equal to both Z_MatrixKmin and Z_MatrixKflood")
}
############################################
# End module 5: BaM environment
############################################
