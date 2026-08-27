############################################
# Module 5 : BaM environment
############################################

# First layer: XR
############################################
# All reaches in the model (MR) defined within the same coordinate reference

Input_XR <- list(
    main_channel = c(1)
)

# At this level, all information is already available to link BaM and HM
# Objective: common information to share in BaM and HM environment
# This step is so important to ensure the transmission of information

# Define the grid of interpolation by XR
Key_Info_XR_MR <- vector(mode = "list", length = length(Input_XR))
names(Key_Info_XR_MR) <- names(Input_XR)
for (i in seq_along(Key_Info_XR_MR)) {
    # Get information from MH by XR
    mask_MR_by_XR <- Input_MR[which(Input_MR$reach %in% Input_XR[[i]]), ]

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
    reach_MR_nodes <- assign_reach_from_a_grid(
        reach_KP_boundaries = mask_MR_by_XR,
        grid = KP_MR_nodes
    )
    # These information allow build RUGFile
    RUGFile <- data.frame(
        id_reach = reach_MR_nodes[-length(reach_MR_nodes)],
        KP_start = KP_MR_nodes[-length(KP_MR_nodes)],
        KP_end = KP_MR_nodes[-1],
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
RUGFile <- do.call(rbind, all_RUGFiles)
# Remove row names
rownames(RUGFile) <- NULL
#######################################
# End RUGFile
#######################################

# Second layer: Parameters for vertical estimation of friction
####################################################################

# In MAGE context:
# Constant piecewise function : Kmin and Kflood

# In Dashflow context:
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
    main_channel =
        list(
            SR1 = list(
                # KP boundary points
                KP_boundaries_points = c(0, 18),
                # Function to apply at this SR
                function_SR = getCovariate_Legendre,
                # Arguments of this SR
                max_polynomial_degree = 0,
                prior = list(
                    RBaM::parameter(name = "MC_Kmin_SR1_a0", init = 1 / 0.010, prior.dist = "FlatPrior+")
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
# Flatten the nested list of matrices
flat_Kmin_SR <- unlist(Z_all_Kmin_SR, recursive = FALSE)

# Pass the flattened list to block_diagonal_matrix
Z_MatrixKmin <- do.call(block_diagonal_matrix, flat_Kmin_SR)

# Check size between RUGFile and Z_file
if (nrow(RUGFile) != nrow(Z_MatrixKmin)) stop("RugFile must have the same size as Z file (spatialisation)")

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
    main_channel =
        list(
            SR1 = list(
                # KP boundary points
                KP_boundaries_points = c(0, 18),
                # Function to apply at this SR
                function_SR = getCovariate_Legendre,
                # Arguments of this SR
                max_polynomial_degree = 0,
                prior = list(
                    RBaM::parameter(name = "MC_Kflood_SR1_a0", init = ((33 + 1 / 0.013) / 2 + 1 / 0.010) / 2, prior.dist = "FlatPrior+")
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
# Flatten the nested list of matrices
flat_Kflood_SR <- unlist(Z_all_Kflood_SR, recursive = FALSE)

# Pass the flattened list to block_diagonal_matrix
Z_MatrixKflood <- do.call(block_diagonal_matrix, flat_Kflood_SR)

# Check size between RUGFile and Z_file
if (nrow(RUGFile) != nrow(Z_MatrixKflood)) stop("RugFile must have the same size as Z file (spatialisation)")

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
