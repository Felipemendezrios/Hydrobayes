############################################
# Module 5 : BaM environment
############################################

# 1st layer: Parameters for vertical estimation of friction
####################################################################

# In MAGE context:
# Constant piecewise function : Kmin and Kflood

# In Dassflow context:
# K = alpha * (H) ^ beta
# Exponential function : so, alpha and beta are the parameters

############################################
# Kmin environment (encapsulated in Kmin_SU)
############################################
# Second layer: Set all SUs to each Typology by Kmin or Kflood
# Third layer: Set properties of each SU in terms of KP and spatialisation function
####################################################################
# Key to relate all SU by XR
Input_Kmin_Key_SU_MR <- list(
    # XR information
    MR =
        list(
            SU1 = list(
                # KP boundary points
                KP_boundaries_points = c(0, 25),
                # Function to apply at this SU
                function_SU = getCovariate_Legendre,
                # Arguments of this SU
                max_polynomial_degree = 6,
                prior = list(
                    RBaM::parameter(
                        name = "Kmin_MC_Up_a0",
                        init = 34
                    ),
                    RBaM::parameter(
                        name = "Kmin_MC_Up_a1",
                        init = 0
                    ),
                    RBaM::parameter(
                        name = "Kmin_MC_Up_a2",
                        init = 0
                    ),
                    RBaM::parameter(
                        name = "Kmin_MC_Up_a3",
                        init = 0
                    ),
                    RBaM::parameter(
                        name = "Kmin_MC_Up_a4",
                        init = 0
                    ),
                    RBaM::parameter(
                        name = "Kmin_MC_Up_a5",
                        init = 0
                    ),
                    RBaM::parameter(
                        name = "Kmin_MC_Up_a6",
                        init = 0
                    )
                )
            )
        ),
    # XR information
    TR =
        list(
            SU1 = list(
                # KP boundary points
                KP_boundaries_points = c(0, 20),
                # Function to apply at this SU
                function_SU = getCovariate_Legendre,
                # Arguments of this SU
                max_polynomial_degree = 6,
                prior = list(
                    RBaM::parameter(
                        name = "Kmin_TR_a0",
                        init = 28
                    ),
                    RBaM::parameter(
                        name = "Kmin_TR_a1",
                        init = 0
                    ),
                    RBaM::parameter(
                        name = "Kmin_TR_a2",
                        init = 0
                    ),
                    RBaM::parameter(
                        name = "Kmin_TR_a3",
                        init = 0
                    ),
                    RBaM::parameter(
                        name = "Kmin_TR_a4",
                        init = 0
                    ),
                    RBaM::parameter(
                        name = "Kmin_TR_a5",
                        init = 0
                    ),
                    RBaM::parameter(
                        name = "Kmin_TR_a6",
                        init = 0
                    )
                )
            )
        )
)
############################################
# End Kmin environment
############################################

############################################
# Kflood environment (encapsulated in Kflood_SU)
############################################
# Second layer: Set all SUs to each Typology by Kmin or Kflood
# Third layer: Set properties of each SU in terms of KP and spatialisation function
####################################################################
# Key to relate all SU by XR
Input_Kflood_Key_SU_MR <- list(
    # XR information
    MR =
        list(
            SU1 = list(
                # KP boundary points
                KP_boundaries_points = c(0, 25),
                # Function to apply at this SU
                function_SU = getCovariate_Legendre,
                # Arguments of this SU
                max_polynomial_degree = 0,
                prior = list(
                    RBaM::parameter(
                        name = "Kf_MC_Up_a0",
                        init = 18,
                        prior.dist = "FIX"
                    )
                )
            )
        ),
    # XR information
    TR =
        list(
            SU1 = list(
                # KP boundary points
                KP_boundaries_points = c(0, 20),
                # Function to apply at this SU
                function_SU = getCovariate_Legendre,
                # Arguments of this SU
                max_polynomial_degree = 0,
                prior = list(
                    RBaM::parameter(
                        name = "Kf_TR_a0",
                        init = 18,
                        prior.dist = "FIX"
                    )
                )
            )
        )
)

############################################
# End Kflood environment
############################################

############################################
# End module 5: BaM environment
############################################
