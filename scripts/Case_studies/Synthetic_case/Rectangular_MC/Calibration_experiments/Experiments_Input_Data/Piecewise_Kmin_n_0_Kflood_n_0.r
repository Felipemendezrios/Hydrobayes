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

# Second layer: Set all SUs to each XR by Kmin or Kflood
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
                function_SU = getCovariate_piecewise,
                # Arguments of this SU
                shiftPoints = c(2.571, 7.714, 10.286, 15.429, 19.4),
                prior = lapply(
                    1:6,
                    function(i) {
                        RBaM::parameter(
                            name = paste0("Kmin_MC_a", i),
                            init = 34
                        )
                    }
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
                function_SU = getCovariate_piecewise,
                # Arguments of this SU
                shiftPoints = c(4.211, 6.316, 10.526, 13.684, 17.895),
                prior = lapply(
                    1:6,
                    function(i) {
                        RBaM::parameter(
                            name = paste0("Kmin_TR_a", i),
                            init = 28
                        )
                    }
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

# Third layer: Set all SUs to each XR by Kflood
# Fourth layer: Set properties of each SU in terms of KP and spatialisation function
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
