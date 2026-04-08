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
                KP_boundaries_points = c(0, 1079.063436),
                # Function to apply at this SU
                function_SU = getCovariate_Legendre,
                # Arguments of this SU
                max_polynomial_degree = 1,
                prior = list(
                    RBaM::parameter(
                        name = "Kmin_a0",
                        init = 30
                    ),
                    RBaM::parameter(
                        name = "Kmin_a1",
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
                KP_boundaries_points = c(0, 1079.063436),
                # Function to apply at this SU
                function_SU = getCovariate_Legendre,
                # Arguments of this SU
                max_polynomial_degree = 0,
                prior = list(
                    RBaM::parameter(
                        name = "Kf_a0",
                        init = 30,
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
# Qinflow (perturbations):
# Perturbation could be different throught events
# The order and the number of Q(t) does not change between events
############################################

# Fifth layer:
# Set the prior information of mulitiplifative factor
# be careful with the order of the event and the nodes
# Prior associated to each dischage time series
############################################

# Default values :
# mean = 1 -> no perturbation (centered to the value)
# sd = 0.1 -> 10% of perturbation
mult_factor <- list(
    # Get order of the event from all_events defined in main file
    event_1 = list(
        # Get order of the nodes inside the .HYD file, must be sure that nodes positions are the same trhougth the events
        node_1 = RBaM::parameter(
            name = "Q_mult_e1_n1",
            init = 1,
            prior.dist = "FIX"
            # prior.par = c(0, 0.1)
        )
    )
)
############################################
# End Qinflow (perturbations)
############################################

############################################
# End module 5: BaM environment
############################################
