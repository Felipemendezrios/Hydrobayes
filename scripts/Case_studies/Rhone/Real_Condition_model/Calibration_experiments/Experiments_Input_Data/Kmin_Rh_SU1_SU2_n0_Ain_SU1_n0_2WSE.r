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

# Second layer: Set all SUs to each MR by Kmin or Kflood
# Third layer: Set properties of each SU in terms of KP and spatialisation function
####################################################################
# Key to relate all SU by MR
Input_Kmin_Key_SU_MR <- list(
    # MR information
    Rhone =
        list(
            SU1 = list(
                # KP boundary points
                KP_boundaries_points = c(55900, 34500),
                # Function to apply at this SU
                function_SU = getCovariate_Legendre,
                # Arguments of this SU
                max_polynomial_degree = 0,
                prior = list(
                    RBaM::parameter(
                        name = "Km_SU1_Rh_a0",
                        init = 36,
                        prior.dist = "FlatPrior+"
                    )
                )
            ),
            SU2 = list(
                # KP boundary points
                KP_boundaries_points = c(34500, 26750),
                # Function to apply at this SU
                function_SU = getCovariate_Legendre,
                # Arguments of this SU
                max_polynomial_degree = 0,
                prior = list(
                    RBaM::parameter(
                        name = "Km_SU2_Rh_a0",
                        init = 31,
                        prior.dist = "FlatPrior+"
                    )
                )
            )
        ),
    # MR information
    Ain =
        list(
            SU1 = list(
                # KP boundary points
                KP_boundaries_points = c(22333, 41461),
                # Function to apply at this SU
                function_SU = getCovariate_Legendre,
                # Arguments of this SU
                max_polynomial_degree = 0,
                prior = list(
                    RBaM::parameter(
                        name = "Km_SU1_Ai_a0",
                        init = 19,
                        prior.dist = "FlatPrior+"
                    )
                )
            )
        ),
    Cassier =
        list(
            SU1 = list(
                # KP boundary points
                KP_boundaries_points = c(0, 40),
                # Function to apply at this SU
                function_SU = getCovariate_Legendre,
                # Arguments of this SU
                max_polynomial_degree = 0,
                prior = list(
                    RBaM::parameter(
                        name = "Km_SU1_Cas_a0",
                        init = 35,
                        prior.dist = "FIX"
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

# Third layer: Set all SUs to each MR by Kflood
# Fourth layer: Set properties of each SU in terms of KP and spatialisation function
####################################################################
# Key to relate all SU by MR
Input_Kflood_Key_SU_MR <- list(
    # MR information
    Rhone =
        list(
            SU1 = list(
                # KP boundary points
                KP_boundaries_points = c(55900, 26750),
                # Function to apply at this SU
                function_SU = getCovariate_Legendre,
                # Arguments of this SU
                max_polynomial_degree = 0,
                prior = list(
                    RBaM::parameter(
                        name = "Kf_SU1_Rh_a0",
                        init = 20,
                        prior.dist = "FIX"
                    )
                )
            )
        ),
    # MR information
    Ain =
        list(
            SU1 = list(
                # KP boundary points
                KP_boundaries_points = c(22333, 41461),
                # Function to apply at this SU
                function_SU = getCovariate_Legendre,
                # Arguments of this SU
                max_polynomial_degree = 0,
                prior = list(
                    RBaM::parameter(
                        name = "Kf_SU1_Ai_a0",
                        init = 15,
                        prior.dist = "FIX"
                    )
                )
            )
        ),
    Cassier =
        list(
            SU1 = list(
                # KP boundary points
                KP_boundaries_points = c(0, 40),
                # Function to apply at this SU
                function_SU = getCovariate_Legendre,
                # Arguments of this SU
                max_polynomial_degree = 0,
                prior = list(
                    RBaM::parameter(
                        name = "Kf_SU1_Cas_a0",
                        init = 15,
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
# Set the prior information of mulitiplicative factor
# be careful with the order of the event and the nodes
# Prior associated to each dischage time series
############################################

# Default values :
# mean = 1 -> no perturbation (centered to the value)
# sd = 0.1 -> 10% of perturbation
mult_factor <- list(
    # Get order of the event from all_events defined in main file
    event_1 = list(
        # Get order of the nodes inside the .HYD file, must be sure that nodes positions are the same througth the events
        node_1 = RBaM::parameter(
            name = "Q_e1_CAIN",
            init = 1,
            prior.dist = "FIX"
            # prior.par = c(0, 0.1)
        ),
        node_2 = RBaM::parameter(
            name = "Q_e1_PCH",
            init = 1,
            prior.dist = "FIX"
        ),
        node_3 = RBaM::parameter(
            name = "Q_e1_LGN",
            init = 1,
            prior.dist = "FIX"
        ),
        node_4 = RBaM::parameter(
            name = "Q_e1_BOU",
            init = 1,
            prior.dist = "FIX"
        )
    ),
    event_2 = list(
        node_1 = RBaM::parameter(
            name = "Q_e2_CAIN",
            init = 1,
            prior.dist = "FIX"
        ),
        node_2 = RBaM::parameter(
            name = "Q_e2_PCH",
            init = 1,
            prior.dist = "FIX"
        ),
        node_3 = RBaM::parameter(
            name = "Q_e2_LGN",
            init = 1,
            prior.dist = "FIX"
        ),
        node_4 = RBaM::parameter(
            name = "Q_e2_BOU",
            init = 1,
            prior.dist = "FIX"
        )
    )
)
############################################
# End Qinflow (perturbations)
############################################

############################################
# End module 5: BaM environment
############################################
