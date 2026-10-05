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
    Durance =
        list(
            SU1 = list(
                # KP boundary points
                KP_boundaries_points = c(0, 1079.063436),
                # Function to apply at this SU
                function_SU = getCovariate_Legendre,
                # Arguments of this SU
                max_polynomial_degree = 8,
                prior = list(
                    name_init = data.frame(
                        # Name to appear into calculations
                        name = c(
                            "Km_a0",
                            "Km_a1",
                            "Km_a2",
                            "Km_a3",
                            "Km_a4",
                            "Km_a5",
                            "Km_a6",
                            "Km_a7",
                            "Km_a8"
                        ),
                        # Initial guess
                        init = c(
                            30,
                            rep(0, 8)
                        )
                    ),
                    config = list(
                        # Prior distribution: either Gaussian or FIX
                        distribution = "Gaussian",
                        # Specific spatially coordinates from the SU to apply the prior information. NULL indicates that algorithm will be used to spread as best as possible
                        x_spatial = c(
                            NULL
                        ),
                        # A only value is accepted and assigned to all the x_spatial.
                        param_values = data.frame(
                            mu = c(
                                30
                            ),
                            sigma = c(
                                12
                            )
                        )
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
    Durance =
        list(
            SU1 = list(
                # KP boundary points
                KP_boundaries_points = c(0, 1079.063436),
                # Function to apply at this SU
                function_SU = getCovariate_Legendre,
                # Arguments of this SU
                max_polynomial_degree = 0,
                prior = list(
                    name_init = data.frame(
                        name = c(
                            "Kf_SU1_Rh_a0"
                        ),
                        init = 30
                    ),
                    config = list(
                        # Prior distribution: either Gaussian or FIX
                        distribution = "FIX", # FIX, no prior is given. All information below is skipped
                        # Specific spatially coordinates from the SU to apply the prior information. NULL indicates that algorithm will be used to spread as best as possible
                        x_spatial = c(
                            NULL
                        ),
                        # A only value is accepted and assigned to all the x_spatial.
                        param_values = data.frame(
                            mu = c(
                                30
                            ),
                            sigma = c(
                                8
                            )
                        )
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
