rm(list = ls())
graphics.off()

# Set directory
dir_workspace <- here::here()

# Load functions
function_list <- list.files("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/Functions/", full.names = TRUE)
for (i in function_list) {
    source(i)
}
# Processed data
load("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/processed_data/HHLab/compound_transition_friction/data_HHLab_uniform_case.RData")

# Libraries
library(RBaM)
library(dplyr)
library(patchwork)
library(tidyr)
library(ggplot2)
library(stringr)


# Set specific paths
dir_exe_BaM <- "/home/famendezrios/Documents/Git/BaM/makefile/"
MAGE_executable <- "/home/famendezrios/Documents/Softwares/pamhyr2/mage8/mage"

# Main path is specific to each study
file_main_path <- "/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/Case_studies/HHLab/Compound_channel_spatial_friction/Calibration_experiments"
mage_projet_name <- "CWMQ"

# Logical:
# Do or not the calibration
do_calibration <- FALSE
# Assign manually uncertainty values in meters
do_manual_uncertainty <- TRUE

############################################
# Observations data treatment:
############################################
# Measurements for calibration by event!
# Ensure that order should be in coherence with mage model set up!
Cal_measures <-
    list(
        all_events = c( # Give the names of the colonnes of observed data to compute. Order is important
            "CWMQ18", "CWMQ12"
        ),
        CWMQ18 =
            data.frame(
                event = 1,
                reach = 1,
                x = c(2.225, 4.310, 5.270, 6.240, 8.315, 9.115, 9.250, 10.250, 11.250, 12.250, 13.250),
                t = rep(3420, 11)
            ),
        CWMQ12 =
            data.frame(
                event = 2,
                reach = 1,
                x = c(2.225, 2.975, 3.345, 4.095, 5.060, 6.240, 6.990, 8.030, 8.795, 9.065, 9.115, 9.545, 10.030, 11.370, 12.120, 12.325, 13.290),
                t = rep(3420, 17)
            )
    )

if (!any(Cal_measures$all_events %in% names(Cal_measures))) stop(paste0("The all_events in Cal_measures should be the same in the other components of the list. Here, the event ", Cal_measures$all_events, " are not the same as the names of the list ", names(Cal_measures)[-1]))

# All available data
WSE_data_temp <- all_data_calibration$WSE

WSE_data_temp <- data.frame(WSE_data_temp %>%
    group_by(ID_experiment, x) %>%
    summarise(
        z_mean = mean(z_mean),
        Yu = mean(Yu),
        z_riverbed = mean(z_riverbed),
        .groups = "drop"
    ))

WSE_data_temp <- WSE_data_temp[str_detect(WSE_data_temp$ID_experiment, pattern = mage_projet_name), ]

# This part help to handle when observed data is presented as a list and user need to specify the the order of the events with the same name for matching with all_events, where the order is important

# Only giving the name keeping the data order and creating a link with the names of the event and the Cal_measures defined previously
match_case_event <- data.frame(
    idx_data = unique(levels(factor(WSE_data_temp$ID_experiment))),
    idx_event = unique(levels(factor(WSE_data_temp$ID_experiment)))
)

Link_x_t_ind_event <- lapply(names(Cal_measures)[-1], function(event) {
    # Find the row in match_case_event where idx_event matches the current event
    row_idx <- which(match_case_event$idx_event == event)

    list(
        idx_data = match_case_event$idx_data[row_idx],
        idx_event = event,
        Cal_measure = Cal_measures[[event]]
    )
})

## Initialize an empty list to store the results
observed_data <- data.frame()
X <- c()

# Loop over each element in Link_x_t_ind_event
for (i in seq_along(Link_x_t_ind_event)) {
    idx_data <- Link_x_t_ind_event[[i]]$idx_data
    cal_x <- Link_x_t_ind_event[[i]]$Cal_measure$x

    # Get the corresponding data frame from WSE_data_temp
    # wse_df <- WSE_data_temp[[idx_data]]
    wse_df <- WSE_data_temp[WSE_data_temp$ID_experiment == idx_data, ]

    # Filter rows where x is in cal_x
    observed_data <- rbind(
        observed_data,
        data.frame(
            event = Link_x_t_ind_event[[i]]$Cal_measure$event,
            wse_df[wse_df$x %in% cal_x, c("z_riverbed", "z_mean", "Yu")]
        )
    )
    # Extract the Cal_measure data frame
    X <- rbind(
        X,
        Link_x_t_ind_event[[i]]$Cal_measure
    )
}

if (!all(observed_data$event == X$event)) stop("The order of the event are not the same between the X data frame and observed_data data frame")


if (do_manual_uncertainty) {
    uncertainty_in_observation <- rep(0.005, length(observed_data$Yu))
} else {
    uncertainty_in_observation <- observed_data$Yu
}

Y <- data.frame(WSE = observed_data$z_mean)
Y$Discharge <- -9999
Y$Velocity <- -9999
Y$Y_Kmin <- -9999
Y$Y_Kmoy <- -9999

Yu <- data.frame(Yu_z = uncertainty_in_observation)
Yu$Yu_Q <- -9999
Yu$Yu_V <- -9999
Yu$Yu_Kmin <- -9999
Yu$Yu_Kmoy <- -9999


## Plot
CalData_plot_df <- data.frame(
    x = X$x,
    z_mean = Y$WSE * 1000,
    Yu = Yu$Yu_z * 1000,
    z_riverbed = observed_data$z_riverbed * 1000,
    ID = factor(X$event)
)

CalData_plot_export <- CalData_plot(
    data = CalData_plot_df,
    scales_free = ,
    y_label = "Water surface elevation (mm)",
    title_label = "WSE observations",
    col_label = "Events",
    plot_water_depth = FALSE,
    wrap = FALSE
)
############################################
# End of Observations data treatment:
############################################


############################################
# Calibration space
############################################


############################################
# Input parameters
############################################
# Reach 1 : main channel
a0_min_reach_1 <- RBaM::parameter(
    name = "a0_min_glass",
    init = 1 / 0.010, # Initial guess
    prior.dist = "FlatPrior+", # Prior distribution
    prior.par = NULL
) # Parameters of the prior distribution

# reach 1 : floodplain
a0_flood_reach_1_mean <- RBaM::parameter(
    name = "a0_flood",
    init = (33 + 1 / 0.013) / 2,
    prior.dist = "FlatPrior+",
    prior.par = NULL
)

a0_flood_reach_1_wood <- RBaM::parameter(
    name = "a0_flood_wood",
    init = 33,
    prior.dist = "FlatPrior+",
    prior.par = NULL
)

a0_flood_reach_1_grass <- RBaM::parameter(
    name = "a0_flood_grass",
    init = 1 / 0.013,
    prior.dist = "FlatPrior+",
    prior.par = NULL
)

a1_flood_reach_1 <- RBaM::parameter(
    name = "a1_flood",
    init = 0,
    prior.dist = "FlatPrior",
    prior.par = NULL
)

a2_flood_reach_1 <- RBaM::parameter(
    name = "a2_flood",
    init = 0,
    prior.dist = "FlatPrior",
    prior.par = NULL
)

a3_flood_reach_1 <- RBaM::parameter(
    name = "a3_flood",
    init = 0,
    prior.dist = "FlatPrior",
    prior.par = NULL
)
a4_flood_reach_1 <- RBaM::parameter(
    name = "a4_flood",
    init = 0,
    prior.dist = "FlatPrior",
    prior.par = NULL
)
############################################
# End input parameters
############################################

############################################
# Create the metadata to link reaches from MAGE to BaM! and inversely
############################################
# Remarks:
# BaM! could aggregate reaches from MAGE or slim that are higher than those from MAGE.
# BaM! reaches, called hereafter spatialisation unit (SU), are the referents for discretization of the points. However some specific points must be added to respect MAGE restrictions
# Metadata is created by SU, classifying in MAGE, Kmin and Kflood specific points. The last two components useful for spatiliasation

# metadata_settings description:
# number_points_interpolation: number of point to perform the interpolation


metadata_settings <- list(
    number_points_interpolation = 100
)

# How many main_SU exist in the MAGE projet?
# A main_SU separate the main channel to the tributaries
# A main_SU is composed by local_SU differently in minor bed and floodplain
number_of_all_main_SU <- 1

# Give information of each local_SU
# local_SU_1 <- list(
#     MAGE = list(
#         spatial_specific_points = data.frame(
#             reach = c(1),
#             KP_start = c(0),
#             KP_end = c(18)
#         )
#     ),
#     Kmin = list(
#         spatial_specific_points = data.frame(
#             reach = c(1),
#             KP_start = c(0),
#             KP_end = c(18)
#         ),
#         spatialisation_function = "constant_piecewise_function()",
#         n_degree = 0,
#         prior = list(
#             a0_min_reach_1
#         )
#     ),
#     Kflood = list(
#         spatial_specific_points = data.frame(
#             reach = c(1, 2),
#             KP_start = c(0, 9.05),
#             KP_end = c(9.05, 18)
#         ),
#         spatialisation_function = "constant_piecewise_function()",
#         n_degree = 0,
#         prior = list(
#             a0_flood_reach_1_wood,
#             a0_flood_reach_1_grass
#         )
#     )
# )


# Constructor for SU structure
new_su <- function(
    specific_points_MAGE,
    specific_points_BaM,
    spatialisation_function = NULL,
    n_degree = NULL,
    prior = NULL) {
    # Validate specific_points format
    required_cols <- c("reach", "KP_start", "KP_end")
    if (!all(required_cols %in% colnames(specific_points_MAGE))) {
        stop("specific_points_MAGE must contain columns: reach, KP_start, KP_end")
    }
    if (!all(required_cols %in% colnames(specific_points_BaM))) {
        stop("specific_points_BaM must contain columns: reach, KP_start, KP_end")
    }

    structure(
        list(
            MAGE = list(specific_points = specific_points_MAGE),
            BaM = list(
                specific_points = specific_points_BaM,
                spatialisation_function = spatialisation_function,
                n_degree = n_degree,
                prior = prior
            )
        ),
        class = "su_structure"
    )
}



ID_Reaches_SU <- list(
    "1" = list(
        Kmin = list(
            SU_1 = list(
                MAGE = list(
                    specific_points = data.frame(1, 2, 3)
                ),
                BaM = list(
                    specific_points = data.frame(1, 2, 3),
                    spatialisation_function = custom_function(),
                    n_degree = 0,
                    prior = list(
                        RBaM::parameter(
                            name = "a0_min_glass",
                            init = 1 / 0.010, # Initial guess
                            prior.dist = "FlatPrior+", # Prior distribution
                            prior.par = NULL
                        ),
                        RBaM::parameter(
                            name = "a0_flood",
                            init = (33 + 1 / 0.013) / 2,
                            prior.dist = "FlatPrior+",
                            prior.par = NULL
                        )
                    )
                )
            )
        ),
        Kmoy =
            list(
                SU_1 = list(
                    MAGE = list(
                        specific_points = data.frame(1, 2, 3)
                    ),
                    BaM = list(
                        specific_points = data.frame(1, 2, 3),
                        spatialisation_function = custom_function(),
                        n_degree = 0,
                        prior = list(
                            RBaM::parameter(
                                name = "a0_min_glass",
                                init = 1 / 0.010, # Initial guess
                                prior.dist = "FlatPrior+", # Prior distribution
                                prior.par = NULL
                            ),
                            RBaM::parameter(
                                name = "a0_flood",
                                init = (33 + 1 / 0.013) / 2,
                                prior.dist = "FlatPrior+",
                                prior.par = NULL
                            )
                        )
                    )
                )
            )
    )
)



Kmin <- list(
    Main_SU_1 = list(
        BaM_specific_points = data.frame(
            reach = c(1),
            KP_start = c(0),
            KP_end = c(18)
        ),
        MAGE_specific_points = data.frame(
            reach = c(1),
            KP_start = c(0),
            KP_end = c(18)
        )
    ),
    spatialisation_function = "constante_piecewise_function()"
)

metadata_Kflood <- list(
    Main_SU_1 = list(
        BaM_specific_points = data.frame(
            reach = c(1, 2, 3),
            KP_start = c(0, 9.05, 15),
            KP_end = c(9.05, 15, 18)
        ),
        MAGE_specific_points = data.frame(
            reach = c(1),
            KP_start = c(0),
            KP_end = c(18)
        )
    ),
    spatialisation_function = "constante_piecewise_function()"
)


all_SUs <- list(
    SU_1 = SU_1
)

if (number_of_SU != length(all_SUs)) stop("Number of SU are not the same at all_SUs declared")

# CHECKS:
for (i in 1:length(all_SUs)) {
    if (length(all_SUs[[i]]$Kmin) != length(all_SUs[[i]]$Kflood)) {
        stop(paste0("The number of elements in Kmin and Kflood must be the same for the SU named ", names(all_SUs)[[i]]))
    }

    check_kp_first_reach(
        MAGE_specific_points = all_SUs[[i]]$MAGE$spatial_specific_points,
        Kmin_specific_ponits = all_SUs[[i]]$Kmin$spatial_specific_points,
        Kflood_specific_ponits = all_SUs[[i]]$Kflood$spatial_specific_points
    )
    check_metadata(metadata = all_SUs[[i]])

    check_prior(metadata = all_SUs[[i]]$Kmin)
    check_prior(metadata = all_SUs[[i]]$Kflood)
}
############################################
# End metadata creation
############################################

############################################
# Set calibration error model setting
############################################
remant_error_list <- list(
    remnantErrorModel(
        fname = "Config_RemnantSigma.txt",
        funk = "Constant",
        par = list(parameter(
            name = "intercept",
            init = 0.0005,
            prior.dist = "FlatPrior",
        ))
    ),
    remnantErrorModel(
        fname = "Config_RemnantSigma2.txt",
        funk = "Constant",
        par = list(parameter(
            name = "intercept",
            init = 10,
            prior.dist = "LogNormal",
            prior.par = c(log(10), 0.2)
        ))
    ),
    remnantErrorModel(
        fname = "Config_RemnantSigma3.txt",
        funk = "Constant",
        par = list(parameter(
            name = "intercept",
            init = 10,
            prior.dist = "LogNormal",
            prior.par = c(log(10), 0.2)
        ))
    ),
    remnantErrorModel(
        fname = "Config_RemnantSigma4.txt",
        funk = "Constant",
        par = list(parameter(
            name = "intercept",
            init = 15,
            prior.dist = "LogNormal",
            prior.par = c(log(15), 0.2)
        ))
    ),
    remnantErrorModel(
        fname = "Config_RemnantSigma5.txt",
        funk = "Constant",
        par = list(parameter(
            name = "intercept",
            init = 15,
            prior.dist = "LogNormal",
            prior.par = c(log(15), 0.2)
        ))
    )
)

############################################
# Grid
############################################
# Remarks: A common RUG file is necessary to perform the calculation, it does not matter if there is minor bed or floodplain
# Spatialisation file could be different between kmin and kflood, respecting only the condition that number of rows must be n + 1 number of rows in RUG file
# Number of covariate depends on the number of reaches in BaM specifications


for (su in all_SUs) {
    # Import all points from the same SU in the main channel, floodplain and MAGE model
    all_points <-
        c(
            as.vector(unlist(su$MAGE_specific_points[-1])), # Remove the reach, get only KP start and end
            as.vector(unlist(su$Kmin$spatial_specific_points[-1])),
            as.vector(unlist(su$Kflood$spatial_specific_points[-1]))
        )

    all_points <- all_points[-which(duplicated(all_points))]
    all_points <- sort(all_points)

    # Calculate the interpolated grid using specific points from MAGE, Kmin and Kflood
    grid_covariate_streamwise <-
        sort(c(
            interpolation_specific_points(
                total_points = metadata_settings$number_points_interpolation,
                all_specific_points = all_points
            ),
            all_points[-c(1, length(all_points))] # Remove boundaries values, only values between must be added
        ))

    # Assign the reach number according to the calculated grid
    su$Kflood$reaches <- assign_reach_from_a_grid(
        metadata_specific_points = su$Kflood$spatial_specific_points,
        grid = grid_covariate_streamwise
    )
    su$Kmin$reaches <- assign_reach_from_a_grid(
        metadata_specific_points = su$Kmin$spatial_specific_points,
        grid = grid_covariate_streamwise
    )

    su$MAGE$reaches <- assign_reach_from_a_grid(
        metadata_specific_points = su$MAGE$spatial_specific_points,
        grid = grid_covariate_streamwise
    )

    #
}

############################################
# End grid
############################################


############################################
# General settings
############################################
## Common for all run of the case of study.
# Harcode values, it does not matter, it just useful for give the size of the .RUG File
RUG_Kmin <- rep(20, length(su$MAGE$reaches))
RUG_Kflood <- rep(10, length(su$MAGE$reaches))

# Setting jump standard deviation for MCMC sampling
jump_MCMC_theta_param_user <- 8
jump_MCMC_error_model_user <- 0.001
threshold_jump_MCMC_error_model <- 0.5
prior_error_model <- get_init_prior(remant_error_list)

mod_polynomials <- list()

############################################
# End general settings
############################################

############################################
# Start calibration
############################################










### To be added
Experiment_id <- all_experiments

path_experiment <- file.path(file_main_path, Experiment_id)
counter_model <- 1

ggsave(file.path(path_experiment, "CalData_plot.png"),
    plot = CalData_plot_export,
    width = 20,
    height = 20,
    units = "cm"
)
