rm(list = ls())
graphics.off()

# Set directory
dir_workspace <- here::here()

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

# Common for observation and calibration
all_experiments <- c(
    "piecewise_function"
)
# Folder related to the polynomial degree calibration
if (all_experiments == "piecewise_function") {
    all_polynomial_degree <- c(
        "n_min_0_n_flood_0"
    )
} else {

}


# Folder related to the observations (careful with the order!)
all_events <- c(
    "CWMQ18",
    "CWMQ12"
)

dir_exe_BaM <- "/home/famendezrios/Documents/Git/BaM/makefile/"
MAGE_executable <- "/home/famendezrios/Documents/Softwares/pamhyr2/mage8/mage"
file_main_path <- "/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/Case_studies/HHLab/Compound_channel_spatial_friction/Calibration_experiments"
do_calibration <- FALSE

# Observations data input:
# Measurements for calibration by event!
# Ensure that order should be in coherence with mage model set up!
Cal_measures <-
    list(
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

mage_projet_name <- "CWMQ"
WSE_data_temp <- WSE_data_temp[str_detect(WSE_data_temp$ID_experiment, pattern = mage_projet_name), ]

# This part help to handle when observed data is presented as a list and user need to specify the the order of the events with the same name for matching with all_all_events, where the order is important

# Only giving the name keeping the data order and creating a link with the names of the event and the Cal_measures defined previously
match_case_event <- data.frame(
    idx_data = unique(levels(factor(WSE_data_temp$ID_experiment))),
    idx_event = unique(levels(factor(WSE_data_temp$ID_experiment)))
)

Link_x_t_ind_event <- lapply(names(Cal_measures), function(event) {
    # Find the row in match_case_event where idx_event matches the current event
    row_idx <- which(match_case_event$idx_event == event)

    list(
        idx_data = match_case_event$idx_data[row_idx],
        idx_event = event,
        Cal_measure = Cal_measures[[event]]
    )
})

# # Initialize an empty list to store the results
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

# Assign manually uncertainty values in meters
do_manual_uncertainty <- TRUE
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
########


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

###################
# Create the metadata to link reaches from MAGE to BaM! and inversely
###################
# Remark: reaches could be equal in BaM! and MAGE, but it is not mandatory

# Two cases are supported:

# 1. In case of spatially distributed friction: BaM! integrates several reaches from MAGE to keep a same spatialisation avoiding discontinuities at nodes of two reaches in MAGE model. So the number reaches in BaM! could be lower than MAGE. The extrem case is when the number of reaches in BaM! are the same from MAGE. A different set of reaches could be done to the main channel and to the floodplain

# E.g: Here the reach in BaM number one (1) in the main channel contains 2 consecutive reaches (1,2) from MAGE perspective
# toto <- list(
#     Kmin = list(
#         link_1 = list(reach_BaM = 1, reach_MAGE = c(1, 2)),
#         link_2 = list(reach_BaM = 2, reach_MAGE = c(3, 5, 6)),
#         link_3 = list(reach_BaM = 3, reach_MAGE = c(4))
#     ),
#     Kflood = list(
#         link_1 = list(reach_BaM = 1, reach_MAGE = c(1, 2, 3, 5)),
#         link_2 = list(reach_BaM = 2, reach_MAGE = c(4, 6))
#     )
# )
# # Validate uniqueness of reach_MAGE
# all_reaches_MAGE <- unlist(lapply(toto, function(x) x$reach_MAGE))
# if (length(all_reaches_MAGE) != length(unique(all_reaches_MAGE))) {
#  stop("Duplicate reach_MAGE values detected!")
# }
# all_reaches_BaM <- unlist(lapply(toto, function(x) x$reach_BaM))
# if (length(all_reaches_BaM) != length(unique(all_reaches_BaM))) {
#   stop("Duplicate reach_BaM values detected!")
# }



# 2. In case of piecewise function: the number of the reaches in MAGE is lower than BaM!. Inversely at the previous case, in a reach of Mage, it could be more than one reach from BaM! perspective. This method let us estimated a piecewise function in a reach of MAGE.

# E.g: Here the reach in MAGE number one (1) in the main channel contains 2 consecutive reaches (1,2) from BaM perspective to estimate a piecewise function
#  list(
#         Kmin = list(
#             link_1 = list(reach_BaM = c(1, 2), reach_MAGE = 1),
#             link_2 = list(reach_BaM = c(3, 5, 6), reach_MAGE = 2),
#             link_3 = list(reach_BaM = c(4), reach_MAGE = 3)
#         ),
#         Kflood = list(
#             link_1 = list(reach_BaM = c(1, 2, 3, 5), reach_MAGE = 1),
#             link_2 = list(reach_BaM = c(4, 6), reach_MAGE = 2)
#         )
#     )



# Both of them are implemented here:
if (all_experiments == "piecewise_function") {
    metadata <- list(
        Kmin = list(
            link_1 = list(reach_BaM = c(1), reach_MAGE = 1)
        ),
        Kflood = list(
            link_1 = list(reach_BaM = c(1, 2), reach_MAGE = 1)
        )
    )
} else {
    stop("not supported yet")
}
# Validate uniqueness of reach_MAGE
all_reaches_MAGE <- unlist(lapply(toto, function(x) x$reach_MAGE))
if (length(all_reaches_MAGE) != length(unique(all_reaches_MAGE))) {
    stop("Duplicate reach_MAGE values detected!")
}
all_reaches_BaM <- unlist(lapply(toto, function(x) x$reach_BaM))
if (length(all_reaches_BaM) != length(unique(all_reaches_BaM))) {
    stop("Duplicate reach_BaM values detected!")
}

# Give the reaches connected
# Reaches for spatially distributed the friction : the order here will be important (be careful depending on the mage model set up). That will impact the way that positions are interpretated.

# Number of reaches in MAGE
Number_reaches_Mage_Model <- list(1)

## Definition of the number of reaches must be adapted to the larger number of segment during calibration either in tha main channel or in the floodplain. A value must be done by number of reaches in MAGE
Number_segment_Cal_Kmin <- list(c(1))
Number_segment_Cal_Kflood <- list(c(1, 2))

if (length(Number_reaches_Mage_Model) != length(Number_segment_Cal_Kmin) | length(Number_reaches_Mage_Model) != length(Number_segment_Cal_Kflood) | length(Number_segment_Cal_Kmin) != length(Number_segment_Cal_Kflood)
) {
    stop("Size of Number_reaches_Mage_Model, Number_segment_Cal_Kmin and Number_segment_Cal_Kflood must be equal")
}

# # Example data
# Number_reaches_Mage_Model <- list(1, 2)  # IDs of the reaches
# Number_segment_Cal_Kmin <- list(c(1), c(1,2,3,4))  # Segments for main channel
# Number_segment_Cal_Kflood <- list(c(1,2), c(1))  # Segments for floodplain


# Mage needs the maximal number of division between the main channel and the floodplain to create the .RUG file. But the estimation will be performed separately.
Maximal_Number_segment_Cal_MAGE <- ifelse(Number_segment_Cal_Kmin >= Number_segment_Cal_Kflood, "Kmin", "Kflood")

# case <- "independant"
# cases:
# follow: two reaches are consecutive
# independant: in a single reach, a piecewise function will be estimated either in the main channel or the floodplain
#  monoreach: a single reach

# if (case == "independant") {
#     reaches_user_Kmin <- list(
#         1
#     )
#     reaches_user_Kflood <- list(
#         1,
#         2
#     )
# } else if (case == "follow") {
#     reaches_user_Kmin <- list(
#         c(1)
#     )
#     reaches_user_Kflood <- list(
#         c(1, 2)
#     )
# } else if (case == "monoreach") {
#     reaches_user_Kmin <- list(
#         1
#     )
#     reaches_user_Kflood <- list(
#         1
#     )
# } else {
#     stop("ERROR")
# }

if (length(reaches_user_Kmin) != Number_segment_Cal_Kmin) stop("reaches_user_Kmin must be coherent to the number of segment during calibration for the main channel")
if (length(reaches_user_Kflood) != Number_segment_Cal_Kflood) stop("reaches_user_Kflood must be coherent to the number of segment during calibration for the floodplain")

# Number of different reaches
Nb_reaches_estimation_Kmin <- length(reaches_user_Kmin)
Nb_reaches_estimation_Kflood <- length(reaches_user_Kflood)



## Calibration setting

# Input parameters in main channel:
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


## Cross-section treatment: BaM perspective
### interpolation

specific_points_Model_Kmin <- data.frame(
    reach = c(1),
    KP_start = c(0),
    KP_end = c(18)
)

specific_points_Model_Kflood <- data.frame(
    reach = c(1, 2),
    KP_start = c(0, 9.05),
    KP_end = c(9.05, 18)
)


if (Nb_reaches_estimation_Kmin != nrow(specific_points_Model_Kmin)) stop("The number of segment to assign during calibration (Nb_reaches_estimation_Kmin) must be equal to the number of reaches in specific_points_Model_Kmin")
if (Nb_reaches_estimation_Kflood != nrow(specific_points_Model_Kflood)) stop("The number of segment to assign during calibration (Nb_reaches_estimation_Kflood) must be equal to the number of reaches in specific_points_Model_Kflood")

if (nrow(specific_points_Model_Kmin) > 1) {
    for (i in 1:(nrow(specific_points_Model_Kmin) - 1)) {
        if (specific_points_Model_Kmin$KP_end[i] != specific_points_Model_Kmin$KP_start[i + 1]) stop("specific_points_Model_Kmin must respect that KP_start of reach i must be equal to the KP_end of the reach i + 1 ")
    }
}
if (nrow(specific_points_Model_Kflood) > 1) {
    for (i in 1:(nrow(specific_points_Model_Kflood) - 1)) {
        if (specific_points_Model_Kflood$KP_end[i] != specific_points_Model_Kflood$KP_start[i + 1]) stop("specific_points_Model_Kflood must respect that KP_start of reach i must be equal to the KP_end of the reach i + 1 ")
    }
}

if (specific_points_Model_Kmin$KP_start[1] != specific_points_Model_Kflood$KP_start[1]) stop("First cross-section is not the same between the main channel and the floodplain. Check specific_points_Model_Kmin and specific_points_Model_Kflood")

if (last(specific_points_Model_Kmin$KP_end) != last(specific_points_Model_Kflood$KP_end)) stop("Last cross-section is not the same between the main channel and the floodplain. Check specific_points_Model_Kmin and specific_points_Model_Kflood")


# Fit to either the size of data frames of the main channel or floodplain to pass into .RUG file
if (Maximal_Number_segment_Cal_MAGE == "Kmin") {
    specific_points_Model <- specific_points_Model_Kmin
    reaches_user <- reaches_user_Kmin
} else {
    specific_points_Model <- specific_points_Model_Kflood
    reaches_user <- reaches_user_Kflood
}

# Real reaches from Mage projet (interpolation by reach)
grid_covariant_discretized_by_reach <- list()
for (i in 1:nrow(specific_points_Model)) {
    boundaries_points <- c(specific_points_Model$KP_start[i], specific_points_Model$KP_end[i])

    grid_covariant_discretized_by_reach[[i]] <- data.frame(
        real_reach = specific_points_Model$reach[i],
        interpol_values = interpolation_specific_points(
            total_points = 100,
            all_specific_points = boundaries_points
        )
    )
}

grid_covariant_discretized_real_bief <- bind_rows(grid_covariant_discretized_by_reach)

# position_direction
# increasing: positions are increasing in flow direction
# decreasing: positions are decreasing in flow direction
position_direction <- "increasing"
# For example: position_direction = increasing and reaches_user_Kmin is c(1,2)

# Only a unit will be considered, composed of two biefs, one and two. A check will be performed to check if all positions are increasing as defined in input data

# Second example: position_direction = decreasing and reaches_user_Kmin is c(2,1)

# In this case, reach 2 is downstream and reach 1 is upstream. A check will be performed to be sure that position will decrease by reach


grid_covariant_discretized <- c()
# Change to reach for spatialisation
for (i in 1:length(reaches_user)) {
    mask <- which(grid_covariant_discretized_real_bief$real_reach %in% reaches_user[[i]])

    if (length(mask) == 0) stop(paste0("Real reaches are ", specific_points_Model$reach, " and the reaches defined by the user are not matching: ", reaches_user))

    df_covariant_by_reach <- grid_covariant_discretized_real_bief[mask, ]

    if (length(reaches_user[[i]]) != 1) {
        # Case union of real bief for spatialisation
        if (position_direction == "increasing") {
            # Check if all interpol_values are positive for each real_reach
            check <- df_covariant_by_reach %>%
                group_by(real_reach) %>%
                summarize(check = all(diff(interpol_values) > 0, na.rm = TRUE))


            if (!any(check$check)) stop("All position values should increase in flow direction")
        } else if (position_direction == "decreasing") {
            # Check if all interpol_values are positive for each real_reach
            check <- df_covariant_by_reach %>%
                group_by(real_reach) %>%
                summarize(check = all(diff(interpol_values) < 0, na.rm = TRUE))

            if (!all(check$check)) stop(paste0("All position values should decrease in flow direction. Please check reach(es) ", check$real_reach[check$check == FALSE]))
        } else {
            stop(paste0("The position_direction variable should be either increasing or decreasing in flow direction all the time"))
        }
    }
    #### Save data independently of the case
    if (i == 1) {
        grid_covariant_discretized <- data.frame(
            Reach = i,
            Covariate = df_covariant_by_reach$interpol_values
        )
    } else {
        grid_covariant_discretized_temp <- data.frame(
            Reach = i,
            Covariate = df_covariant_by_reach$interpol_values
        )
        grid_covariant_discretized <- rbind(
            grid_covariant_discretized, grid_covariant_discretized_temp
        )
    }
}

# Adapt the values of the reach separately from the main channel to the floodplain, but keeping the same covariate discretization
if (Maximal_Number_segment_Cal_MAGE == "Kmin") {
    grid_covariant_discretized_Kmin <- grid_covariant_discretized
    grid_covariant_discretized_Kflood <- assign_reach(
        specific_points_Model = specific_points_Model_Kflood,
        grid_covariant_discretized$Covariate
    )
} else {
    grid_covariant_discretized_Kflood <- grid_covariant_discretized
    grid_covariant_discretized_Kmin <- assign_reach(
        specific_points_Model = specific_points_Model_Kmin,
        covariate_value = grid_covariant_discretized$Covariate
    )
}

##################################
## Common for all run of the case of study.
# Harcode values, it does not matter, it just useful for give the size of the .RUG File
RUG_Kmin <- rep(20, nrow(specific_points_Model))
RUG_Kflood <- rep(10, nrow(specific_points_Model))

# Setting jump standard deviation for MCMC sampling
jump_MCMC_theta_param_user <- 8
jump_MCMC_error_model_user <- 0.001
threshold_jump_MCMC_error_model <- 0.5
prior_error_model <- get_init_prior(remant_error_list)
############################

mod_polynomials <- list()

# for (Experiment_id in all_experiments) {
# }

Experiment_id <- all_experiments

path_experiment <- file.path(file_main_path, Experiment_id)
counter_model <- 1

ggsave(file.path(path_experiment, "CalData_plot.png"),
    plot = CalData_plot_export,
    width = 20,
    height = 20,
    units = "cm"
)

wood_first <- mage_projet_name == "CWMQ"

for (polynomial_id in all_polynomial_degree) {
    if (Experiment_id == "piecewise_function") {
        if (polynomial_id == "n_min_0_n_flood_0") {
            n_degree_Kmin <- 0
            n_degree_Kflood <- 0

            Kmin_prior <- list(
                a0_min_reach_1
            )

            if (wood_first) {
                Kflood_prior <- list(
                    a0_flood_reach_1_wood,
                    a0_flood_reach_1_grass
                )
            } else {
                Kflood_prior <- list(
                    a0_flood_reach_1_grass,
                    a0_flood_reach_1_wood
                )
            }
        } else {
            stop(paste0(polynomial_id, " is not supported in the experiment_id = pieceswise_function. Only case n_min_0_n_flood_0 is supported"))
        }
    } else {
        # That will work only for exploring several polynomial degrees in the main chanil, fixing the floodplain
        # if (polynomial_id == "n_0") {
        #     n_degree_Kmin <- 0
        #     Kmin_prior <- list(
        #         a0_min_reach_1
        #     )
        # } else if (polynomial_id == "n_1") {
        #     n_degree_Kmin <- 1
        #     Kmin_prior <- list(
        #         a0_min_reach_1,
        #         a1_min_reach_1
        #     )
        # } else if (polynomial_id == "n_2") {
        #     n_degree_Kmin <- 2
        #     Kmin_prior <- list(
        #         a0_min_reach_1,
        #         a1_min_reach_1,
        #         a2_min_reach_1
        #     )
        # } else if (polynomial_id == "n_3") {
        #     n_degree_Kmin <- 3
        #     Kmin_prior <- list(
        #         a0_min_reach_1,
        #         a1_min_reach_1,
        #         a2_min_reach_1,
        #         a3_min_reach_1
        #     )
        # } else {
        #     stop("polynomial_id not supported yet")
        # }
    }

    path_polynomial <- file.path(path_experiment, polynomial_id)
    # Path to save results
    workspace_user <- file.path(path_polynomial, "BaM")
    path_post_traitement <- file.path(workspace_user, "post_traitement")
    path_post_traitement_data <- file.path(path_post_traitement, "RData")

    if (!dir.exists(workspace_user)) {
        dir.create(workspace_user)
    } else {
        if (do_calibration) {
            file.remove(list.files(workspace_user, full.names = TRUE, recursive = TRUE))
        }
    }
    if (!dir.exists(path_post_traitement)) {
        dir.create(path_post_traitement)
    }

    if (!dir.exists(path_post_traitement_data)) {
        dir.create(path_post_traitement_data)
    }
    if (length(Kmin_prior) != (Nb_reaches_estimation_Kmin * max((n_degree_Kmin + 1)))) stop(paste0("More prior information (", length(Kmin_prior), ") than the number of reaches for estimation (", Nb_reaches_estimation_Kmin, ")"))
    if (length(Kflood_prior) != (Nb_reaches_estimation_Kflood * (n_degree_Kflood + 1))) stop(paste0("More prior information (", length(Kflood_prior), ") than the number of reaches for estimation (", Nb_reaches_estimation_Kflood, ")"))

    theta_param <- c(Kmin_prior, Kflood_prior)

    # Write Legendre vector:
    write.table(
        grid_covariant_discretized_Kmin,
        file = file.path(workspace_user, "legendre_covariate_non_normalized_Kmin.txt"), row.names = F
    )
    # Write Legendre vector:
    write.table(
        grid_covariant_discretized_Kflood,
        file = file.path(workspace_user, "legendre_covariate_non_normalized_Kflood.txt"), row.names = F
    )


    mageDir <- c(paste0(path_polynomial, "/model_mage/", all_events, "/"))

    RUGFiles <- paste0(mageDir, mage_projet_name, ".RUG")

    # j is a index for multi-events
    for (j in seq_along(RUGFiles)) {
        # Calculate original covariant considering that grid_covariant_discretized is in the middle of the real covariant value except on the boundaries
        grid_covariant_meshed_Model <- RUG_KP_end_by_reach <- RUG_KP_start_by_reach <- RUG_id_reach_by_reach <- RUG_Kmin_by_reach <- RUG_Kflood_by_reach <- list()

        for (reach in 1:nrow(specific_points_Model)) {
            covariant_discretized_by_reach <- grid_covariant_discretized_real_bief$interpol_values[which(grid_covariant_discretized_real_bief$real_reach == specific_points_Model$reach[reach])]

            grid_covariant_meshed_Model[[reach]] <-
                data.frame(
                    reach = specific_points_Model$reach[reach],
                    real_values = get_covariant_meshed(
                        covariant_discretized = covariant_discretized_by_reach,
                        specific_points_Model = specific_points_Model[reach, ]
                    )
                )
            RUG_id_reach_by_reach[[reach]] <- rep(specific_points_Model$reach[reach], length(covariant_discretized_by_reach))

            RUG_KP_start_by_reach[[reach]] <- grid_covariant_meshed_Model[[reach]]$real_values[-length(grid_covariant_meshed_Model[[reach]]$real_values)]

            RUG_KP_end_by_reach[[reach]] <- grid_covariant_meshed_Model[[reach]]$real_values[-1]

            RUG_Kmin_by_reach[[reach]] <- rep(RUG_Kmin[reach], length(covariant_discretized_by_reach))
            RUG_Kflood_by_reach[[reach]] <- rep(RUG_Kflood[reach], length(covariant_discretized_by_reach))
        }

        RUG_KP_start <- unlist(RUG_KP_start_by_reach)
        RUG_KP_end <- unlist(RUG_KP_end_by_reach)
        RUG_id_reach <- unlist(RUG_id_reach_by_reach)
        RUG_Kmin_args <- unlist(RUG_Kmin_by_reach)
        RUG_Kflood_args <- unlist(RUG_Kflood_by_reach)


        write_RUGFile(
            RUG_path = RUGFiles[j],
            RUG_id_reach = RUG_id_reach,
            RUG_KP_start = RUG_KP_start,
            RUG_KP_end = RUG_KP_end,
            RUG_Kmin = RUG_Kmin_args,
            RUG_Kmoy = RUG_Kflood_args,
            RUG_format = "%1s%3d      %10.0f%10.0f%10.2f%10.2f"
        )
    }

    # Spatial distributed friction in the main channel
    matrix_zFileKmin <- KFile_spatial(
        reaches = reaches_user_Kmin,
        max_degree = n_degree_Kmin,
        grid_covariant_discretized = grid_covariant_discretized_Kmin
    )

    # Spatial distributed friction in the floodplain
    matrix_zFileKmoy <- KFile_spatial(
        reaches = reaches_user_Kflood,
        max_degree = n_degree_Kflood, grid_covariant_discretized = grid_covariant_discretized_Kflood
    )

    zFileKmin <- file.path(workspace_user, "Zfile_Kmin.txt")
    zFileKmoy <- file.path(workspace_user, "Zfile_Kflood.txt")

    xtra <- xtraModelInfo(
        fname = "Config_setup.txt",
        object = list(
            exeFile = MAGE_executable,
            version = "8",
            mageDir = mageDir,
            repFile = paste0(mage_projet_name, ".REP"),
            zKmin = matrix_zFileKmin,
            zFileKmin = zFileKmin,
            doExpKmin = FALSE,
            zKmoy = matrix_zFileKmoy,
            zFileKmoy = zFileKmoy,
            doExpKmoy = FALSE
        )
    )

    mod_polynomials[[counter_model]] <- model(
        ID = "MAGE_ZQV",
        nX = 4,
        nY = 5, ,
        par = theta_param,
        xtra = xtra
    )
    mod <- mod_polynomials[[counter_model]]
    # Re-write RUG file : keep the position of the cross sections as original

    prior_theta_param <- get_init_prior(theta_param)
    prior_all_dist_theta_param <- get_init_prior(theta_param, FIX_dist = TRUE)

    data <- dataset(X = X, Y = Y, Yu = Yu, data.dir = file.path(workspace_user))

    jump_MCMC_theta_param <- ifelse(prior_theta_param != 0,
        prior_theta_param[(prior_theta_param != 0)] * 0.1,
        jump_MCMC_theta_param_user
    )

    jump_MCMC_error_model <- list()
    for (ind in 1:length(prior_error_model)) {
        jump_MCMC_error_model[[ind]] <- ifelse(prior_error_model[[ind]] > threshold_jump_MCMC_error_model,
            prior_error_model[[ind]][(prior_error_model[[ind]] != 0)] * 0.1,
            jump_MCMC_error_model_user
        )
    }

    mcmcOptions_user <- mcmcOptions(
        nCycles = 50,
        nAdapt = 50,
        manualMode = TRUE,
        thetaStd = jump_MCMC_theta_param,
        gammaStd = jump_MCMC_error_model
    )
    mcmcCooking_user <- mcmcCooking(
        burn = 0.2,
        nSlim = 1
    )
    mcmcSummary_user <- mcmcSummary(xtendedMCMC.fname = "Results_xtendedMCMC.txt")

    # Save all data used during calibration for prediction
    save(mod, data, remant_error_list,
        mcmcOptions_user, mcmcCooking_user,
        mcmcSummary_user, workspace_user,
        file = file.path(path_post_traitement_data, "BaM_objects.RData")
    )
    BaM(
        mod = mod,
        data = data,
        remnant = remant_error_list,
        mcmc = mcmcOptions_user,
        cook = mcmcCooking_user,
        summary = mcmcSummary_user,
        residuals = residualOptions(),
        pred = NULL,
        doCalib = TRUE,
        doPred = FALSE,
        na.value = -9999,
        run = FALSE,
        preClean = FALSE,
        workspace = workspace_user,
        # dir.exe = file.path(find.package("RBaM"), "bin"),
        # name.exe = "BaM",
        predMaster_fname = "Config_Pred_Master.txt"
    )
    counter_model <- counter_model + 1

    dir_cf <- file.path(workspace_user, "Config_BaM.txt")

    if (do_calibration) {
        system2(
            command = file.path(dir_exe_BaM, "BaM"),
            args = c("-cf", file.path(workspace_user, "Config_BaM.txt")),
            wait = FALSE
        )
    }
}
