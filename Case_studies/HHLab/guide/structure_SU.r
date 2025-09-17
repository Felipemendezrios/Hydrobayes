rm(list = ls())
graphics.off()

# Set directory
dir_workspace <- here::here()



#############################################
# Information obtained from hydraulic model
#############################################
# Reaches and KP

Input_hydraulic_model <- data.frame(
    reach = c(1, 2, 3, 4),
    KP_start = c(0, 5, 10, 25),
    KP_end = c(10, 15, 20, 35)
)


#############################################
# Functions
#############################################

# Interpolation function passing by some specific nodes:
interpolation_specific_points <- function(total_points = 100,
                                          specific_nodes) {
    if (total_points <= length(specific_nodes)) stop("Total points defined is lower or equal to the specific points required. Please either increase the number of total_points or not used interpolation function")
    # Calculate the intervals between the specific points
    intervals <- diff(specific_nodes)
    # Distribute the number of points proportionally to each interval
    total_local_points <- total_points - 1 # last value during diff function
    points_per_interval_float <- (intervals / sum(intervals)) * total_local_points

    points_per_interval <- round(points_per_interval_float)
    # Due to round function, we need to check the size
    # Check if more points than expected, need to remove some of them
    if (sum(points_per_interval) > total_local_points) {
        id_order <- order(points_per_interval_float %% 1)
        i <- 1
        while (sum(points_per_interval) > total_local_points) {
            id_position <- which(i == id_order)
            points_per_interval[id_position] <- points_per_interval[id_position] - 1
            i <- i + 1
        }
        # Check if less points than expected, need to add some of them
    } else if (sum(points_per_interval) < total_local_points) {
        id_order <- order(points_per_interval_float %% 1)
        i <- 1
        while (sum(points_per_interval) < total_local_points) {
            id_position <- which(i == id_order)
            points_per_interval[id_position] <- points_per_interval[id_position] + 1
            i <- i + 1
        }
    }

    # Generate sub-sequences for each interval, excluding the duplicate endpoint
    specific_seq_temp <- unlist(mapply(
        function(start, end, n) seq(start, end, length.out = n + 1)[-(n + 1)], # Exclude the endpoint
        specific_nodes[-length(specific_nodes)], # Start of each interval
        specific_nodes[-1], # End of each interval
        points_per_interval # Number of points per interval
    ))

    # Add the final endpoint manually
    specific_seq <- c(specific_seq_temp, specific_nodes[length(specific_nodes)])

    if (!all(specific_nodes %in% specific_seq)) stop("Distance so close, increase number of interpolate KP to get all points required")

    if (length(specific_seq) != total_points) stop("Something is wrong with the interpolation")
    return(specific_seq)
}

# Function to assign reach to a grid
assign_reach_from_a_grid <- function(reach_KP_boundaries, grid) {
    reaches <- rep(NA, length(grid))

    if (nrow(reach_KP_boundaries) == 1) {
        all_set <- dplyr::between(grid, reach_KP_boundaries$KP_start, reach_KP_boundaries$KP_end)
        if (!all(all_set)) stop("reach_KP_boundaries must contain all grid")
        reaches[which(all_set)] <- reach_KP_boundaries$reach[1]
    } else {
        duplicated_numbers <- grid[duplicated(grid)]
        if ((nrow(reach_KP_boundaries) - 1) != length(duplicated_numbers)) stop("Number of rows of specific_points_Model minus 1 must coincide with the number of repeated values in the grid")

        for (i in 1:nrow(reach_KP_boundaries)) {
            idx_position_reaches <- which(dplyr::between(grid, reach_KP_boundaries$KP_start[i], reach_KP_boundaries$KP_end[i]))

            idx_duplicated_values <- which(grid[idx_position_reaches] %in% duplicated_numbers)

            if (length(idx_duplicated_values) != 0) {
                if (length(idx_duplicated_values) != 2 & length(idx_duplicated_values) != 4) stop("Duplicated values are not found, the number of duplicated values must be two or four")

                # First duplicated values
                if (i == 1 & length(idx_duplicated_values) == 2) {
                    reaches[idx_position_reaches] <- reach_KP_boundaries$reach[i]

                    reaches[idx_position_reaches[idx_duplicated_values[2]]] <- reach_KP_boundaries$reach[i + 1]
                } else if (length(idx_duplicated_values) == 4) {
                    idx_position_reaches <- idx_position_reaches[-idx_duplicated_values[c(1, 2)]]

                    reaches[idx_position_reaches] <- reach_KP_boundaries$reach[i]

                    reaches[idx_position_reaches[(idx_duplicated_values[4] - 2)]] <- reach_KP_boundaries$reach[i + 1] # -2 because two values have been removed previously
                } else if (length(idx_duplicated_values) == 2) {
                    idx_position_reaches <- idx_position_reaches[-idx_duplicated_values[c(1, 2)]]
                    reaches[idx_position_reaches] <- reach_KP_boundaries$reach[i]
                }
            } else {
                stop("It should never get here. Something is wrong")
                reaches[idx_position_reaches] <- reach_KP_boundaries$reach[i]
            }
        }
    }

    if (any(is.na(reaches))) {
        return(stop("Some values in the vector reaches have not been assigned"))
    }

    return(reaches)
}
assign_Flat_Prior <- function(prior) {
    assign_Flat_properties <- list()
    i <- 1
    for (id in prior) {
        assign_Flat_properties[[i]] <- RBaM::parameter(
            name = id$name,
            init = id$init
        )
        i <- i + 1
    }
    return(assign_Flat_properties)
}
#############################################
# Step 1: SU constructor
#############################################
# Information needed to build a SU :
# Key between MAGE reaches and SU
Input_SU_key_reachesMage <- list(
    c(1, 3),
    2,
    4
)

# Declare how many covariate do you have and give the list of these covariates?
Input_n_covariate <- 2
# piecewise function or Legendre (polynomial)

# Specify prior by SU
Input_SU_prior <- list(
    SU_1 = list(
        data.frame(
            name = "a0_Kmin_wood",
            init = 1 / 0.010
        ),
        data.frame(
            name = "a0_Kmin_grass",
            init = 1 / 0.010
        )
    )
)

SU <- list()
for (i in seq_along(Input_SU_key_reachesMage)) {
    # 1. KP grid
    hydraulic_model_SU_i <- Input_hydraulic_model[which(Input_hydraulic_model$reach %in% Input_SU_key_reachesMage[[i]]), ]

    nodes <- c(
        hydraulic_model_SU_i$KP_start,
        hydraulic_model_SU_i$KP_end
    )

    nodes_to_interpolate <- unique(nodes)

    # I need to duplicated the node to keep in mind RUG specifications
    SU_KP <- sort(
        c(
            interpolation_specific_points(
                total_points = 100,
                specific_nodes = nodes_to_interpolate
            ),
            nodes_to_interpolate[-c(1, length(nodes_to_interpolate))] # Remove intermediate nodes of the interpolated nodes
        )
    )
    # 3. ID reach MAGE
    SU_reach <- assign_reach_from_a_grid(
        reach_KP_boundaries = hydraulic_model_SU_i,
        grid = SU_KP
    )

    if (length(SU_reach) != length(SU_KP)) stop("Lengths of SU_reach and SU_KP must be equal")

    # 3. Covariate matrix
    covariate_list <- list(
        SU_KP
    )

    if (length(covariate_list) != Input_n_covariate) stop("List of the covariates must equal to the number of covariates declared")
    if (length(covariate_list) > 1) {
        # Get lengths of each component
        lengths_list <- lengths(covariate_list)

        # Check if all lengths are equal
        if (!all(lengths_list == lengths_list[1])) stop("Length of the component of the covariates list must be equal, same discretization is mandatory")
    }

    SU_covariate_matrix <- matrix(
        nrow = length(covariate_list[[1]]),
        ncol = Input_n_covariate
    )

    # Fill covariate matrix
    for (i in Input_n_covariate) {
        SU_covariate_matrix[, i] <- covariate_list[[i]]
    }

    # 4. Prior
    SU_Prior <- assign_Flat_Prior(prior = Input_SU_prior[[i]])

    # Final step: build SU

    SU[[i]] <- list(
        KP = SU_KP,
        reach = SU_reach,
        covariate_Matrix = SU_covariate_matrix,
        prior = SU_Prior
    )
}

#############################################
# End SU constructor
#############################################

#############################################
# Covariantes utilities
#############################################

#############################################
# Constant piecewise function
#############################################

# constant_piecewise_function <- function(

# )
