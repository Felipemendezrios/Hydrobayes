# Interpolation function passing by some specific nodes:
interpolation_specific_points <- function(total_points = 100,
                                          specific_nodes) {
    if (total_points <= length(specific_nodes)) stop("Total points defined is lower or equal to the specific points required. Please either increase the number of total_points or avoid to used interpolation function")

    # Order specific_nodes
    specific_nodes <- sort(specific_nodes)

    # Calculate the intervals between the specific points
    intervals <- diff(specific_nodes)
    # Distribute the number of points proportionally to each interval
    total_local_points <- total_points - 1 # last value during diff function
    points_per_interval_float <- (intervals / sum(intervals)) * total_local_points

    points_per_interval <- round(points_per_interval_float)

    # Generate sub-sequences for each interval, excluding the duplicate endpoint
    specific_seq_temp <- unlist(mapply(
        function(start, end, n) seq(start, end, length.out = n + 1)[-(n + 1)], # Exclude the endpoint
        specific_nodes[-length(specific_nodes)], # Start of each interval
        specific_nodes[-1], # End of each interval
        points_per_interval # Number of points per interval
    ))

    # Add the final endpoint manually
    specific_seq <-
        c(
            specific_seq_temp,
            last(specific_nodes)
        )


    if (!all(specific_nodes %in% specific_seq)) stop("Distance so close, increase number of interpolate KP to get all points required")

    return(specific_seq)
}


# Function to assign reach to a grid
assign_reach_from_a_grid <- function(reach_KP_boundaries, grid, Logical_decreasing) {
    grid_with_reaches <- data.frame(
        reaches = rep(NA, length(grid)),
        grid = rep(NA, length(grid))
    )

    if (nrow(reach_KP_boundaries) == 1) {
        all_set <- dplyr::between(grid, reach_KP_boundaries$KP_start, reach_KP_boundaries$KP_end)
        if (!all(all_set)) stop("reach_KP_boundaries must contain all grid")
        grid_with_reaches$reaches[which(all_set)] <- reach_KP_boundaries$reach[1]
        grid_with_reaches$grid[which(all_set)] <- grid
    } else {
        if (Logical_decreasing) {
            # Include end (right) and exclude start (left)
            for (i in 1:nrow(reach_KP_boundaries)) {
                if (i != nrow(reach_KP_boundaries)) {
                    idx_position_reaches <- which(
                        grid > reach_KP_boundaries$KP_start[i] &
                            grid <= reach_KP_boundaries$KP_end[i]
                    )
                } else {
                    idx_position_reaches <- which(
                        grid >= reach_KP_boundaries$KP_start[i] &
                            grid <= reach_KP_boundaries$KP_end[i]
                    )
                }

                grid_with_reaches$reaches[idx_position_reaches] <- reach_KP_boundaries$reach[i]
                grid_with_reaches$grid[idx_position_reaches] <- grid[idx_position_reaches]
            }
        } else {
            # Include start (left) and exclude end (right)
            for (i in 1:nrow(reach_KP_boundaries)) {
                if (i != nrow(reach_KP_boundaries)) {
                    idx_position_reaches <- which(
                        grid >= reach_KP_boundaries$KP_start[i] &
                            grid < reach_KP_boundaries$KP_end[i]
                    )
                } else {
                    idx_position_reaches <- which(
                        grid >= reach_KP_boundaries$KP_start[i] &
                            grid <= reach_KP_boundaries$KP_end[i]
                    )
                }


                grid_with_reaches$reaches[idx_position_reaches] <- reach_KP_boundaries$reach[i]
                grid_with_reaches$grid[idx_position_reaches] <- grid[idx_position_reaches]
            }
        }
    }

    if (any(is.na(grid_with_reaches))) {
        return(stop("Some values in the vector reaches have not been assigned"))
    }
    return(grid_with_reaches)
}

# Block-diagonal binding of matrices
# source: https://stackoverflow.com/questions/17495841/block-diagonal-binding-of-matrices
block_diagonal_matrix <- function(...) {
    d <- list(...)
    nrows <- sum(sapply(d, NROW))
    ncols <- sum(sapply(d, NCOL))
    ans <- matrix(0, nrows, ncols)
    i1 <- 1
    j1 <- 1
    for (m in d) {
        if (!is.matrix(m)) {
            return(stop("All arguments must be matrices"))
        }
        i2 <- i1 + NROW(m) - 1
        j2 <- j1 + NCOL(m) - 1
        ans[i1:i2, j1:j2] <- m
        i1 <- i2 + 1
        j1 <- j2 + 1
    }
    return(ans)
}
