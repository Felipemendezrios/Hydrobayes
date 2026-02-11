# Function to get covariate in piecewise function
getCovariate_piecewise <- function(KP_grid, shiftPoints) {
    if (any(diff(KP_grid) < 0)) {
        return(stop("KP grid must be increasing all the time"))
    }
    if (KP_grid[1] == min(shiftPoints)) {
        return(stop("Shift point must be different to the first KP_grid"))
    }
    if (last(KP_grid) == max(shiftPoints)) {
        return(stop("Shift point must be different to the last KP_grid"))
    }
    if (!any(KP_grid < min(shiftPoints))) {
        return(stop(paste0("Shift point ", min(shiftPoints), "is out of the KP grid")))
    }
    if (!any(KP_grid > max(shiftPoints))) {
        return(stop(paste0("Shift point ", max(shiftPoints), " is out of the KP grid")))
    }
    if (any(duplicated(shiftPoints))) {
        return(stop("Non duplicated values in shiftPoints"))
    }

    all_shiftPoints <- sort(c(
        shiftPoints,
        KP_grid[1],
        last(KP_grid)
    ))

    # Assign each value in KP_grid to an interval
    intervals <- cut(KP_grid, breaks = all_shiftPoints, include.lowest = TRUE, right = TRUE)

    # Create binary matrix version 0:
    covariate_matrix_v0 <- stats::model.matrix(~ intervals - 1)
    covariate_matrix_v1 <- cbind(KP_grid, covariate_matrix_v0)
    # Get duplicated values from grid respresenting switch of reach
    id_shiftSReaches <- which(duplicated(KP_grid) | duplicated(KP_grid, fromLast = TRUE))

    if (length(id_shiftSReaches) != 0) {
        # Identify which of these are the second (or later) occurrences
        is_second_occurrence <- duplicated(KP_grid[id_shiftSReaches])
        # If TRUE, matrix should be corrected to consider shift point coincide with switch of reach
        # Correct the matrix for the second (or later) occurrences
        for (i in which(is_second_occurrence)) {
            id <- id_shiftSReaches[i]
            # Find the column where the current 1 is set
            current_col <- which(covariate_matrix_v1[id, ] == 1)
            # If the current column is not the last one, set the next column to 1 and the current to 0
            if (current_col < ncol(covariate_matrix_v1)) {
                covariate_matrix_v1[id, current_col] <- 0
                covariate_matrix_v1[id, current_col + 1] <- 1
            }
        }
    }
    covariate_matrix <- covariate_matrix_v1[, -c(1)]
    colnames(covariate_matrix) <- NULL
    return(covariate_matrix)
}
