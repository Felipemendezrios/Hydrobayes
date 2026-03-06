# Function to get covariate in piecewise function
getCovariate_piecewise <- function(KP_grid, shiftPoints) {
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

    if (!(all(diff(KP_grid) > 0) || all(diff(KP_grid) < 0))) {
        return(stop("KP_grid must be strictly increasing or strictly decreasing"))
    }

    increasing <- ifelse(all(diff(KP_grid) > 0), TRUE, FALSE)

    all_shiftPoints <- c(
        shiftPoints,
        KP_grid[1],
        last(KP_grid)
    ) %>%
        sort(decreasing = !increasing)

    # Assign each value in KP_grid to an interval
    intervals <- cut(KP_grid, breaks = all_shiftPoints, include.lowest = TRUE, right = TRUE)
    # Keep the same order
    intervals <- factor(intervals, levels = unique(intervals))

    # Create binary matrix:
    covariate_matrix <- as.matrix(as.data.frame(stats::model.matrix(~ intervals - 1)))

    colnames(covariate_matrix) <- NULL
    return(covariate_matrix)
}

# Function to get covariate in Legendre polynomial
getCovariate_Legendre <- function(max_polynomial_degree, ...) {
    if (!(max_polynomial_degree >= 0 && max_polynomial_degree == floor(max_polynomial_degree))) stop("max_polynomial_degree must be a positive integer")

    # Capture the first argument passed via ...
    args <- list(...)
    # Check if required arguments are present
    if (!("covariate_discretization" %in% names(args))) {
        stop("covariate_discretization must be provided as a named argument")
    }
    # Extract covariate_discretization
    covariate_discretization <- args$covariate_discretization

    Z <- matrix(
        data = 0,
        nrow = length(covariate_discretization),
        ncol = max_polynomial_degree + 1 # Because degrees start at 0
    )
    for (polynomial_degree in seq(0, max_polynomial_degree)) {
        Z[, polynomial_degree + 1] <- getlegendre(
            polynomial_degree = polynomial_degree,
            covariate_discretization = covariate_discretization
        )
    }

    return(Z)
}
