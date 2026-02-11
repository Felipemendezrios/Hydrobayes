# Function to estimate covariate's values of a polynomial degree fixed
getlegendre <- function(
    polynomial_degree,
    covariate_discretization) {
    if (length(unique(covariate_discretization)) == 1) {
        stop("covariate_discretization must contain more than one unique value")
    }
    normalized_values <- 2 * (covariate_discretization - min(covariate_discretization)) / (max(covariate_discretization) - min(covariate_discretization)) - 1

    # Write co variant matrix for main channel
    # Create an empty data frame and fill it with polynomial values
    legendre_df_covariate <- data.frame(x = normalized_values) # Start with x values
    legendre_df <- legendre_df_covariate

    if (polynomial_degree < 0) stop("The polynomial degree must be positive")
    if (any(!dplyr::between(normalized_values, -1, 1))) stop("Range of normalized_values should be between [-1,1]")

    if (polynomial_degree == 0) {
        return(rep(1, length(normalized_values)))
    }
    if (polynomial_degree == 1) {
        return(normalized_values)
    }

    # P_0(normalized_values) = 1
    Pn_1 <- rep(1, length(normalized_values))
    # P_1(normalized_values) = normalized_values
    Pn <- normalized_values

    for (k in 1:(polynomial_degree - 1)) {
        P_next <- ((2 * k + 1) * normalized_values * Pn - k * Pn_1) / (k + 1)
        Pn_1 <- Pn
        Pn <- P_next
    }
    return(Pn)
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
