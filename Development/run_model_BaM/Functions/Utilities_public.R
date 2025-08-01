#' Legendre polynomial recursive function for a only covariate
#'
#' @param degree integer, polynomial degree starting from 1
#' @param normalized_values numeric vector, the covariate normalized values between [-1,1]
#'
#' @return numeric vector, value of polynomial of degree  `n`
#'
#' @export
getlegendre <- function(degree, normalized_values) {
    if (degree < 0) stop("Degree must be positive")
    if (any(!dplyr::between(normalized_values, -1, 1))) stop("Range of normalized_values should be between [-1,1]")

    if (degree == 0) {
        return(rep(1, length(normalized_values)))
    }
    if (degree == 1) {
        return(normalized_values)
    }

    # P_0(normalized_values) = 1
    Pn_1 <- rep(1, length(normalized_values))
    # P_1(normalized_values) = normalized_values
    Pn <- normalized_values

    for (k in 1:(degree - 1)) {
        P_next <- ((2 * k + 1) * normalized_values * Pn - k * Pn_1) / (k + 1)
        Pn_1 <- Pn
        Pn <- P_next
    }
    return(Pn)
}
