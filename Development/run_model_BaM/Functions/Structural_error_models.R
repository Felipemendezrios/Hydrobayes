#' Structural error model
#'
#' Linear depending on space covariate for WSE and Gaussian in discharge and velocity
#'
#' @param Ysim numeric vector, simulated values
#' @param Yobs numeric vector, observed values
#' @param Yu numeric vector, uncertainty in observed values
#' @param gamma data frame, parameters of structural error model
#' @param x numeric vector, space covariate
#'
#' @return real value, loglikelihood computation
#'
#' @export
llfunk_Linear_X_in_Z_Gaussian_QV <- function(Ysim, Yobs, Yu, gamma, x) {
    p <- NCOL(Ysim)
    out <- 0
    for (i in 1:p) {
        if (i == 1) {
            g0 <- gamma[2 * i - 1]
            g1 <- gamma[2 * i]
            s <- abs(g0) + abs(g1) * abs(x)
        } else {
            s <- gamma[1 + i] # Case Q and V
        }
        m <- Ysim[, i]
        ps <- dnorm(
            x = Yobs[, i],
            mean = m,
            sd = sqrt(s^2 + Yu[, i]^2),
            log = TRUE
        )
        out <- out + sum(ps, na.rm = TRUE)
    }
    return(out)
}
