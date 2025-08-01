llfunk_Felipe <- function(Ysim, Yobs, Yu, gamma, x) {
    p <- NCOL(Ysim)
    out <- 0
    for (i in 1:p) {
        if (i == 1) {
            g0 <- gamma[2 * i - 1]
            g1 <- gamma[2 * i]
            s <- g0 + g1 * abs(x)
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
