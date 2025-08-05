logPrior_Verification <- function(prior_param_case) {
    out <- 0
    nb_param <- length(prior_param_case)
    for (i in 1:nb_param) {
        initial_guess <- prior_param_case[[i]]$init
        distribution <- prior_param_case[[i]]$prior$dist
        distribution_param <- prior_param_case[[i]]$prior$par

        out_temp <-
            if (distribution == "FIX") {
                # Default values: only for not avoiding to accept all samples
                dlnorm(initial_guess,
                    meanlog = log(initial_guess),
                    sdlog = 0.2,
                    log = TRUE
                )
            } else if (
                distribution == "FlatPrior" | distribution == "FlatPrior+" | distribution == "FlatPrior-"
            ) {
                0
            } else if (
                distribution == "Gaussian"
            ) {
                dnorm(initial_guess,
                    mean = distribution_param[1],
                    sd = distribution_param[2],
                    log = TRUE
                )
            } else if (
                distribution == "Uniform"
            ) {
                dunif(initial_guess,
                    min = distribution_param[1],
                    max = distribution_param[2],
                    log = TRUE
                )
            } else if (
                distribution == "LogNormal"
            ) {
                dlnorm(initial_guess,
                    meanlog = distribution_param[1],
                    sdlog = distribution_param[2],
                    log = TRUE
                )
            } else {
                stop(paste0("The Prior distribution: ", distribution, " is not supported"))
            }
        out <- out + out_temp
    }
    if (is.na(out) | is.infinite(out)) {
        out <- -Inf
    }
    return(out)
}
