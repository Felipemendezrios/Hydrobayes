logPrior <- function(prior_param_case) {
    n_type_param <- length(prior_param_case)
    out <- 0
    for (j in 1:n_type_param) {
        nb_param <- length(prior_param_case[[j]][[1]])
        for (i in 1:nb_param) {
            initial_guess <- prior_param_case[[j]][[1]]
            distribution <- prior_param_case[[j]][[2]]
            distribution_param <- prior_param_case[[j]][[3]]

            out_temp <-
                if (distribution[[i]] == "FIX") {
                    # Default values: only for not avoiding to accept all samples
                    dlnorm(initial_guess[[i]],
                        meanlog = log(initial_guess[[i]]),
                        sdlog = 0.2,
                        log = TRUE
                    )
                } else if (
                    distribution[[i]] == "FlatPrior" | distribution[[i]] == "FlatPrior+" | distribution[[i]] == "FlatPrior-"
                ) {
                    0
                } else if (
                    distribution[[i]] == "Gaussian"
                ) {
                    dnorm(initial_guess[[i]],
                        mean = distribution_param[[i]][1],
                        sd = distribution_param[[i]][2],
                        log = TRUE
                    )
                } else if (
                    distribution[[i]] == "Uniform"
                ) {
                    dunif(initial_guess[[i]],
                        min = distribution_param[[i]][1],
                        max = distribution_param[[i]][2],
                        log = TRUE
                    )
                } else if (
                    distribution[[i]] == "LogNormal"
                ) {
                    dlnorm(initial_guess[[i]],
                        meanlog = distribution_param[[i]][1],
                        sdlog = distribution_param[[i]][2],
                        log = TRUE
                    )
                } else {
                    stop(paste0("The Prior distribution: ", distribution[[i]], " is not supported"))
                }
            out <- out + out_temp
        }
    }
    if (is.na(out) | is.infinite(out)) {
        out <- -Inf
    }
    return(out)
}
