Input_SU_key_reachesMage <- list(
    c(1, 3),
    2,
    4
)

SU <- vector(mode = "list", length = length(Input_SU_key_reachesMage))
# SU1
pk <- c(1:7)
bief <- c(1, 1, 3, 3, 3, 3, 3)
Z <- cbind(c(1, 1, 1, 1, 0, 0, 0), c(0, 0, 0, 0, 1, 1, 1)) # or getCovariate_piecewise(pk,shiftPoints)
prior <- list(
    RBaM::parameter(name = "alpha1", init = 35),
    RBaM::parameter(name = "alpha2", init = 55)
)
SU[[1]] <- list(pk = pk, bief = bief, Z = Z, prior = prior)
