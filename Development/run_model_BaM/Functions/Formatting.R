#' Set the prior on the model's parameters in BaM format
#'
#' @param n_degree_Kmin integer, Legendre polynomial degree in the main channel
#' @param n_degree_Kflood integer, Legendre polynomial degree in the floodplain
#' @param param_model_Prior_specification list, prior information on the model's parameters : see`Details`
#'
#' @return list, prior information on the parameters in BaM formatting
#'
#' @details
#' param_model_Prior_specification must contain the while information about parameters in the main channel and floodplain:
#' - Initial guess
#' - Prior distribution
#' - Parameter of each prior distribution
#'
#' @export
set_param_model <- function(
    n_degree_Kmin,
    n_degree_Kflood,
    param_model_Prior_specification) {
    theta_param <- vector(
        mode = "list",
        length = n_degree_Kmin + n_degree_Kflood
    )

    for (i in 1:(n_degree_Kmin + 1)) {
        theta_param[[i]] <- RBaM::parameter(
            name = paste0("a", i - 1, "_min"),
            init = param_model_Prior_specification$Kmin.init[i],
            prior.dist = param_model_Prior_specification$Kmin.distri[i],
            prior.par = param_model_Prior_specification$Kmin.prior.par[[i]]
        )
    }

    for (i in 1:(n_degree_Kflood + 1)) {
        theta_param[[(n_degree_Kmin + 1) + i]] <- RBaM::parameter(
            name = paste0("a", i - 1, "_flood"),
            init = param_model_Prior_specification$Kflood.init[i],
            prior.dist = param_model_Prior_specification$Kflood.distri[i],
            prior.par = param_model_Prior_specification$Kflood.prior.par[[i]]
        )
    }
    return(theta_param)
}
