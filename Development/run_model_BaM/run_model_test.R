rm(list = ls())


library(RBaM)
library(stringr)
library(reshape2)
library(ggplot2)

# setPathToBaM("/home/famendezrios/Documents/Git/BaM_dev/makefile/")

# path_results = '/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/Development/run_model_RBaM/HHLab_manually_structure_estimation/'
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Example of using function runModel() ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

n_degree <- 1

# Do not touch: BaM model specification
nY <- 3
nX <- 4
param_theta <- list(
    Kmin.init = 1 / 0.010,
    kmin.distri = "FlatPrior+",
    Kmoy.init = 20,
    kmoy.distri = "FIX"
)
# param_error_model <- list(
#     ## Stage model error prior information
#     z.ini = 0.5,
#     z.prior.dist = "Uniform",
#     z.prior.par = c(0, 1)
# )
# Defining theta parameter list
n_Covariates_flood <- 1
n_Covariates_min <- (n_degree + 1)
theta_param <- vector(
    mode = "list",
    length = n_Covariates_min + n_Covariates_flood
)

for (i in 1:n_Covariates_min) {
    if (i == 1) {
        Param_min <- RBaM::parameter(
            name = paste0("a", i - 1, "_min"),
            init = param_theta$Kmin.init,
            prior.dist = param_theta$kmin.distri,
            prior.par = param_theta$Kmin.prior.par
        )
    } else {
        Param_min <- RBaM::parameter(
            name = paste0("a", i - 1, "_min"),
            init = 0,
            prior.dist = "FlatPrior"
        )
    }
    theta_param[[i]] <- Param_min
}

for (i in 1:n_Covariates_flood) {
    Param_flood <- RBaM::parameter(
        name = paste0("a", i - 1, "_flood"),
        init = param_theta$Kmoy.init,
        prior.dist = param_theta$kmoy.distri
    )
    theta_param[[n_Covariates_min + i]] <- Param_flood
}

ID <- "MAGE_ZQV"
MAGE_executable <- "/home/famendezrios/Documents/Softwares/pamhyr2/mage8/mage"
MageVersion <- "8"

path_model_mage <- "/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/Development/run_model_RBaM/run_model_BaM/"

mageDir <- c(
    file.path(path_model_mage, "model_mage_event_main/"),
    file.path(path_model_mage, "model_mage_event_2/"),
    file.path(path_model_mage, "model_mage_event_3/"),
    file.path(path_model_mage, "model_mage_event_4/")
)

dir_REPFile <- list.files(mageDir[1], pattern = "\\.REP$", ignore.case = TRUE)

# Spatialisation:

# Total number of discretization points required for Legendre polynomial calculations
total_points <- 200

# Boundary points !!! need to be adapt in multi-reach case
start_KP <- 0
end_KP <- 16.504

boundaries_points <- c(start_KP, end_KP)

interpolation_specific_points <- function(total_points = 100,
                                          all_specific_points) {
    if (total_points <= length(all_specific_points)) stop("Total points defined is lower or equal to the specific points required. Please either increase the number of total_points or not used interpolation function")
    # Calculate the intervals between the specific points
    intervals <- diff(all_specific_points)
    # Distribute the number of points proportionally to each interval
    total_local_points <- total_points - 1 # last value during diff function
    points_per_interval_float <- (intervals / sum(intervals)) * total_local_points

    points_per_interval <- round(points_per_interval_float)
    # Due to round function, we need to check the size
    # Check if more points than expected, need to remove some of them
    if (sum(points_per_interval) > total_local_points) {
        id_order <- order(points_per_interval_float %% 1)
        i <- 1
        while (sum(points_per_interval) > total_local_points) {
            id_position <- which(i == id_order)
            points_per_interval[id_position] <- points_per_interval[id_position] - 1
            i <- i + 1
        }
        # Check if less points than expected, need to add some of them
    } else if (sum(points_per_interval) < total_local_points) {
        id_order <- order(points_per_interval_float %% 1)
        i <- 1
        while (sum(points_per_interval) < total_local_points) {
            id_position <- which(i == id_order)
            points_per_interval[id_position] <- points_per_interval[id_position] + 1
            i <- i + 1
        }
    }

    # Generate sub-sequences for each interval, excluding the duplicate endpoint
    specific_seq_temp <- unlist(mapply(
        function(start, end, n) seq(start, end, length.out = n + 1)[-(n + 1)], # Exclude the endpoint
        all_specific_points[-length(all_specific_points)], # Start of each interval
        all_specific_points[-1], # End of each interval
        points_per_interval # Number of points per interval
    ))

    #####  reviser le fait que j'ai à la fin 101 valeurs dans la séquance, alors que je demande 100 points, mais les boundary font le bordel
    # Add the final endpoint manually
    specific_seq <- c(specific_seq_temp, all_specific_points[length(all_specific_points)])

    if (!all(all_specific_points %in% specific_seq)) stop("Distance so close, increase number of interpolate KP to get all points required")

    if (length(specific_seq) != total_points) stop("Something is wrong with the interpolation")
    return(specific_seq)
}

# Get the covariant discretization
grid_covariant_discretized <- interpolation_specific_points(
    total_points = total_points,
    all_specific_points = boundaries_points
)

# Set Z file : cov_bar for the legendre polynomials
min_cov_legendre <- min(grid_covariant_discretized)
max_cov_legendre <- max(grid_covariant_discretized)
cov_bar_legendre <- 2 * (grid_covariant_discretized - min_cov_legendre) / (max_cov_legendre - min_cov_legendre) - 1
# Write Legendre vector:
write.table(
    data.frame(
        KP_original = grid_covariant_discretized,
        KP_normalized = cov_bar_legendre
    ),
    file = file.path(tempdir(), "vector_legendre.txt"), row.names = F
)
# Write co variant matrix for main channel
# Create an empty data frame and fill it with polynomial values
legendre_df <- data.frame(x = cov_bar_legendre) # Start with x values
########
# Generate files
########
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

for (n in 0:n_degree) {
    legendre_df[[paste0("P", n)]] <- getlegendre(
        degree = n,
        normalized_values = cov_bar_legendre
    )
}

zFileKmin <- file.path(tempdir(), "Zfile_Kmin.txt")
zFileKmoy <- file.path(tempdir(), "Zfile_Kflood.txt")

if (ncol(legendre_df) == 2) { # if n degree is equal to 0
    matrix_zFileKmin <- legendre_df["P0"]
} else {
    matrix_zFileKmin <- legendre_df[, -1]
}
matrix_zFileKmoy <- legendre_df["P0"]


xtra <- xtraModelInfo(
    fname = "Config_setup.txt",
    object = list(
        exeFile = MAGE_executable,
        version = MageVersion,
        mageDir = mageDir,
        repFile = dir_REPFile,
        zKmin = matrix_zFileKmin,
        zFileKmin = zFileKmin,
        doExpKmin = FALSE,
        zKmoy = matrix_zFileKmoy,
        zFileKmoy = zFileKmoy,
        doExpKmoy = FALSE
    )
)
# Assemble Model object
mod <- model(
    ID = ID,
    nX = nX,
    nY = nY, ,
    par = theta_param,
    xtra = xtra
)

n_events <- 4
names_CalData <- c(
    "case_60_60",
    "case_30_30",
    "case_15_15",
    "case_07_07"
)

# Information to be given by event!
n_obs <- c(8, 8, 8, 8)
specific_reach <- c(1, 1, 1, 1)

# Space-temporal coordinates
x_coordinates <- list(
    c(0.060, 0.650, 1.650, 3.650, 5.650, 7.650, 9.650, 11.650),
    c(0.060, 0.650, 1.650, 3.650, 5.650, 7.650, 9.650, 11.650),
    c(0.060, 0.650, 1.650, 3.650, 5.650, 7.650, 9.650, 11.650),
    c(0.060, 0.650, 1.650, 3.650, 5.650, 7.650, 9.650, 11.650)
)
t_coordinates <- c(3420, 3420, 3420, 3420)

if (any(c(
    length(n_obs),
    length(specific_reach),
    length(x_coordinates),
    length(t_coordinates),
    length(names_CalData)
)
!= n_events)) {
    stop(paste0("Information must be given by event, donc all the calibration data must have a length of ", n_events))
}
events <- rep(seq_len(n_events), times = n_obs)
reaches <- rep(specific_reach, times = n_obs)

if (is.list(x_coordinates)) {
    x <- unlist(x_coordinates)
} else {
    x <- rep(x_coordinates, times = n_obs)
}

if (is.list(t_coordinates)) {
    t <- unlist(t_coordinates)
} else {
    t <- rep(t_coordinates, times = n_obs)
}

X <- data.frame(
    events = events,
    reaches = reaches,
    x = x,
    t = t
)
###################################

# Run model
Ysim <- runModel(
    workspace = tempdir(),
    mod = mod,
    X = X,
    stout = NULL
)

results_sim <- cbind(X, Ysim)

# Plot WSE residuals
ggplot(results_sim, aes(x = x, color = factor(events))) +
    geom_point(aes(y = Y1)) +
    labs(
        title = "WSE simulation",
        x = "Position (m)",
        y = "Water surface elevation (m)",
        color = "Events"
    ) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), )

# Plot Q residuals
ggplot(results_sim, aes(x = x, color = factor(events))) +
    geom_point(aes(y = Y2)) +
    labs(
        title = "Discharge simulations",
        x = "Position (m)",
        y = "Discharge (m3/s)",
        color = "Events"
    ) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Example of computing the log-posterior ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Define a prior function, which takes as input a numeric vector (parameter values)
# and return a value (the prior log-density, or -Inf if prior is 0 or unfeasible)

# logPrior_Flat()

#### Get calibration data
# User input:
## Give the information to use during calibration

# Load calibration data:
load(file.path(
    "/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/processed_data/HHLab/smooth_bed/",
    "data_HHLab_all_cases.RData"
))
all_data_calibration <- list(WS_profiles = WS_profiles)

# Extract the WS_profiles list from your data
ws_profiles <- all_data_calibration$WS_profiles

metadata <- cbind(name_event = rep(names_CalData, times = n_obs), X)

columns_interested <- c("z_mean", "Yu")

# Initialize an empty list to store extracted rows
CalibrationData <- data.frame(X,
    Y_Z = NA,
    Y_Q = NA,
    Y_V = NA,
    Yu_Z = NA,
    Yu_Q = NA,
    Yu_V = NA
)

for (i in seq_len(nrow(metadata))) {
    event_name <- metadata$name_event[i]
    x_value <- metadata$x[i]

    # Access the correct data frame inside WS_profiles
    df <- ws_profiles[[event_name]]

    # Filter rows where x matches x_value (exact match)
    # If floating point comparison is an issue, consider a tolerance
    filtered_df <- df[df$x == x_value, columns_interested]

    # Replace NA values
    CalibrationData[i, c("Y_Z", "Yu_Z")] <- filtered_df
}

# Harcode for WSE
Yobs <- CalibrationData[, c("Y_Z", "Y_Q", "Y_V")]
Yu <- CalibrationData[, c("Yu_Z", "Yu_Q", "Yu_V")]

# Select the log-likelihood function from the available RBaM catalog.
# You can also define your own log-likelihood function, which takes as inputs Ysim, Yobs, Yu
# (3 data frames with identical sizes), gamma (vector of structural error parameters),
# and returns a value (the log-likelihood, -Inf if the likelihood is zero or unfeasible).
# options : llfunk_iid_Gaussian ; llfunk_iLinear_Gaussian
# inference functions and model

# The logPosterior function returns a list containing log-post, log-prior, log-lkh and simulated values
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

logLikelihood <- llfunk_Felipe

n_param_gamma <- if (identical(logLikelihood, llfunk_iid_Gaussian)) {
    1
} else if (identical(logLikelihood, llfunk_iLinear_Gaussian)) {
    2
} else {
    stop("Unknown log-likelihood function")
}

########### Adapt where observational data are WSE,Q, V with different structural error model. For the moment, constant for all of them is assigned

# number_parameters <- length(mod$par) + nY * n_param_gamma
number_parameters <- length(mod$par) + 4


# Add check to give a value for each parameter (using number_parameters)
logPrior <- function(parvector) {
    out <-
        # Priors for model parameters
        # dlnorm(parvector[1], meanlog = log(100), sdlog = 0.1, log = TRUE) + # a0_min
        0 + # a0_min
        0 + # a1_min
        dlnorm(parvector[3], meanlog = log(20), sdlog = 0.2, log = TRUE) + # a0_flood
        # Priors for structural error parameters for each output variable !
        # WSE simulation
        0 + # gamma0
        0 + # gamme1
        # Discharge simulation
        dlnorm(parvector[6], meanlog = log(10), sdlog = 0.2, log = TRUE) + # gamma0
        # Velocity simulation
        dlnorm(parvector[7], meanlog = log(10), sdlog = 0.2, log = TRUE) # gamma0
    if (is.na(out) | is.infinite(out)) {
        out <- -Inf
    }
    return(out)
}

# Compute the posterior log-density at some parameter values
parvector <- c(
    RBaM::getInitPar(mod$par), # initial parameter values for the model parameters
    c(0.001, 0, 10, 10) # Initial values for the structural error parameters (gamma) for each output variable
)

logpost <- logPosterior_BaM(
    parvector = parvector, # parameter values at which the posterior is evaluated
    X = X, Yobs = Yobs, Yu = Yu, # calibration data
    lpfunk = logPrior, llfunk = logLikelihood, mod = mod, llargs = X$x
)
logpost

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Example of maximizing the log-posterior ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Just use function optim() coming with R base installation
maxpost <- stats::optim(
    par = parvector, # starting point
    fn = logPosterior_BaM_wrapped, # function to be optimized. It's a wrapped version of function logPosterior, returning the log-post value only. See ?logPosterior_wrapped
    control = list(fnscale = -1), # tells optim to maximize rather than minimize
    X = X, Yobs = Yobs, Yu = Yu, lpfunk = logPrior, llfunk = logLikelihood, mod = mod,
    llargs = X$x
) # arguments passed to logPosterior_wrapped
maxpost$par
maxpost$convergence

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Example of MCMC-sampling the log-posterior ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
x0 <- maxpost$par
names(x0) <- c(
    RBaM::getNames(mod$par),
    "g0_Y1", "g1_Y1", "g0_Y2", "g0_Y3"
)
# Adaptive Metropolis
ptm <- proc.time()
mcmc <- MCMC_AM(
    logPdf = logPosterior_BaM, x0 = x0,
    X = X, Yobs = Yobs, Yu = Yu, lpfunk = logPrior, llfunk = logLikelihood, mod = mod, llargs = X$x
) # arguments passed to logPosterior
proc.time() - ptm
pairs(cbind(mcmc$samples, mcmc$components$logPosterior))

library(patchwork)

plots <- tracePlot(mcmc$samples)

wrap_plots(plots, ncol = 3)


plots <- tracePlot(mcmc$components)
wrap_plots(plots, ncol = 3)





# One-At-A-Time Metropolis, no speed-up
# ptm <- proc.time()
# mcmc_OAAT <- MCMC_OAAT(
#     logPdf = logPosterior_BaM, x0 = x0,
#     X = X, Yobs = Yobs, Yu = Yu, lpfunk = logPrior, llfunk = logLikelihood, mod = mod
# ) # arguments passed to logPosterior
# proc.time() - ptm
# pairs(cbind(mcmc_OAAT$samples, mcmc_OAAT$components$logPosterior))

# tracePlot(mcmc$samples)
# tracePlot(mcmc$components)

# # One-At-A-Time Metropolis, with speed-up
# ptm <- proc.time()
# mcmc_OAAT_speed_up <- MCMC_OAAT(
#     logPdf = logPosterior_BaM, x0 = x0, nTheta = length(x0) - 2,
#     X = X, Yobs = Yobs, Yu = Yu, lpfunk = logPrior, llfunk = logLikelihood, mod = mod
# ) # arguments passed to logPosterior
# proc.time() - ptm
# pairs(cbind(mcmc_OAAT_speed_up$samples, mcmc_OAAT_speed_up$components$logPosterior))
