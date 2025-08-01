rm(list = ls())
graphics.off()

# ~~~~~~~~~~~~~~~~~~~ #
# Do not touch:
# ~~~~~~~~~~~~~~~~~~~ #

# Libraries
packages <- c("RBaM", "stringr", "reshape2", "ggplot2", "here")

# Set path if RBaM is updated or first installation of RBaM
# setPathToBaM("/home/famendezrios/Documents/Git/BaM_dev/makefile/")

# Source functions
# Set the path to your folder
dir_functions <- "Development/run_model_BaM/Functions"

# List all .R files in the folder
folder_functions <- list.files(path = dir_functions, pattern = "\\.R$", full.names = TRUE)

# Source each file
invisible(lapply(folder_functions, source))

invisible(lapply(packages, install_and_load))

# Root git folder:
repo_root <- here::here()

# Do not touch: MAGE model specification
ID <- "MAGE_ZQV"
MAGE_executable <- "/home/famendezrios/Documents/Softwares/pamhyr2/mage8/mage"
MageVersion <- "8"
nY <- 3
nX <- 4
# ~~~~~~~~~~~~~~~~~~~ #
# End hardcode data
# ~~~~~~~~~~~~~~~~~~~ #

# ~~~~~~~~~~~~~~~~~~~ #
# User settings:
# ~~~~~~~~~~~~~~~~~~~ #

# ~~~~~~~~~~~~~~~~~~~ #
# Define case of study and calibration data:


# Number of event to be computed and observed
n_events <- 4

# Catalog:
# HHLab/smooth_bed and dataset: data_HHLab_all_cases.RData
case_study <- "HHLab/smooth_bed"
dataset <- "data_HHLab_all_cases.RData"
# CalData_measures: WSE , Discharge, Velocity
CalData_measures <- c("WSE")

# Calibration Data: Case-specific

if (case_study == "HHLab/smooth_bed") {
    variable_observed <- c(
        "WSE",
        "WSE",
        "WSE",
        "WSE"
    )
    names_CalData <- c(
        "case_60_60",
        "case_30_30",
        "case_15_15",
        "case_07_07"
    )

    columns_interested <- c("z_mean", "Yu")
}
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

# ~~~~~~~~~~~~~~~~~~~ #
# Define simulations settings:

# Set MAGE model directory
# Must be part of the structure of the git repository and end by '\'
models <- "Development/run_model_BaM/Model_Mage_events/"
path_model_mage <- file.path(repo_root, models)

# Do not forget to put '/' at the end of the path
mage_folders_events <- c(
    "model_mage_event_main/",
    "model_mage_event_2/",
    "model_mage_event_3/",
    "model_mage_event_4/"
)

# ~~~~~~~~~~~~~~~~~~~ #
# Define estimation settings:

# Legendre polynomial degree in the main channel
n_degree_Kmin <- 1
n_degree_Kflood <- 0

# Input parameters in main channel:
Kmin_prior <- list(
    # Initial guess by parameter in main channel
    Kmin.init = c(
        1 / 0.010, # a0_min
        0 # a1_min
    ),
    # Prior distribution by parameter in main channel
    Kmin.distri = c(
        "FlatPrior+",
        "FlatPrior"
    ),
    # Parameters of the prior distribution by parameter in main channel
    Kmin.prior.par = list(
        NULL,
        NULL
    )
)

# Input parameters in the floodplain:
Kflood_prior <- list(
    # Initial guess by parameter in the floodplain
    Kflood.init = c(20),
    # Prior distribution by parameter in the floodplain
    Kflood.distri = c("FIX"),
    # Parameters of the prior distribution by parameter in the floodplain
    Kflood.prior.par = list(
        NULL
    )
)

# Spatialisation:
# Total number of discretization points required for Legendre polynomial calculations
total_points <- 200

# Boundary points !!! need to be adapt in multi-reach case
start_KP <- 0
end_KP <- 16.504

# ~~~~~~~~~~~~~~~~~~~ #
# Save temporal results in a temporal internal folder?
# If TRUE: internal temporal folder is used
# Otherwise, a folder will be created in the workspace to save the results
temporal_SaveResults <- FALSE
remove_oldFiles <- TRUE

if (temporal_SaveResults) {
    path_temporal_results <- tempdir()
} else {
    # Be careful that new temporal folder is root in repo_root. In other words, the root is the git repository
    path_temporal_results <- file.path(repo_root, "Development/temporal_results")

    # If the folder exists, remove all files inside
    if (dir.exists(path_temporal_results)) {
        old_files <- list.files(path_temporal_results, full.names = remove_oldFiles)
        file.remove(old_files)
    } else {
        dir.create(path_temporal_results, recursive = TRUE)
    }
}

do.plot <- FALSE

# ~~~~~~~~~~~~~~~~~~~ #
# End user setting
# ~~~~~~~~~~~~~~~~~~~ #

# ~~~~~~~~~~~~~~~~~~~ #
# Calibration Data processing
# ~~~~~~~~~~~~~~~~~~~ #

file_case_study <- file.path(
    "data/processed_data",
    case_study,
    dataset
)

if (!file.exists(file_case_study)) stop(paste0("Case of study:", case_study, " is not supported yet. Ensure that the case study is saved in data/processed_data"))


# Load into a temporary environment
tmp_env <- new.env()
loaded_names <- load(file_case_study, envir = tmp_env)

# Check if the variable name matches what you expect
if (loaded_names != "all_data_calibration") {
    stop("Unexpected variable name in file. Expected: all_data_calibration, but found: ", paste(loaded_names, collapse = ", "))
}
load(file_case_study)
sample_data_calibration <- list(
    WSE = NA,
    Discharge = NA,
    Velocity = NA
)

for (i in which(CalData_measures %in% names(all_data_calibration))
) {
    sample_data_calibration[[i]] <- all_data_calibration[[i]]
}

# Mage directory
mageDir <- file.path(path_model_mage, mage_folders_events)

# Check model projects
if (!any(dir.exists(mageDir))) stop(paste0("The MAGE directories do not exist"))

# Check if number of event is correct with defined data
if (any(c(
    length(mageDir),
    length(n_obs),
    length(specific_reach),
    length(x_coordinates),
    length(t_coordinates),
    length(names_CalData),
    length(variable_observed)
)
!= n_events)) {
    stop(paste0("Information must be given by event, donc all the calibration data must have a length of ", n_events))
}

# ~~~~~~~~~~~~~~~~~~~ #
# Estimation processing
# ~~~~~~~~~~~~~~~~~~~ #

# Check number of prior given according to n_degree chosen
if (!all(lengths(Kmin_prior) == (n_degree_Kmin + 1))) stop(paste0(n_degree_Kmin + 1, " prior data must be specified by parameter in main channel"))

# Check number of prior given according to n_degree chosen
if (!all(lengths(Kflood_prior) == (n_degree_Kflood + 1))) stop(paste0(n_degree_Kflood + 1, " prior data must be specified by parameter in main channel"))

param_model_Prior_specification <- c(Kmin_prior, Kflood_prior)

param_model_BaMPrior <- set_param_model(
    n_degree_Kmin = n_degree_Kmin,
    n_degree_Kflood = n_degree_Kflood,
    param_model_Prior_specification = param_model_Prior_specification
)

# Check
if (length(param_model_BaMPrior) != ((n_degree_Kflood + 1) + (n_degree_Kmin + 1))) stop("The number of polynomial degree in main channel and floodplain must be coincide with the number of prior data given")

###
dir_REPFile <- list.files(mageDir[1], pattern = "\\.REP$", ignore.case = TRUE)

# Boundary points
boundaries_points <- c(start_KP, end_KP)

# Get the covariant discretization
grid_covariant_discretized <- interpolation_specific_points(
    total_points = total_points,
    all_specific_points = boundaries_points
)
# Spatialisation:
# Set Z file : cov_bar for the Legendre polynomials
min_cov_legendre <- min(grid_covariant_discretized)
max_cov_legendre <- max(grid_covariant_discretized)
cov_bar_legendre <- 2 * (grid_covariant_discretized - min_cov_legendre) / (max_cov_legendre - min_cov_legendre) - 1
# Write Legendre vector:
write.table(
    data.frame(
        KP_original = grid_covariant_discretized,
        KP_normalized = cov_bar_legendre
    ),
    file = file.path(path_temporal_results, "vector_legendre.txt"), row.names = F
)

# Write co variant matrix for main channel
# Create an empty data frame and fill it with polynomial values
legendre_df_covariate <- data.frame(x = cov_bar_legendre) # Start with x values
legendre_df_Kmin <- legendre_df_Kflood <- legendre_df_covariate

# Spatial distributed friction in the main channel
for (n in 0:n_degree_Kmin) {
    legendre_df_Kmin[[paste0("P", n)]] <- getlegendre(
        degree = n,
        normalized_values = cov_bar_legendre
    )
}
# Spatial distributed friction in the floodplain
for (n in 0:n_degree_Kflood) {
    legendre_df_Kflood[[paste0("P", n)]] <- getlegendre(
        degree = n,
        normalized_values = cov_bar_legendre
    )
}

zFileKmin <- file.path(path_temporal_results, "Zfile_Kmin.txt")
zFileKmoy <- file.path(path_temporal_results, "Zfile_Kflood.txt")

if (ncol(legendre_df_Kmin) == 2) { # if n degree is equal to 0
    matrix_zFileKmin <- legendre_df_Kmin["P0"]
} else {
    matrix_zFileKmin <- legendre_df_Kmin[, -1]
}

if (ncol(legendre_df_Kflood) == 2) { # if n degree is equal to 0
    matrix_zFileKmoy <- legendre_df_Kflood["P0"]
} else {
    matrix_zFileKmoy <- legendre_df_Kflood[, -1]
}

# ~~~~~~~~~~~~~~~~~~~ #
# Xtra settings
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
# ~~~~~~~~~~~~~~~~~~~ #
# Assemble Model object
mod <- model(
    ID = ID,
    nX = nX,
    nY = nY, ,
    par = param_model_BaMPrior,
    xtra = xtra
)

# ~~~~~~~~~~~~~~~~~~~ #
# Calibration data processing
# ~~~~~~~~~~~~~~~~~~~ #
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

# Run model
Ysim <- runModel(
    workspace = path_temporal_results,
    mod = mod,
    X = X,
    stout = NULL
)

results_sim <- cbind(X, Ysim)

# ~~~~~~~~~~~~~~~~~~~ #
# Some plots
# ~~~~~~~~~~~~~~~~~~~ #
if (do.plot) {
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
}

# ~~~~~~~~~~~~~~~~~~~ #
# Computing structural error model
# ~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Example of computing the log-posterior ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Define a prior function, which takes as input a numeric vector (parameter values) and return a value (the prior log-density, or -Inf if prior is 0 or unfeasible)

metadata <- cbind(
    variable_observed = rep(variable_observed, times = n_obs),
    name_event = rep(names_CalData, times = n_obs),
    X
)

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
    variable_observed_event <- metadata$variable_observed[i]
    event_name <- metadata$name_event[i]
    x_value <- metadata$x[i]

    # Access the correct data frame inside WS_profiles
    df <- sample_data_calibration[[variable_observed_event]][[event_name]]

    # Filter rows where x matches x_value (exact match)
    # If floating point comparison is an issue, consider a tolerance
    filtered_df <- df[df$x == x_value, columns_interested]

    # Replace NA values
    CalibrationData[i, c("Y_Z", "Yu_Z")] <- filtered_df
}

# Get output variable with their uncertainties
Yobs <- CalibrationData[, c("Y_Z", "Y_Q", "Y_V")]
Yu <- CalibrationData[, c("Yu_Z", "Yu_Q", "Yu_V")]

# Select the log-likelihood function from the available RBaM catalog.

# You can also define your own log-likelihood function, which takes as inputs Ysim, Yobs, Yu

# (3 data frames with identical sizes), gamma (vector of structural error parameters), and returns a value (the log-likelihood, -Inf if the likelihood is zero or unfeasible).

# options : llfunk_iid_Gaussian ; llfunk_iLinear_Gaussian
# inference functions and model

# The logPosterior function returns a list containing log-post, log-prior, log-lkh and simulated values

# TESTS ERROR MODELS : INPUT
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
