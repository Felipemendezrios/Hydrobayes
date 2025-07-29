rm(list = ls())
library(RBaM)
library(stringr)
library(reshape2)
library(ggplot2)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Example of using function runModel() ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

calibration_type <- "Q" # 'WSE' or 'Q' (water surface elevation or discharge)

path_calibration_data <- ifelse(calibration_type == "WSE",
    "/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/Development/run_model_RBaM/Calibration_water_surface_profiles/BaM_Model/n_2",
    "/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/Development/run_model_RBaM/Calibration_time_series/BaM_Model/n_2"
)
path_model_mage_global <- "/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/Development/run_model_RBaM/model_mage"

# Define parameters - only name and initial value are needed, no need to specify priors
n_degree <- 2

if (n_degree != 2 | basename(path_calibration_data) != "n_2") {
    stop("This example is only for n_2 case, please change the path_calibration_data to n_2")
}
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


# Xtra config:

check_cal_WS_profiles <- calibration_type == "WSE"
# Create a sequence of numbers: time fixed to extract simulation
if (check_cal_WS_profiles) {
    # Create a sequence of numbers: time fixed to extract simulation (Y)
    sequence <- seq(0.95, 3.95, by = 1)
    # Time mapping representing the order of constant discharge (id_case) simulation and the dataset used. All consistent with the MAGE model!!
    # For example, 'case_60_60' means Q = 120 L/s in [Dataset Proust et al., (2022)](https://doi.org/10.57745/EQURJN)
    time_mapping_temp <- data.frame(
        # id from dataset order (flow steps configuration)
        id_case = c(
            "case_60_60",
            "case_30_30",
            "case_15_15",
            "case_07_07"
        )
    )

    sequence_all <- c(
        sequence
    )


    # Find the maximum number of decimal places
    max_decimals <- max(sapply(sequence_all, function(x) nchar(sub(".*\\.", "", as.character(x)))))

    sequence_all <- sort(round(sequence_all, max_decimals))

    # Mage extraire arguments
    mage_extraire_args <- lapply(
        sequence_all,
        function(x) {
            paste0(
                "ZdX 1 h",
                trimws(format(x,
                    nsmall = max_decimals
                ))
            )
        }
    )

    # Position to put the observed data in the grid
    id_fixed <- data.frame(
        id_fixed = c(
            0.060,
            0.650,
            1.650,
            3.650,
            5.650,
            7.650,
            9.650,
            11.650
        )
    )
} else {
    # Create a sequence: position fixed to extract simulation (Y)
    # If bug would detect, it could be the number of decimals imposed by mage ?
    sequence_all <- c(
        0.060,
        0.650,
        1.650,
        3.650,
        5.650,
        7.650,
        9.650,
        11.650
    )
    # Find the maximum number of decimal places
    max_decimals <- max(sapply(sequence_all, function(x) nchar(sub(".*\\.", "", as.character(x)))))

    sequence_all <- sort(round(sequence_all, max_decimals))

    # Mage extraire arguments
    mage_extraire_args <- lapply(
        sequence_all,
        function(x) {
            paste(
                "ZdT 1",
                format(x, nsmall = max_decimals)
            )
        }
    )


    # Time to introduce calibration data (in hours)
    # Time mapping between cases with constant discharge (id_case) and simulation time (id_fixed) for extraction, all consistent with the MAGE model!!
    time_mapping_temp <- data.frame(
        # id from dataset order (flow steps configuration)
        id_case = c(
            "case_60_60",
            "case_30_30",
            "case_15_15",
            "case_07_07"
        ),
        # time to put the observed data in the grid
        id_fixed = c(
            0.95,
            1.95,
            2.95,
            3.95
        )
    )
}


# Assign the number of output for WS profiles with the time fixed
# Variable extraction
first_char <- substr(unlist(mage_extraire_args), 1, 3)

# Position or time extraction
last_number <- stringr::str_extract(unlist(mage_extraire_args), "\\d+\\.\\d+$")

extraction_data_mage_extraire_args <- paste0(first_char, "_", last_number)

path_kmin <- file.path(path_calibration_data, "Zfile_Kmin.txt")
path_kflood <- file.path(path_calibration_data, "Zfile_Kflood.txt")

matrix_zFileKmin <- read.table(path_kmin, header = TRUE, stringsAsFactors = FALSE)
matrix_zFileKmoy <- read.table(path_kflood, header = TRUE, stringsAsFactors = FALSE)
FilesNames <- list.files(path_model_mage_global, recursive = T)
dir_REPFile <- FilesNames[grep(FilesNames, pattern = ".REP")]

xtra <- xtraModelInfo(
    fname = "Config_setup.txt",
    object = list(
        exeFile = "/home/famendezrios/Documents/Softwares/pamhyr2/mage8/mage",
        version = "8",
        mageDir = paste0(path_model_mage_global, "/"),
        repFile = dir_REPFile,
        zKmin = matrix_zFileKmin,
        zFileKmin = path_kmin,
        doExpKmin = FALSE,
        zKmoy = matrix_zFileKmoy,
        zFileKmoy = path_kflood,
        doExpKmoy = FALSE,
        mage_extraire_file = "/home/famendezrios/Documents/Softwares/pamhyr2/mage8/mage_extraire",
        mage_extraire_args = unlist(mage_extraire_args)
    )
)
nY <- length(mage_extraire_args)

# Assemble Model object
mod <- model(
    ID = "MAGE_TEMP",
    nX = 1,
    nY = nY, ,
    par = theta_param,
    xtra = xtra
)
# Define input data frame
REPFile <- read.table(file.path(path_model_mage_global, "smooth_bed.REP"),
    comment.char = "*"
)

if (check_cal_WS_profiles) {
    time_mapping_temp2 <- cbind(time_mapping_temp,
        id_position = extraction_data_mage_extraire_args
    )

    time_mapping <- merge(id_fixed, time_mapping_temp2, by = NULL)

    ############################
    # Module mage
    ############################


    ### Affect input variable for BaM!: spatial grid
    #-----------------------------------------------------------#
    #                            NET FILE :                     #
    #-----------------------------------------------------------#

    REPFile <- read.table(file.path(path_model_mage_global, dir_REPFile),
        comment.char = "*"
    )

    dir.NETFile <- REPFile["NET", ]

    NETFile <- read.table(file.path(path_model_mage_global, dir.NETFile),
        comment.char = "*"
    )

    # Read dir ST file
    dir.STFiles <- NETFile[, 4]
    # Initialize the size of the list depending on the number of reaches
    Cross_sections_by_reach <- rep(list(0), nrow(NETFile))

    for (j in 1:length(dir.STFiles)) {
        STFile <- readLines(file.path(path_model_mage_global, dir.STFiles[j]))

        # Collect the KP of the ST file
        KP_elements <- c()

        # Get the first line
        first_line <- STFile[grep("^#", STFile, invert = TRUE)[1]]
        # Split values
        split_line <- strsplit(first_line, "\\s+")[[1]]
        split_line <- split_line[split_line != ""]

        KP_elements <- c(KP_elements, split_line[5])

        # Identify which lines contain "999.9990"
        marker_lines <- grep("999\\.9990", STFile)

        for (i in marker_lines) {
            if (i + 1 <= length(STFile)) {
                next_line <- STFile[i + 1]
                split_line <- strsplit(next_line, "\\s+")[[1]]
                split_line <- split_line[split_line != ""] # Remove empty entries
                if (length(split_line) == 6 || length(split_line) == 5) {
                    KP_elements <- c(KP_elements, split_line[5])
                }
            }
        }
        Cross_sections_by_reach[[j]] <- as.numeric(KP_elements)
    }

    # When a single reach is analyzed, the code run.
    # However, if it is a multiple reach case, we need to think how to handle it
    if (length(Cross_sections_by_reach) > 1) stop("Must be set for multiple reaches cases, at the moment only a single reach case is available")

    X <- data.frame(id_order = Cross_sections_by_reach[[1]])
} else {
    # Repeat and bind rows with extraction labels
    time_mapping <- do.call(rbind, lapply(extraction_data_mage_extraire_args, function(val) {
        df <- time_mapping_temp
        df$id_position <- val
        df
    }))


    ############################
    # Module mage
    ############################
    #-----------------------------------------------------------#
    #                            PAR FILE :                     #
    #-----------------------------------------------------------#
    # Read time of simulation from MAGE model (.PAR)
    SimulationTimeModel <- read.table(
        file = file.path(file.path(
            path_model_mage_global,
            FilesNames[grep(FilesNames, pattern = ".PAR")]
        )),
        sep = "",
        header = F,
        skip = 1,
        nrows = 3 # be careful timestep could be different to timestep_bin
    )
    SimulationTimeModelDF <- sapply(
        sapply(SimulationTimeModel[-3, 2],
            str_split,
            pattern = ":",
            simplify = TRUE
        ),
        as.numeric
    )
    # Get the simulation time in seconds
    for (i in 1:ncol(SimulationTimeModelDF)) {
        sim_time_temp <- SimulationTimeModelDF[, i] * c(86400, 3600, 60, 1)
        if (i == 1) {
            start_run <- sum(sim_time_temp)
        } else {
            end_run <- sum(sim_time_temp)
        }
    }
    # Get time step of the simulation
    TimeStep <- as.numeric(SimulationTimeModel[3, 2])
    # Size of observation data
    SizeObs <- (end_run - start_run) / TimeStep + 1
    # Define input data for BaM! model
    X <- data.frame(time_hours = round(
        seq(
            from = 0,
            by = TimeStep / 3600, # Get time in hours
            length.out = SizeObs
        ),
        2
    ))

    Cross_section_calibration_dT <- sapply(mage_extraire_args, function(arg) {
        if (grepl(arg, pattern = "dT")) {
            split_text <- strsplit(arg, " ")[[1]]
            return(as.numeric(split_text[length(split_text)]))
        } else {
            return(NA) # Retourne NA si le motif "dT" n'est pas trouvé
        }
    })
    # Remove NA (dX detected)
    Cross_section_calibration_dT <- Cross_section_calibration_dT[!is.na(Cross_section_calibration_dT)]
    # Remove repeated positions and get specific points defined in mage_extraire_args
    specific_points_CalStations <- unique(Cross_section_calibration_dT)
}



# Run model
Ysim <- runModel(
    workspace = tempdir(),
    mod = mod,
    X = X
    # stout = NULL
)

###################################### CHECK OK UNTIL HERE ######################################
Ysim_long <- reshape2::melt(Ysim, variable.name = "Ysim_id", value.name = "Ysim_value")
Ysim_long$id_order <- rep(X$id_order, times = ncol(Ysim))


results <- merge(X, Ysim_long, by = "id_order")


ggplot(results, aes(x = id_order, y = Ysim_value, color = Ysim_id)) +
    geom_line() +
    labs(
        title = "Simulation of water surface elevation",
        x = "Streamwise Position (m)",
        y = "Simulated Value (mm)"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(color = "Discharge cases")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Example of computing the log-posterior ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Define a prior function, which takes as input a numeric vector (parameter values)
# and return a value (the prior log-density, or -Inf if prior is 0 or unfeasible)

logPrior_Flat()

logPrior <- function(parvector) {
    out <-
        # Priors for model parameters
        # dlnorm(parvector[1], meanlog = log(100), sdlog = 0.1, log = TRUE) + # a0_min
        0 + # a0_min
        0 + # a1_min
        0 + # a2_min
        dlnorm(parvector[4], meanlog = log(20), sdlog = 0.2, log = TRUE) + # a0_flood
        # Priors for structural error parameters for each observation !
        # observation 0.95
        0 + # gamma0
        0 + # gamma1
        # observation 1.95
        0 + # gamma0
        0 + # gamma1
        # observation 2.95
        0 + # gamma0
        0 + # gamma1
        # observation 3.95
        0 + # gamma0
        0 # gamma1
    if (is.na(out) | is.infinite(out)) {
        out <- -Inf
    }
    return(out)
}

# Define calibration data
CalData <- read.table(
    file = "/home/famendezrios/Documents/These/Developpement/run_model_RBaM/CalibrationData.txt",
    header = TRUE,
    stringsAsFactors = FALSE
)

X <- CalData["id_order"]
Yobs <- CalData[c("ZdX_0.95", "ZdX_1.95", "ZdX_2.95", "ZdX_3.95")]
Yu <- CalData[c("Yu_ZdX_0.95", "Yu_ZdX_1.95", "Yu_ZdX_2.95", "Yu_ZdX_3.95")]

# Select the log-likelihood function from the available RBaM catalog.
# You can also define your own log-likelihood function, which takes as inputs Ysim, Yobs, Yu
# (3 data frames with identical sizes), gamma (vector of structural error parameters),
# and returns a value (the log-likelihood, -Inf if the likelihood is zero or unfeasible).
# options : llfunk_iid_Gaussian ; llfunk_iLinear_Gaussian

logLikelihood <- llfunk_iLinear_Gaussian

# Compute the posterior log-density at some parameter values
parvector <- c(RBaM::getInitPar(mod$par), rep(c(0.001, 0.1), 4)) # initial parameter values for the model parameters. Then, the initial values for the structural error parameters (gamma) for each observation
logpost <- logPosterior_BaM(
    parvector = parvector, # parameter values at which the posterior is evaluated
    X = X, Yobs = Yobs, Yu = Yu, # calibration data
    lpfunk = logPrior, llfunk = logLikelihood, mod = mod
) # inference functions and model

# The logPosterior function returns a list containing log-post, log-prior, log-lkh and simulated values
logpost

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Example of maximizing the log-posterior ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Just use function optim() coming with R base installation
maxpost <- stats::optim(
    par = parvector, # starting point
    fn = logPosterior_BaM_wrapped, # function to be optimized. It's a wrapped version of function logPosterior, returning the log-post value only. See ?logPosterior_wrapped
    control = list(fnscale = -1), # tells optim to maximize rather than minimize
    X = X, Yobs = Yobs, Yu = Yu, lpfunk = logPrior, llfunk = logLikelihood, mod = mod
) # arguments passed to logPosterior_wrapped
maxpost$par
maxpost$convergence

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Example of MCMC-sampling the log-posterior ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
x0 <- maxpost$par
names(x0) <- c(
    RBaM::getNames(mod$par),
    "g1_Y1", "g2_Y1", "g1_Y2", "g2_Y2",
    "g1_Y3", "g2_Y3", "g1_Y4", "g2_Y4"
)
# Adaptive Metropolis
ptm <- proc.time()
mcmc <- MCMC_AM(
    logPdf = logPosterior_BaM, x0 = x0,
    X = X, Yobs = Yobs, Yu = Yu, lpfunk = logPrior, llfunk = logLikelihood, mod = mod
) # arguments passed to logPosterior
proc.time() - ptm
pairs(cbind(mcmc$samples, mcmc$components$logPosterior))









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
