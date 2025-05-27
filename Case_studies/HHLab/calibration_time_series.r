cat("\014")
rm(list = ls())

library(RBaM)
library(dplyr)
library(stringr)
library(tidyr)

workspace <- "/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/Case_studies/HHLab/"
setwd(workspace)

path_results <- "Calibration_time_series"
path_model_mage_global <- "model_mage"
path_data <- "data_smooth_bed"

# General setting options:
# Filter observation data
nSlim_ZData <- 1

# BaM! model
nCycles <- 10 # Number of cycles
burn <- 0.5 # Percentage of data burned
nSlim <- 5 # Slim factor: after burning, only one iteration each Nslim is kept.

# Define Mage executable path, hard to define because depending on the machine. Avoid blank spaces!
MAGE_executable <- "/home/famendezrios/Documents/Softwares/pamhyr2/mage8/mage"
Mage_extraire_executable <- "/home/famendezrios/Documents/Softwares/pamhyr2/mage8/mage_extraire"
############################
# Module_Specific_Case
############################

# Input_CaseStudy
FilesNames <- list.files(path_model_mage_global, recursive = T)
dir_REPFile <- FilesNames[grep(FilesNames, pattern = ".REP")]
# Define the fixed-width format for each line : Result of Fortran coding format
col_widths <- c(1, 3, 6, 10, 10, 10, 10)

# Load calibration data:
load(file.path(
    "/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/processed_data/HHLab/smooth_bed/",
    "data_HHLab_all_cases.RData"
))

# User input:
## Give the information to use during calibration
all_data_calibration <- list(WS_profiles = WS_profiles)

u_h_obs <- 0.2 / 1000 # ultrasonic sensor manufactured by Baumer (UNDK 20IG903/S35A), with a standard measurement error of approximately 0.2 mm.


# Total number of discretization points required for Legendre polynomial calculations
total_points <- 200
doInterpolation <- TRUE
# Create a sequence of numbers
sequence <- seq(0.504, 16.504, by = 0.100)

sequence_all <- c(
    # sequence,
    0.060,
    0.650,
    1.650,
    3.650,
    5.650,
    7.650,
    9.650,
    11.650,
    15.650
)
# Find the maximum number of decimal places
max_decimals <- max(sapply(sequence_all, function(x) nchar(sub(".*\\.", "", as.character(x)))))

sequence_all <- sort(round(sequence_all, max_decimals))

# Prefix each number with "ZdT 1 "
mage_extraire_args <- lapply(sequence_all, function(x) paste("ZdT 1", format(x, nsmall = max_decimals)))

# Time to introduce calibration data (in hours)
# Time mapping between cases with constant discharge (id_case) and simulation time (id_fixed) for extraction, all consistent with the MAGE model!!
time_mapping_temp <- data.frame(
    id_case = c(
        "case_60_60",
        "case_30_30",
        "case_15_15",
        "case_07_07"
    ),
    id_fixed = c(
        0.95,
        1.95,
        2.95,
        3.95
    )
)

# Assign the number of output for WS profiles with the time fixed
# Variable extraction
first_char <- substr(unlist(mage_extraire_args), 1, 3)

# Position or time extraction
last_number <- stringr::str_extract(unlist(mage_extraire_args), "\\d+\\.\\d+$")

extraction_data_mage_extraire_args <- paste0(first_char, "_", last_number)

# Repeat and bind rows with extraction labels
time_mapping <- do.call(rbind, lapply(extraction_data_mage_extraire_args, function(val) {
    df <- time_mapping_temp
    df$id_position <- val
    df
}))

# Number of output variables : defined in mage_extraire_args
nY <- length(mage_extraire_args)

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

# Declaration of the variable
Data_Measured <- setNames(vector("list", length(all_data_calibration)), names(all_data_calibration))
# Check if WS_profiles are given for calibration
logical_WS_Calibration <- names(all_data_calibration) %in% "WS_profiles"
# Treatment for WS_profiles
if (logical_WS_Calibration) {
    id_WS_calibration <- which(names(all_data_calibration) == "WS_profiles")

    WS_profiles_calData <- all_data_calibration[[id_WS_calibration]]
    for (i in 1:length(WS_profiles_calData)) {
        id_case <- names(WS_profiles_calData)[i]
        WS_profile <- WS_profiles_calData[[i]]

        Data_Measured_temp <- cbind(WS_profile, id_case)

        Data_Measured_temp <- Data_Measured_temp %>%
            arrange(Data_Measured_temp$x)

        Data_Measured[[id_WS_calibration]] <- rbind(Data_Measured[[id_WS_calibration]], Data_Measured_temp)
    }
    ######
    # Get the "X" or 'PK' of the model defined in mage_extraire_args (
    # case WSL, need to be adapted for time series!)
    data_obs_calibration_temp <- Data_Measured[[id_WS_calibration]] %>%
        filter(Data_Measured[[id_WS_calibration]]$x %in% specific_points_CalStations)

    # # Find the maximum number of decimal places
    # max_decimals <- max(sapply(data_obs_calibration_temp$x, function(x) nchar(sub(".*\\.", "", as.character(x)))))

    # Fusionner avec les données originales
    data_obs_calibration_temp <- data_obs_calibration_temp %>%
        mutate(id_position = paste0(first_char[1], "_", trimws(format(x, nsmall = max_decimals, scientific = FALSE))))

    data_obs_calibration <- time_mapping %>%
        left_join(data_obs_calibration_temp, by = c("id_case", "id_position"))

    ####################################
    # Assign data
    ####################################
    # 1. Prepare `data_obs_calibration` for merging
    data_obs_calibration_merged_wide_all <- data.frame(data_obs_calibration %>%
        select(id_fixed, z_mean, id_position, Yu) %>%
        pivot_wider(names_from = id_position, values_from = c(z_mean, Yu)))

    data_obs_calibration_merged_wide <- data_obs_calibration_merged_wide_all[, c(1:(nY + 1))]

    clean_names_data <- sub("^z_mean_", "", colnames(data_obs_calibration_merged_wide)[-1])

    colnames(data_obs_calibration_merged_wide)[-1] <- clean_names_data

    u_data_obs_calibration_merged_wide <- data_obs_calibration_merged_wide_all[, c(1, (nY + 2):ncol(data_obs_calibration_merged_wide_all))]

    clean_names_u_data <- sub("^Yu_", "", colnames(data_obs_calibration_merged_wide)[-1])

    colnames(u_data_obs_calibration_merged_wide)[-1] <- clean_names_u_data
    # 2. Convert `CalData` from wide to long format
    CalData_temp <- data.frame(matrix(NA, nrow = nrow(X), ncol = nY))

    # Get column names from mage_extraire_args
    column_names_mage_extraire_args <- sapply(mage_extraire_args, function(arg) {
        # Get the first letter (e.g., "Z" or "Q")
        prefix <- substr(arg, 1, 3)

        # Split by spaces and take the last element (the number)
        split_text <- strsplit(arg, " ")[[1]]
        # value <- split_text[length(split_text)]
        value <- stringr::str_extract(arg, "\\d+\\.\\d+$")
        # Combine prefix and value to create column name
        paste0(prefix, "_", value)
    })

    colnames(CalData_temp) <- column_names_mage_extraire_args
    CalData <- cbind(X, CalData_temp)

    id_column <- which(colnames(data_obs_calibration_merged_wide) %in% colnames(CalData)[-1])
    id_replace <- which(CalData[, 1] %in% data_obs_calibration_merged_wide[, 1])
    CalData[id_replace, id_column] <- data_obs_calibration_merged_wide[, -1]


    # Calibration data which will be given to BaM! for calibration
    Y <- CalData[, -1]
    if (is.vector(Y)) {
        Y <- data.frame(Y)
        colnames(Y) <- colnames(CalData)[2]
    }

    Yu <- data.frame(matrix(NA, nrow = nrow(X), ncol = nY))
    colnames(Yu) <- colnames(Y)

    id_column <- which(colnames(u_data_obs_calibration_merged_wide)[-1] %in% colnames(CalData)[-1])
    id_replace <- which(CalData[, 1] %in% u_data_obs_calibration_merged_wide[, 1])
    Yu[id_replace, id_column] <- u_data_obs_calibration_merged_wide[, -1]
    colnames(Yu) <- paste0("Yu_", colnames(Y))
}
####


#-----------------------------------------------------------#
#                            RUG FILE :                     #
#-----------------------------------------------------------#
# Read RUG file
read_fortran_data <- function(file_path, col_widths, skip = 0) {
    # Read the file with the fixed-width format
    data <- read.fwf(file_path, widths = col_widths, header = FALSE, skip = skip)
    return(data)
}

Read_RUGFile <- read_fortran_data(
    file_path = file.path(
        path_model_mage_global,
        FilesNames[grep(FilesNames, pattern = ".RUG")]
    ),
    col_widths = col_widths,
    skip = 1
)

# Rename columns for better clarity
colnames(Read_RUGFile) <- c("K", "Reach", "Blank", "KP_start", "KP_end", "Kmin", "Kmoy")

# Drop the "Blank" column if not needed
RUGFile <- Read_RUGFile[, -3]
Nb_reaches <- unique(RUGFile$Reach)
# Be careful with RUG file format (talk to Theophile) to adapt the code
specific_points_Model <- data.frame()
for (i in 1:length(Nb_reaches)) {
    reach_temp <- RUGFile[which(RUGFile$Reach == Nb_reaches[i]), ]
    if (nrow(reach_temp) == 1) {
        specific_points_Model <- rbind(
            specific_points_Model,
            data.frame(
                reach = reach_temp$Reach,
                KP_start = reach_temp$KP_start,
                KP_end = reach_temp$KP_end
            )
        )
    } else {
        specific_points_Model <- rbind(
            specific_points_Model,
            data.frame(
                reach = reach_temp$Reach[1],
                KP_start = reach_temp$KP_start[1],
                KP_end = reach_temp$KP_end[nrow(reach_temp)]
            )
        )
    }
}

specific_points <- unique(sort(c(
    specific_points_Model$KP_start,
    specific_points_Model$KP_end,
    specific_points_CalStations
)))
##### Return module_mage

all_specific_points <- specific_points
specific_points_CalStations <- sort(specific_points_CalStations)

#######################################
# Come back to the specific case module
#######################################


# Get values interpolated respecting defined values: these values will be given to estimate Legendre polynomials
#-----------------------------#
# !!!!!!! To improve : when reaches are not consecutive, it must be handled differently. Because a interpolation must be done reach by reach
#-----------------------------#

# Function to interpolate specifying some points
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

    if (!all(all_specific_points %in% specific_seq)) stop("Distance so close, increase number of interpolate KP to get all
                                                    points required")

    if (length(specific_seq) != total_points) stop("Something is wrong with the interpolation")
    return(specific_seq)
}

if (doInterpolation) {
    grid_covariant_discretized <- interpolation_specific_points(
        total_points = total_points,
        all_specific_points = all_specific_points
    )
} else {
    grid_covariant_discretized <- all_specific_points
}


# Function to get the boundaries points from the middle point.
# After using a interpolation, the values are at the middle of two points.
# In this case, we have the interpolated values, we are searching the boundaries points
get_covariant_meshed <- function(covariant_discretized, specific_points_Model) {
    # Reconstruction of covariant_meshed_temps :
    covariant_meshed_temps <- numeric()

    # 1. Find the limit points for the first two data, then for the rest of the data
    covariant_meshed_temps[1] <- sum(covariant_discretized[c(1, 2)]) / 2

    for (i in 2:(length(covariant_discretized) - 1)) {
        # Loop starts at 2 because first data is already given
        # Loop ends at n-1 number of covariant_discretized because
        # to find limits points, at least two points are required
        covariant_meshed_temps[i] <- sum(covariant_discretized[c(i + 1, i)]) / 2
    }
    # Add model's specific points to the mesh. Necessary condition for running Mage model
    covariant_meshed <- unique(sort(c(
        specific_points_Model$KP_start[1],
        specific_points_Model$KP_end,
        covariant_meshed_temps
    )))
    return(covariant_meshed)
}

# Calculate original covariant considering that grid_covariant_discretized is in the middle of the real covariant value except on the boundaries

grid_covariant_meshed_Model <- get_covariant_meshed(
    covariant_discretized = grid_covariant_discretized,
    specific_points_Model = specific_points_Model
)


##############################################
#####  Polynomial de legendre
##############################################

###############################
# n iteration : methodology
###############################
n_degree_max <- 4
n_degree_seq <- seq(0, n_degree_max, 1)

dir_exe_BaM <- file.path(find.package("RBaM"), "bin")
dir_cf <- file.path(dir_exe_BaM, "Config_BaM.txt")

for (n_degree in n_degree_seq) {
    # Polynomial degree i
    path_polynomial <- file.path(path_results, paste0("n_", n_degree))
    if (!dir.exists(path_polynomial)) {
        dir.create(path_polynomial)
    }

    # Get model mage files
    files_model_mage <- list.files(path_model_mage_global, full.names = TRUE, recursive = TRUE)
    # Copy files for modify at each polynomial degree
    # Loop through each file path
    for (i in seq_along(files_model_mage)) {
        from <- files_model_mage[i]
        to <- file.path(path_polynomial, files_model_mage[i])

        # Create destination directory if needed
        dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)

        # Copy file
        success <- file.copy(from = from, to = to, overwrite = TRUE)
    }

    path_model_mage <- file.path(path_polynomial, path_model_mage_global)

    # dataset object
    data <- dataset(X = X, Y = Y, Yu = Yu, data.dir = file.path(workspace, path_polynomial))
    # #-----------------------------#
    # Write RUG file:
    RUG_id_reach <- vector(mode = "numeric", length = (length(grid_covariant_meshed_Model) - 1))

    for (i in 1:nrow(specific_points_Model)) {
        RUG_id_reach[dplyr::between(
            grid_covariant_meshed_Model[-1], # Remove the first data, because RUG file is given by intervals, so n-1 number of intervals where n is the number of points of positions
            specific_points_Model$KP_start[i], specific_points_Model$KP_end[i]
        )] <- specific_points_Model$reach[i]
    }

    RUG_KP_start <- grid_covariant_meshed_Model[-length(grid_covariant_meshed_Model)]
    RUG_KP_end <- grid_covariant_meshed_Model[-1]
    ######################################################
    # Estimation of the prior on the friction coefficient
    ######################################################
    # ref : https://www.hec.usace.army.mil/confluence/rasdocs/ras1dtechref/6.1/modeling-culverts/culvert-data-and-coefficients/manning-s-roughness-coefficient

    Ks_glass_quantiles <- c(1 / 0.013, 1 / 0.009)

    # Estimate standard deviation
    prior_log_u_ks_mage <- (log(Ks_glass_quantiles[2]) - log(Ks_glass_quantiles[1])) / (2 * qnorm(0.975))

    # Estimate mean
    prior_ks_mage <- 1 / 0.010

    RUG_Kmin <- prior_ks_mage
    RUG_Kmoy <- 10

    RUG_path <- file.path(
        path_model_mage,
        FilesNames[grep(FilesNames, pattern = ".RUG")]
    )
    # Function to write RUG file
    write_RUGFile <- function(RUG_path,
                              RUG_id_reach,
                              RUG_KP_start,
                              RUG_KP_end,
                              RUG_Kmin,
                              RUG_Kmoy,
                              RUG_format) {
        # Open a .RUG file for writing
        fileConn <- file(RUG_path, "w")
        # Write the first line as a comment
        writeLines("* This file is generated by PAMHYR, please don't modify", fileConn)
        # Loop to write each line in Fortran format
        for (i in 1:length(RUG_KP_start)) {
            # Format the line according to Fortran style
            line <- sprintf(
                "%1s%3d      %10.3f%10.3f%10.2f%10.2f",
                "K",
                RUG_id_reach[i],
                RUG_KP_start[i],
                RUG_KP_end[i],
                RUG_Kmin,
                RUG_Kmoy
            )
            # Write to file
            writeLines(line, fileConn)
        }
        # Close the file
        close(fileConn)
    }

    write_RUGFile(
        RUG_path = RUG_path,
        RUG_id_reach = RUG_id_reach,
        RUG_KP_start = RUG_KP_start,
        RUG_KP_end = RUG_KP_end,
        RUG_Kmin = RUG_Kmin,
        RUG_Kmoy = RUG_Kmoy,
        RUG_format = "%1s%3d      %10.0f%10.0f%10.2f%10.2f"
    )

    # Polynomial degree 0 :
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
        file = file.path(path_polynomial, "vector_legendre.txt"), row.names = F
    )
    # Write co variant matrix for main channel
    # Create an empty data frame and fill it with polynomial values
    legendre_df <- data.frame(x = cov_bar_legendre) # Start with x values
    ########
    # Generate files
    ########

    getlegendre <- function(degree, normalized_values) {
        if (degree < 0) stop("Degree must be positive")
        if (any(!dplyr::between(normalized_values, -1, 1))) stop("Range of normaized_values should be between [-1,1]")

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

    zFileKmin <- file.path(getwd(), path_polynomial, "Zfile_Kmin.txt")
    zFileKmoy <- file.path(getwd(), path_polynomial, "Zfile_Kflood.txt")

    if (ncol(legendre_df) == 2) { # if n degree is equal to 0
        matrix_zFileKmin <- legendre_df["P0"]
    } else {
        matrix_zFileKmin <- legendre_df[, -1]
    }
    matrix_zFileKmoy <- legendre_df["P0"]

    #######################
    # Fin legendre
    ########################

    ##########################
    #### BaM!
    ###########################
    param_theta <- list(
        # Kmin.init = log(130),
        Kmin.init = prior_ks_mage,
        kmin.distri = "FlatPrior+",
        Kmoy.init = 10,
        kmoy.distri = "FIX"
    )

    # hist(rlnorm(10000, log(prior_ks_mage), prior_log_u_ks_mage))



    param_error_model <- list(
        ## Stage model error prior information
        z.ini = 0.001,
        z.prior.dist = "FlatPrior"
        # z.prior.par = c(0, 1)
    )

    # prior_setting <- list(
    #     param_theta = param_theta,
    #     param_error_model = param_error_model
    # )

    # Input
    remant_error <- list(
        remnantErrorModel(
            fname = "Config_RemnantSigma.txt",
            funk = "Constant",
            par = list(parameter(
                name = "intercept",
                init = param_error_model$z.ini,
                prior.dist = param_error_model$z.prior.dist,
                prior.par = param_error_model$z.prior.par
            ))
        )
    )

    remant_error_list <- rep(remant_error, nY)

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

    ID <- "MAGE_TEMP"
    MageVersion <- "8"

    xtra <- xtraModelInfo(
        fname = "Config_setup.txt",
        object = list(
            exeFile = MAGE_executable,
            version = MageVersion,
            mageDir = file.path(getwd(), paste0(path_model_mage), ""),
            repFile = dir_REPFile,
            zKmin = matrix_zFileKmin,
            zFileKmin = zFileKmin,
            doExpKmin = FALSE,
            zKmoy = matrix_zFileKmoy,
            zFileKmoy = zFileKmoy,
            doExpKmoy = FALSE,
            mage_extraire_file = Mage_extraire_executable,
            mage_extraire_args = unlist(mage_extraire_args)
        )
    )
    mod <- model(
        ID = ID,
        nX = 1,
        nY = nY,
        par = theta_param,
        xtra = xtra
    )
    # Define MCMC options:
    # Main idea, replace default value for MCMC jump exploration
    # In fact when init prior of a parameter or a error model parameter is equal 0, jump is defined as 0.1 * (init prior), so jump is so small
    # To improve this, we need to pass to manual mode for assigning jumps individually for all parameters (theta and gamma)

    # Function to get prior of the parameters
    get_init_prior <- function(parameter) {
        # Identify if parameter is remnantErrorModel class
        logical_test <- class(parameter[[1]]) == "remnantErrorModel"
        if (logical_test) { # if remnantErrorModel, a list is needed
            init_priors <- list()
        } else { # if not a vector is needed
            init_priors <- numeric(0)
        }

        counter_gamma <- 1
        for (i in parameter) {
            # Handle if parameters is remnantErrorModel
            if (logical_test) {
                param <- i$par
                number_var_error_model <- seq_along(param)

                for (local_counter in number_var_error_model) {
                    if (param[[local_counter]]$prior$dist != "FIX") {
                        init_priors[[counter_gamma]] <- param[[local_counter]]$init
                        counter_gamma <- counter_gamma + 1
                    }
                }
            } else {
                # Handle if parameters comes from theta
                param <- i
                if (param$prior$dist != "FIX") {
                    init_priors <- c(init_priors, param$init)
                }
            }
        }
        return(init_priors)
    }

    prior_error_model <- get_init_prior(remant_error_list)
    prior_theta_param <- get_init_prior(theta_param)

    jump_MCMC_theta_param_user <- 10
    jump_MCMC_error_model_user <- 0.005

    jump_MCMC_theta_param <- ifelse(prior_theta_param != 0,
        prior_theta_param[(prior_theta_param != 0)] * 0.1,
        jump_MCMC_theta_param_user
    )

    jump_MCMC_error_model <- list()
    for (ind in 1:length(prior_error_model)) {
        jump_MCMC_error_model[[ind]] <- ifelse(prior_error_model[[ind]] > 0.5,
            prior_error_model[[ind]][(prior_error_model[[ind]] != 0)] * 0.1,
            jump_MCMC_error_model_user
        )
    }


    mcmcOptions_user <- mcmcOptions(
        nCycles = nCycles,
        manualMode = TRUE,
        thetaStd = jump_MCMC_theta_param,
        gammaStd = jump_MCMC_error_model
    )
    BaM(
        mod = mod,
        data = data,
        remnant = remant_error_list,
        mcmc = mcmcOptions_user,
        cook = mcmcCooking(
            burn = burn,
            nSlim = nSlim
        ),
        summary = mcmcSummary(xtendedMCMC.fname = "Results_xtendedMCMC.txt"),
        residuals = residualOptions(),
        pred = NULL,
        doCalib = TRUE,
        doPred = FALSE,
        na.value = -9999,
        run = FALSE,
        preClean = FALSE,
        workspace = file.path(workspace, path_polynomial),
        dir.exe = dir_exe_BaM,
        name.exe = "BaM",
        predMaster_fname = "Config_Pred_Master.txt"
    )

    # Move Config_BaM.txt file to the result folder
    file.copy(from = dir_cf, to = path_polynomial, overwrite = TRUE)

    system2(
        command = file.path(dir_exe_BaM, "BaM"),
        args = c("-cf", file.path(workspace, path_results, paste0("n_", n_degree, "/Config_BaM.txt"))),
        wait = FALSE
    )
}

# # PREDICTION : not working!

# # Define the grid of stage values on which the rating curve will be computed
# Tgrid <- data.frame(time = seq(0, 4, 0.05))
# # Define a 'prediction' object for total predictive uncertainty
# totalU <- prediction(
#     X = Tgrid, # stage values
#     spagFiles = "totalU.spag", # file where predictions are saved
#     data.dir = file.path(getwd(), path_results), # a copy of data files will be saved here
#     doParametric = TRUE, # propagate parametric uncertainty, i.e. MCMC samples?
#     doStructural = TRUE
# ) # propagate structural uncertainty ?
# # Define a 'prediction' object for parametric uncertainty only - not the doStructural=FALSE
# paramU <- prediction(
#     X = Tgrid, spagFiles = "paramU.spag", data.dir = file.path(getwd(), path_results),
#     doParametric = TRUE, doStructural = FALSE
# )
# # Define a 'prediction' object with no uncertainty - this corresponds to the 'maxpost' rating curve maximizing the posterior
# maxpost <- prediction(
#     X = Tgrid, spagFiles = "maxpost.spag", data.dir = file.path(getwd(), path_results),
#     doParametric = FALSE, doStructural = FALSE
# )
# # Re-run BaM, but in prediction mode
# BaM(
#     mod = mod,
#     data = data,
#     remnant = rep(remant_error, mod$nY),
#     pred = list(totalU, paramU, maxpost), # list of predictions,
#     doCalib = FALSE,
#     doPred = TRUE,
#     workspace = file.path(getwd(), path_results)
# )
