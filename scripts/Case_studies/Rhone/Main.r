rm(list = ls())
graphics.off()

# Set directory (root of the repository)
dir_workspace <- "/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git"
setwd(dir_workspace)

function_list <- list.files("R", full.names = TRUE)
for (i in function_list) {
    source(i)
}

#############################################
# Load libraries
#############################################

library(RBaM)
library(dplyr)
library(patchwork)
library(tidyr)
library(ggplot2)
library(stringr)
library(lubridate)
library(utils)
#############################################
# End load libraries
#############################################

############################################
# Module 0 : General settings
############################################

# Logical
do_calibration <- FALSE
do_plot_calibration <- FALSE
do_prediction <- FALSE


# Setting jump standard deviation for MCMC sampling if initial guess = 0
jump_MCMC_theta_param_user_regression <- 5
jump_MCMC_theta_param_user_coeff <- 0.1
jump_MCMC_error_model_user <- 0.001
threshold_jump_MCMC_error_model <- 0.5

############################################
# Module 1 : set directory paths
############################################

# Name of the experiment. All scenarios will be used the same calibration data

# Common for observation and calibration folders: names of the experiment
Experiment_id <- c(
    "1_WSE_AIN_90_2_WSE_RHONE_525_750"
)


# Experiments input data to be used during calibration setting
all_cal_case <- c(
    "Kmin_n_0_Q0.r",
    "Kmin_n_4_Q0.r"
)

# Folder related to the observations (careful with the order!)

all_events <- c(
    "AIN_90",
    "RHONE_525",
    "RHONE_750"
)

command_line_MAGE <- "-fp=2 -LC=0 -eps=5"
file_main_path <- file.path(dir_workspace, "scripts/Case_studies/Rhone/Real_Condition_model/Calibration_experiments")
MAGE_main_folder <- "/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/scripts/Case_studies/Rhone/Real_Condition_model/model_mage"
mage_projet_name <- "Rhone_PCH_Ain"

check_mage_folder(name_folder = MAGE_main_folder)
############################################
# End module 1: set directories paths
############################################


############################################
# Module 2 : hydraulic model (HM) environment
############################################
# Input information about reaches and boundaries KP of the HM.
# Model_Reach : Model reach interpreted by the HM.
# This input information must be coherent with HM specifications.

# KP start and end must be reliable to the model. To define KP range, it does not matter the direction, it must be coherent to the HM

Input_Model_Reach <- data.frame(
    reach = c(1, 2, 3, 4, 5, 6, 7, 8),
    KP_start = c(55900, 36250, 34500, 37491, 41211, 0, 20, 22333),
    KP_end = c(36250, 34500, 26750, 41211, 41461, 20, 40, 37491)
)

# First layer: Typology
# ###########################################
# All reaches in the model (Model_Reach) defined within the same coordinate reference. Useful to separate main reach to tributary or even a diversion channel

# Must be careful with the order of the reaches, it must be given upstream to downstream
Input_Typology <- list(
    Rhone = c(1, 2, 3),
    Ain = c(8, 4, 5),
    Cassier = c(6, 7)
)

############################################
# End module 2 : hydraulic model (HM) environment
############################################


############################################
# Module 3: calibration data
############################################

# Processed data
if (Experiment_id %in% c("1_WSE_AIN_90_2_WSE_RHONE_525_750")) {
    load("data/processed_data/Rhone/Ain_90_RH_525_750/observed_data.RData")
} else {
    stop("Experiment id in not correct")
}


######################################## UNTIL HERE

X <- observed_data[, c(
    "event",
    "id_reach_CAL",
    "KP",
    "time"
)] %>%
    rename(
        "reach" = "id_reach_CAL",
        "x" = "KP",
        "t" = "time"
    )


results_CalData <- constructor_CalData(
    observed_data = observed_data,
    do_manual_uncertainty = FALSE,
    sd_WSE_fixed = 0.05, # in meters
    sd_Q_fixed = 8, # in %
    sd_V_fixed = 3, # in %
    sd_Kmin_fixed = 5, # in m1/3/s
    sd_Kflood_fixed = 5 # in m1/3/s
)

Y <- results_CalData$Y
Yu <- results_CalData$Yu

CalData <- cbind(X, Y, Yu)


##########################################
# Specificities : add pseudo obs of friction
# obs always has a gaussian distribution
Kmin_prior_u_Yu <- c(30, 10)
Kflood_prior_u_Yu <- c(25, 8)


#############################
# Add x,t coordinates to spatially distributed friction observations
#############################

# Specific of each case: here, i do not want assign pseudo obs at the reservoir
nodes_to_add_pseudo_obs <- unlist(c(Input_Typology$Rhone, Input_Typology$Ain))

# Distance between two consecutive pseudo obs
step_size <- 150

K_pseudo_obs <- data.frame(Input_Model_Reach %>%
    filter(reach %in% nodes_to_add_pseudo_obs) %>%
    rowwise() %>%
    mutate(KP = list(
        seq(KP_start,
            KP_end,
            length.out = ceiling(abs(KP_end - KP_start) / step_size) + 1
        )
    )) %>%
    unnest(KP) %>%
    select(reach, KP) %>%
    rename(x = KP) %>%
    mutate(
        event = 1,
        t = CalData$t[1]
    ) %>%
    select(event, reach, x, t) %>%
    mutate(
        WSE = -9999,
        Q = -9999,
        V = -9999,
        Kmin = Kmin_prior_u_Yu[1],
        Kflood = Kflood_prior_u_Yu[1],
        Yu_WSE = -9999,
        Yu_Q = -9999,
        Yu_V = -9999,
        Yu_Kmin = Kmin_prior_u_Yu[2],
        Yu_Kflood = Kflood_prior_u_Yu[2]
    ))


CalData <- rbind(CalData, K_pseudo_obs)

path_experiment <- file.path(file_main_path, Experiment_id)

if (do_plot_calibration) {
    CalData_ordered <- CalData %>%
        mutate(reach = factor(reach,
            levels = c(1, 2, 3, 8, 4, 5)
        ))

    plots_CalData <- plot_CalData(
        CalData = CalData_ordered,
        scales_free = "free",
        wrap = TRUE
    )
    # Customize the plot
    plots_CalData$plot_WSE <- plots_CalData$plot_WSE +
        facet_wrap(
            ~ event + reach,
            labeller = labeller(
                event = as_labeller(
                    c(
                        "1" = "Event~1:~Q(Ain):~90~m^3/s",
                        "2" = "Event~2:~Q(Rhone):~525~m^3/s",
                        "3" = "Event~3:~Q(Rhone):~750~m^3/s"
                    ),
                    label_parsed
                ),
                reach = as_labeller(c(
                    "1" = "Reach~1:~Lagnieu~-~Bourbe",
                    "2" = "Reach~2:~Bourbe~-~Anthon",
                    "3" = "Reach~3:~Anthon~-~Jons",
                    "4" = "Reach~4:~Port~Galland~-~reservoir",
                    "5" = "Reach~5:~reservoir~-~Anthon",
                    "8" = "Reach~8:~Pont~Chazey~-~Port~Galland"
                ), label_parsed)
            ),
            scales = "free",
            ncol = 3
        )

    plots_CalData$plot_Kmin <- plots_CalData$plot_Kmin +
        facet_wrap(
            ~reach,
            labeller = labeller(
                reach = as_labeller(c(
                    "1" = "Reach~1:~Lagnieu~-~Bourbe",
                    "2" = "Reach~2:~Bourbe~-~Anthon",
                    "3" = "Reach~3:~Anthon~-~Jons",
                    "4" = "Reach~4:~Port~Galland~-~reservoir",
                    "5" = "Reach~5:~reservoir~-~Anthon",
                    "8" = "Reach~8:~Pont~Chazey~-~Port~Galland"
                ), label_parsed)
            ),
            scales = "free",
            ncol = 3
        )

    plots_CalData$plot_Kflood <- plots_CalData$plot_Kflood +
        facet_wrap(
            ~reach,
            labeller = labeller(
                reach = as_labeller(c(
                    "1" = "Reach~1:~Lagnieu~-~Bourbe",
                    "2" = "Reach~2:~Bourbe~-~Anthon",
                    "3" = "Reach~3:~Anthon~-~Jons",
                    "4" = "Reach~4:~Port~Galland~-~reservoir",
                    "5" = "Reach~5:~reservoir~-~Anthon",
                    "8" = "Reach~8:~Pont~Chazey~-~Port~Galland"
                ), label_parsed)
            ),
            scales = "free",
            ncol = 3
        )

    obs_adapted <- observed_data %>%
        rename("reach" = "id_reach_CAL")

    plot_WSE_Thalweg <- plots_CalData$plot_WSE +
        geom_line(data = obs_adapted, aes(x = KP, y = Z_thalweg), color = "black")

    plots_CalData$plot_WSE_Thalweg <- plot_WSE_Thalweg

    if (!dir.exists(path_experiment)) {
        dir.create(path_experiment)
    }

    for (name in names(plots_CalData)) {
        plot <- plots_CalData[[name]]
        if (!is.null(plot)) {
            ggsave(file.path(path_experiment, paste0(name, ".png")),
                plot = plot,
                width = 35,
                height = 20,
                units = "cm"
            )

            save(
                file = file.path(path_experiment, paste0(name, ".RData")),
                plot
            )
        }
    }
}
############################################
# End module 2: calibration data
############################################

############################################
# Module 4: calibration setting
############################################

######################################################
# Set Remnant error model
# Structural error information
# Put remnantErrorModel_default on the output variable without calibration data
remant_error_list <- list(
    # WSE
    RBaM::remnantErrorModel(
        fname = "Config_RemnantSigma.txt",
        funk = "Constant",
        par = list(parameter(
            name = "intercept",
            init = 0.0005,
            prior.dist = "FlatPrior"
        ))
    ),
    # Q
    remnantErrorModel_default(name = "Config_RemnantSigma2.txt"),
    # V
    remnantErrorModel_default(name = "Config_RemnantSigma3.txt"),
    # Both Kmin and Kflood, they must always be forced to 0. Information will be passed by pseudo-obs
    # Kmin
    remnantErrorModel_default(name = "Config_RemnantSigma4.txt"),
    # Kflood
    remnantErrorModel_default(name = "Config_RemnantSigma5.txt")
)

# Get initial prior of structural error
prior_error_model <- get_init_prior(remant_error_list)

############################################
# End module 4: calibration setting
############################################



############################################
# Module 5: BaM environment
############################################

Key_Info_Typology_Model_Reach <- get_Key_Info_Typology_Model_Reach(
    Input_Typology = Input_Typology,
    Input_Model_Reach = Input_Model_Reach,
    total_points_discretization = 100
)

############################################
# Module 6: Calibration
############################################
list_mod_polynomials <- list_Z_MatrixKmin <- list_Z_MatrixKflood <- list_Kmin_prior <- list_Kflood_prior <- list_Kmin_SU <- list_Kflood_SU <- list_summary_SU_Kflood <- list_summary_SU_Kmin <- list()

for (id_cal_case in 1:length(all_cal_case)) {
    # Load experiment
    paths <- load_experiment(
        file_main_path = file_main_path,
        cal_case = all_cal_case[[id_cal_case]],
        path_experiment = path_experiment,
        all_events = all_events
    )

    results_estimation <- Estimation_Mage(
        paths = paths,
        Key_Info_Typology_Model_Reach = Key_Info_Typology_Model_Reach,
        Input_Typology = Input_Typology,
        MAGE_main_folder = MAGE_main_folder,
        do_calibration = do_calibration,
        command_line_MAGE = command_line_MAGE,
        ID_model_BaM = ID_model_BaM,
        nX_BaM = nX_BaM,
        nY_BaM = nY_BaM,
        mage_projet_name = mage_projet_name,
        mcmcCooking = RBaM::mcmcCooking(burn = 0, nSlim = 1),
        mcmcOptions = RBaM::mcmcOptions(nAdapt = 20, nCycles = 10),
        mcmcSummary = RBaM::mcmcSummary(xtendedMCMC.fname = "Results_xtendedMCMC.txt"),
        remant_error_list = remant_error_list
    )
    if (do_calibration) {
        script_path <- file.path(paths$path_BaM_folder, "run_BaM.sh")

        writeLines(
            c(
                "#!/bin/bash",
                paste(
                    shQuote(file.path(RBaM::getPathToBaM(), "BaM")),
                    "-cf",
                    shQuote(file.path(paths$path_BaM_folder, "Config_BaM.txt"))
                )
            ),
            script_path
        )

        Sys.chmod(script_path, "0755")
        # Run outside of Vscodium. To kill a job : pkill -f BaM
        system2(
            "nohup",
            args = c("bash", script_path),
            stdout = FALSE,
            stderr = FALSE,
            wait = FALSE
        )
    }
    list_Z_MatrixKmin[[id_cal_case]] <- results_estimation$Z_MatrixKmin
    list_Z_MatrixKflood[[id_cal_case]] <- results_estimation$Z_MatrixKflood
    list_Kmin_prior[[id_cal_case]] <- results_estimation$Kmin_prior
    list_Kflood_prior[[id_cal_case]] <- results_estimation$Kflood_prior
    list_Kmin_SU[[id_cal_case]] <- results_estimation$Kmin_SU
    list_Kflood_SU[[id_cal_case]] <- results_estimation$Kflood_SU
    list_mod_polynomials[[id_cal_case]] <- results_estimation$mod
    list_summary_SU_Kmin[[id_cal_case]] <- results_estimation$summary_SU_Kmin
    list_summary_SU_Kflood[[id_cal_case]] <- results_estimation$summary_SU_Kflood
}






Kmin_segment_layer <- segment_layer_reference(K_literature = Kmin_literature)
Kflood_segment_layer <- segment_layer_reference(K_literature = Kflood_literature)


#######################
plotDIC <- plot_DIC(path_experiment)
ggsave(
    file.path(
        path_experiment,
        "DIC.png"
    ),
    plotDIC,
    width = 20,
    height = 20,
    units = "cm"
)

# At each eid_cal_case
Q_observed <- c(162, 162)
final_calibraion <- FALSE
counter <- 1
for (id_cal_case in 1:length(all_cal_case)) {
    # Source the input data for the experiments
    path_Experiment_Input_Data <- file.path(file_main_path, "Experiments_Input_Data", all_cal_case[id_cal_case])
    if (!file.exists(path_Experiment_Input_Data)) {
        stop(paste0(
            "The experiment input data named : '",
            all_cal_case[id_cal_case],
            "' does not exist"
        ))
    }

    source(path_Experiment_Input_Data)

    path_polynomial <- file.path(
        path_experiment,
        sub("\\.r$", "", all_cal_case[id_cal_case])
    )
    path_temp_plots <- file.path(path_polynomial, "BaM")
    path_post_traitement <- file.path(path_temp_plots, "post_traitement")
    path_post_traitement_data <- file.path(path_post_traitement, "RData")

    if (!dir.exists(path_post_traitement)) {
        dir.create(path_post_traitement)
    }

    if (!dir.exists(path_post_traitement_data)) {
        dir.create(path_post_traitement_data)
    }

    path_model_mage <- c(paste0(path_polynomial, "/model_mage/", all_events, "/"))

    mcmc_not_cooked <- readMCMC(file.path(path_temp_plots, "Results_MCMC.txt"))
    plots <- tracePlot(mcmc_not_cooked)

    mcmcplot <- wrap_plots(plots, ncol = 3)
    # Save plot and data
    ggsave(
        file.path(
            path_post_traitement,
            "MCMC_not_cooked.png"
        ),
        mcmcplot,
        width = 20,
        height = 20,
        units = "cm"
    )

    # Density plot for each parameter
    plots <- densityPlot(mcmc_not_cooked)
    pdf_plot <- wrap_plots(plots, ncol = 3)

    ggsave(
        file.path(
            path_post_traitement,
            "densityplot_not_cooked.png"
        ),
        pdf_plot,
        width = 20,
        height = 20,
        units = "cm"
    )
    if (final_calibraion) {
        if (!file.exists(file.path(path_temp_plots, "Results_Cooking.txt"))) stop("MCMC is still running or calculation is not going to the end. Verify if calibration is already finished or verify that calibration has not error messages. Please put final_results = FALSE as input")

        mcmc <- readMCMC(file.path(path_temp_plots, "Results_Cooking.txt"))
        plots <- tracePlot(mcmc)

        mcmcplot <- wrap_plots(plots, ncol = 3)


        ggsave(
            file.path(
                path_post_traitement,
                "MCMC_Cooked.png"
            ),
            mcmcplot,
            width = 20,
            height = 20,
            units = "cm"
        )

        # Density plot for each parameter
        plots <- densityPlot(mcmc)
        pdf_plot <- wrap_plots(plots, ncol = 3)

        ggsave(
            file.path(
                path_post_traitement,
                "densityplot_Cooked.png"
            ),
            pdf_plot,
            width = 20,
            height = 20,
            units = "cm"
        )


        png(
            file.path(path_post_traitement, "corelation_cooked.png"),
            width = 800,
            height = 800,
            res = 120
        )
        pairs(mcmc)
        dev.off()

        getSummary <- read.table(
            file = file.path(path_temp_plots, "Results_Summary.txt"),
            header = TRUE,
            stringsAsFactors = FALSE
        )

        # Values of error model in meter for WSE, in m3/s for discharge and m/s for velocity.
        knitr::kable(getSummary,
            align = "c"
        )


        # Zoom into the MAP and standard deviation of the error model
        getSummary_zoom <- getSummary[c(11, 16), ]
        getSummary_zoom[, c("Y1_intercept")] <- getSummary_zoom[, c("Y1_intercept")] # WSE in m

        # Get MAP simulation
        MAP_param_matrix <- as.numeric(getSummary_zoom[2, c(1:(length(Kmin_prior) + length(Kflood_prior)))])
    } else {
        results_MCMC_sampling <- read.table(file.path(path_temp_plots, "Results_MCMC.txt"), header = TRUE)
        mcmc <- data.frame(results_MCMC_sampling[, 1:(length(Kmin_prior) + length(Kflood_prior))])

        # Get MAP simulation: from current analysis without burning and slim
        MAP_param_matrix <- as.numeric(results_MCMC_sampling[which.max(results_MCMC_sampling$LogPost), 1:(length(Kmin_prior) + length(Kflood_prior))])
    }

    SR_reaches <- do.call(
        rbind,
        lapply(names(SR_reaches), function(river) {
            data.frame(
                value = SR_reaches[[river]]$SR1,
                id_river = river
            )
        })
    )


    Kmin_calculations <- K_plot(
        matrix_spatialisation = Z_MatrixKmin,
        mcmc = mcmc,
        n_param_Kmin = length(Kmin_prior),
        n_param_Kflood = length(Kflood_prior),
        covariate_discretization = covariate_grid,
        SR_reaches = SR_reaches,
        MAP_param_vector = MAP_param_matrix,
        main_channel = TRUE
    )
    df_MAP_Kmin <- Kmin_calculations[[1]]
    Kmin_plot <- Kmin_calculations[[2]]
    df_envelope_Kmin <- Kmin_calculations[[3]]

    ls_spatial_friction_Kmin <- list(
        df_envelope = df_envelope_Kmin,
        df_MAP = df_MAP_Kmin
    )
    save(ls_spatial_friction_Kmin,
        file = file.path(path_post_traitement_data, "Data_friction_estimation_ls_spatial_friction_Kmin.RData")
    )

    final_Kmin_plot <- Kmin_plot + Kmin_segment_layer
    scale_linetype_manual(
        name = "Reference\nvalues",
        values = c("Chow (1959)" = "dashed")
    )

    ggsave(
        file.path(
            path_post_traitement,
            "Kmin.png"
        ),
        final_Kmin_plot,
        width = 20,
        height = 20,
        units = "cm"
    )
    save(final_Kmin_plot,
        file = file.path(path_post_traitement_data, "Plot_friction_estimation_plot_Kmin_plot.RData")
    )

    Kflood_calculations <- K_plot(
        matrix_spatialisation = Z_MatrixKflood,
        mcmc = mcmc,
        n_param_Kmin = length(Kmin_prior),
        n_param_Kflood = length(Kflood_prior),
        covariate_discretization = covariate_grid,
        MAP_param_vector = MAP_param_matrix,
        main_channel = FALSE,
    )
    df_MAP_Kflood <- Kflood_calculations[[1]]
    Kflood_plot <- Kflood_calculations[[2]]
    df_envelope_Kflood <- Kflood_calculations[[3]]

    ls_spatial_friction_Kflood <- list(
        df_envelope = df_envelope_Kflood,
        df_MAP = df_MAP_Kflood
    )
    save(ls_spatial_friction_Kflood,
        file = file.path(path_post_traitement_data, "Data_friction_estimation_ls_spatial_friction_Kflood.RData")
    )

    final_Kflood_plot <- Kflood_plot + Kflood_segment_layer +
        scale_linetype_manual(
            name = "Reference\nvalues",
            values = c("Chow (1959)" = "dashed")
        )

    ggsave(
        file.path(
            path_post_traitement,
            "Kflood.png"
        ),
        final_Kflood_plot,
        width = 20,
        height = 20,
        units = "cm"
    )
    save(final_Kflood_plot,
        file = file.path(path_post_traitement_data, "Plot_friction_estimation_plot_Kflood_plot.RData")
    )
    ####################################
    # Residuals
    ###################################
    CalData <- read.table(file.path(path_temp_plots, "CalibrationData.txt"), header = TRUE)
    # Read residuals
    if (final_calibraion) {
        residuals <- read.table(
            file = file.path(path_temp_plots, "Results_Residuals.txt"),
            header = TRUE,
            stringsAsFactors = FALSE
        )

        # Convert residuals to mm and m3/s
        residuals_mm_m3_s <- data.frame(
            residuals[, 1:4],
            residuals[, c(9, 19, 20, 24, 29)]
        )
        sim_event_mm_m3_s <- data.frame(residuals_mm_m3_s, Yu_z = CalData$Yu_z)
    } else {
        # Residuals file is not ready, so I need to create by myself with MAP estimation from sampled data while MCMC is turning
        # Get data format from mage results
        temporal_dir <- tempdir()
        files <- list.files(temporal_dir, full.names = TRUE)

        exclude_file <- file.path(temporal_dir, "vscode-R")

        # List all files in the temp directory
        files <- list.files(temporal_dir, full.names = TRUE)

        # Exclude the specific file
        files_to_delete <- setdiff(files, exclude_file)

        # Remove remaining files if any
        if (length(files_to_delete) > 0) {
            unlink(files_to_delete, recursive = TRUE)
        }
        file.copy(
            from = path_model_mage,
            to = temporal_dir, recursive = TRUE
        )
        temp_path <- file.path(temporal_dir, basename(path_model_mage))

        path_RUGFile <- file.path(
            temp_path,
            list.files(temp_path)[grep(list.files(temp_path), pattern = ".RUG")]
        )
        MAP_RUGFile <- lapply(path_RUGFile, function(file) {
            read_fortran_data(
                file_path = file,
                col_widths_RUGFile = c(1, 3, 6, 10, 10, 10, 10),
                skip = 1
            )
        })
        df_MAP <- list(
            df_MAP_Kmin,
            df_MAP_Kflood
        )
        for (i in seq_along(MAP_RUGFile)) {
            if (nrow(df_MAP[[i]]) != nrow(MAP_RUGFile[[i]])) stop("Mage project has not the same size as the value estimated based on the MAP estimator")
            MAP_RUGFile[[i]]$V6 <- df_MAP[[1]]$Value

            write_RUGFile(
                RUG_path = path_RUGFile[i],
                RUGFile_data = MAP_RUGFile[[i]],
                RUG_format = "%1s%3d      %10.3f%10.3f%10.2f%10.2f"
            )
        }

        REPFile <- list.files(temp_path)[grep(list.files(temp_path), pattern = ".REP")]
        REPFile <- str_remove(REPFile, pattern = ".REP")

        for (i in 1:length(temp_path)) {
            setwd(temp_path[i])
            system2(MAGE_executable, args = REPFile[i], wait = TRUE)
        }

        runModel(
            workspace = temporal_dir,
            mod = mod_polynomials[[counter]],
            X = X,
            stout = NULL
        )
        counter <- counter + 1

        setwd(temporal_dir)

        sim <- read.table("Y.txt", header = TRUE)
        obs <- CalData[, 5:9]
        obs[obs == -9999] <- NA
        u_obs <- CalData[, 10:14]
        u_obs[u_obs == -9999] <- NA
        res <- obs - sim

        sim_event_mm_m3_s <- data.frame(
            X1_obs = X[, 1],
            X2_obs = X[, 2],
            X3_obs = X[, 3],
            X4_obs = X[, 4],
            Y1_obs = obs[, 1],
            Y1_sim = sim[, 1],
            Y2_sim = sim[, 2],
            Y1_res = res[, 1],
            Yu_z = u_obs[, 1]
        )
        setwd(dir_workspace)
    }

    Q_sim_vs_obs <- ZQdX_sim(
        sim_event_mm_m3_s = sim_event_mm_m3_s,
        Q_input = Q_observed,
        Qplot = TRUE,
        title_label = "MAP simulations vs boundary condition \n(Discharge upstream)",
        ylabel = "Q (L/s)"
    )
    ggsave(
        file.path(
            path_post_traitement,
            "Q_residuals.png"
        ),
        Q_sim_vs_obs,
        width = 20,
        height = 20,
        units = "cm"
    )
    save(Q_sim_vs_obs,
        file = file.path(path_post_traitement_data, "Plot_Q_sim_vs_obs.RData")
    )
    Z_sim_vs_obs <- ZQdX_sim(
        sim_event_mm_m3_s = sim_event_mm_m3_s,
        Q_input = NULL,
        Qplot = FALSE,
        title_label = "MAP simulations vs Observations (WSE)",
        ylabel = "Water surface elevation (mm)"
    )
    ggsave(
        file.path(
            path_post_traitement,
            "Plot_Z_residuals.png"
        ),
        Z_sim_vs_obs,
        width = 20,
        height = 20,
        units = "cm"
    )
    save(Z_sim_vs_obs,
        file = file.path(path_post_traitement_data, "Plot_Z_sim_vs_obs.RData")
    )
    save(sim_event_mm_m3_s,
        file = file.path(path_post_traitement_data, "Data_Z_sim_vs_obs.RData")
    )
}
