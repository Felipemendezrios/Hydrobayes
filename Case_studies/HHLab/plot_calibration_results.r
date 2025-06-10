cat("\014")
rm(list = ls())

# Plots:
library(ggplot2)
library(RBaM)
library(dplyr)
library(tidyr)
library(stringr)
library(stats)

workspace <- "/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/Case_studies/HHLab"
setwd(workspace)

# Here, define calibration type
# 'Calibration_time_series'
# 'Calibration_water_surface_profiles'
path_results <- "Calibration_water_surface_profiles"
path_model_mage_global <- "model_mage"
final_results <- TRUE
n_degree_max <- 4

# Check
check_cal_WS_profiles <- path_results == "Calibration_water_surface_profiles"
# Create a sequence of numbers: time fixed to extract simulation
if (check_cal_WS_profiles) {
    sequence <- seq(0.95, 3.95, by = 1)

    sequence_all <- c(
        sequence
    )

    # Find the maximum number of decimal places
    max_decimals <- max(sapply(sequence_all, function(x) nchar(sub(".*\\.", "", as.character(x)))))

    sequence_all <- sort(round(sequence_all, max_decimals))

    # Mage extraire arguments
    mage_extraire_args <- lapply(sequence_all, function(x) paste0("ZdX 1 h", trimws(format(x, nsmall = max_decimals))))
} else {
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
    mage_extraire_args <- lapply(sequence_all, function(x) paste("ZdT 1", format(x, nsmall = max_decimals)))
}

# Adjust size of plot : in pixel!
dimensions_plot_friction <- c(
    1300, # Width
    900 # Height
)
dimensions_plot_sim_vs_obs <- if (check_cal_WS_profiles) {
    c(
        1600, # Width
        1800 # Height
    )
} else {
    c(
        2000, # Width
        2000 # Height
    )
}
legend_name_col <- ifelse(check_cal_WS_profiles,
    "Calibration : water surface profiles",
    "Calibration : time series"
)
legend_wrap_name <- ifelse(check_cal_WS_profiles,
    "Time: ",
    "Position: "
)
legend_x_name <- ifelse(check_cal_WS_profiles,
    "Lengthwise position(meters)",
    "Time (hours)"
)

# Do not touch
n_degree_seq <- seq(0, n_degree_max, 1)
nY <- length(mage_extraire_args)
col_widths <- c(1, 3, 6, 10, 10, 10, 10)
MAGE_executable <- "/home/famendezrios/Documents/Softwares/pamhyr2/mage8/mage"
Mage_extraire_executable <- "/home/famendezrios/Documents/Softwares/pamhyr2/mage8/mage_extraire"

# Get column names from mage_extraire_args
column_names_mage_extraire_args <-
    sapply(mage_extraire_args, function(arg) {
        # Get the first letter (e.g., "Z" or "Q")
        prefix <- substr(arg, 1, 3)

        # Split by spaces and take the last element (the number)
        split_text <- strsplit(arg, " ")[[1]]
        # value <- split_text[length(split_text)]
        value <- stringr::str_extract(arg, "\\d+\\.\\d+$")
        # Combine prefix and value to create column name
        paste0(prefix, "_", value)
    })

for (n_degree in n_degree_seq) {
    path_polynomial <- file.path(path_results, paste0("n_", n_degree))
    path_post_traitement <- file.path(path_polynomial, "post_traitement")
    path_post_traitement_data <- file.path(path_post_traitement, "RData")
    path_model_mage <- file.path(path_polynomial, path_model_mage_global)


    if (!dir.exists(path_post_traitement)) {
        dir.create(path_post_traitement)
    }

    if (!dir.exists(path_post_traitement_data)) {
        dir.create(path_post_traitement_data)
    }
    MCMC <- readMCMC(file.path(path_polynomial, "Results_MCMC.txt"))

    png(
        filename = file.path(path_post_traitement, "MCMC_tracePlot.png"),
        res = 180, width = 1000, height = 1000
    )
    plots <- tracePlot(MCMC)
    gridExtra::grid.arrange(grobs = plots, ncol = 3)
    dev.off()

    png(
        filename = file.path(path_post_traitement, "MCMC_densityPlot.png"),
        res = 180, width = 1000, height = 1000
    )
    # Density plot for each parameter
    plots <- densityPlot(MCMC)
    gridExtra::grid.arrange(grobs = plots, ncol = 3)
    dev.off()

    # Friction coefficient estimation:
    # Get position from the real geometry
    position <- read.table(file.path(path_polynomial, "vector_legendre.txt"), header = TRUE)

    matrix_zFileKmin <- read.table(file.path(path_polynomial, "Zfile_Kmin.txt"), header = TRUE)


    # Get all MCMC sampling for creating the envelope
    if (final_results) {
        if (!file.exists(file.path(path_polynomial, "Results_Cooking.txt"))) stop("MCMC is still running or calculation is not going to the end. Verify if calibration is already finished or verify that calibration has not error messages. Please put final_results = FALSE as input")

        results_MCMC_cooked_read <- read.table(file.path(path_polynomial, "Results_Cooking.txt"), header = TRUE)
        results_MCMC <- results_MCMC_cooked_read[, 1:(n_degree + 1)]

        # Get MAP simulation
        summary_data <- read.table(file.path(path_polynomial, "Results_Summary.txt"), header = TRUE)
        MAP_param_matrix <- as.numeric(summary_data["MaxPost", 1:(n_degree + 1)])
    } else {
        results_MCMC_sampling <- read.table(file.path(path_polynomial, "Results_MCMC.txt"), header = TRUE)
        results_MCMC <- results_MCMC_sampling[, 1:(n_degree + 1)]

        # Get MAP simulation: from current analysis without burning and slim
        MAP_param_matrix <- as.numeric(results_MCMC_sampling[which.max(results_MCMC_sampling$LogPost), 1:(n_degree + 1)])
    }

    k_estimated_all <- as.data.frame(as.matrix(matrix_zFileKmin) %*% as.matrix(t(results_MCMC)))

    k_estimated_all$KP <- position[, 1]

    # Convert to long format
    df_MCMC_sampling <- pivot_longer(
        k_estimated_all,
        cols = -KP,
        values_to = "Value"
    ) %>%
        select(KP, Value) %>%
        mutate(ID = "MCMC Sampling")

    k_estimated_MAP <- as.matrix(matrix_zFileKmin) %*% MAP_param_matrix

    df_MAP <- data.frame(
        KP = position[, 1],
        Value = k_estimated_MAP,
        ID = "MAP"
    )

    # Get 95% uncertainty for envelope curve : create ribbon data from MCMC
    df_envelope <- df_MCMC_sampling %>%
        filter(ID == "MCMC Sampling") %>%
        group_by(KP) %>%
        summarise(
            ymin = quantile(Value, probs = 0.025, na.rm = TRUE),
            ymax = quantile(Value, probs = 0.975, na.rm = TRUE),
            ID = "Parametric\nuncertainty", # so we can map to fill
            .groups = "drop"
        )

    plot_spatial_friction <- ggplot() +
        geom_ribbon(
            data = df_envelope,
            aes(x = KP, ymin = ymin, ymax = ymax, fill = ID)
        ) +
        geom_line(
            data = df_MAP,
            aes(x = KP, y = Value, color = ID)
        ) +
        scale_fill_manual(values = c("Parametric\nuncertainty" = "pink")) +
        scale_color_manual(values = c("MAP" = "black")) +
        labs(
            title = "Friction coefficient estimation \nwith parametric uncertainty",
            x = "Lengthwise position (meters)",
            y = expression("Friction coefficient (m"^
                {
                    1 / 3
                } * "/s)"),
            fill = "95% credibility\ninterval",
            color = NULL
        ) +
        theme_bw() +
        theme(
            plot.title = element_text(hjust = 0.5),
            legend.title = element_text(hjust = 0.5)
        )
    # Save plot and data
    ggsave(
        filename = file.path(path_post_traitement, "friction_estimation.png"),
        plot = plot_spatial_friction,
        dpi = 300, width = dimensions_plot_friction[1], height = dimensions_plot_friction[2],
        units = "px"
    )
    save(plot_spatial_friction,
        file = file.path(path_post_traitement_data, "friction_estimation_plot.RData")
    )

    ls_spatial_friction <- list(
        df_envelope = df_envelope,
        df_MAP = df_MAP
    )
    save(ls_spatial_friction,
        file = file.path(path_post_traitement_data, "friction_estimation_data.RData")
    )

    # Residuals simulation vs observations
    # Save Calibration Data file
    CalData_raw <- read.table(file.path(path_polynomial, "CalibrationData.txt"), header = TRUE)

    if (final_results) {
        if (!file.exists(file.path(path_polynomial, "Results_Residuals.txt"))) stop("MCMC is still running or calculation is not going to the end. Verify if calibration is already finished or verify that calibration has not error messages. Please put final_results = FALSE as input")

        residuals <- read.table(file.path(path_polynomial, "Results_Residuals.txt"), header = TRUE)
        residuals_extraction <- residuals[, c(1, 3:(2 + nY), (nY * 2 + 2 + 1):(nY * 4 + 2))]
        colnames_residuals <- colnames(residuals_extraction)

        # This leds to be sure that we are comparing the same data
        check_df <- residuals[, c(1, 3:(2 + nY))]
        colnames(check_df) <- c("X1_obs", column_names_mage_extraire_args)
        # Check
        if (any(is.na(match(colnames(check_df)[-1], colnames(CalData_raw)[-1])))) stop("mage_extraire_args are not defined in the same order as during calibration step")

        for (i in 1:ncol(residuals_extraction)) {
            residuals_extraction[which(residuals_extraction[, i] == -9999), i] <- NA
        }
    } else {
        # Residuals file is not ready, so I need to create by myself with MAP estimation from sampled data while MCMC is turning

        CalData <- CalData_raw[, c(1:(1 + length(mage_extraire_args)))]

        colnames(CalData)[1] <- "X"
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
        temp_path <- file.path(temporal_dir, path_model_mage_global)

        read_fortran_data <- function(file_path, col_widths, skip = 0) {
            # Read the file with the fixed-width format
            data <- read.fwf(file_path, widths = col_widths, header = FALSE, skip = skip)
            return(data)
        }
        path_RUGFile <- file.path(
            temp_path,
            list.files(temp_path)[grep(list.files(temp_path), pattern = ".RUG")]
        )
        MAP_RUGFile <- read_fortran_data(
            file_path = path_RUGFile,
            col_widths = col_widths,
            skip = 1
        )
        MAP_RUGFile$V6 <- df_MAP$Value

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
            RUG_path = path_RUGFile,
            RUG_id_reach = MAP_RUGFile$V2,
            RUG_KP_start = MAP_RUGFile$V4,
            RUG_KP_end = MAP_RUGFile$V5,
            RUG_Kmin = MAP_RUGFile$V6,
            RUG_Kmoy = MAP_RUGFile$V7,
            RUG_format = "%1s%3d      %10.0f%10.0f%10.2f%10.2f"
        )

        REPFile <- list.files(temp_path)[grep(list.files(temp_path), pattern = ".REP")]
        REPFile <- str_remove(REPFile, pattern = ".REP")

        setwd(temp_path)
        # Run Mage
        system2(
            command = MAGE_executable,
            args = REPFile,
            wait = TRUE
        )
        # Sys.sleep(0.1)
        counter <- 1
        sim <- res <- obs <- data.frame(X1_obs = CalData[, 1])
        # get sim at the grid, estimate residuals
        for (i in mage_extraire_args) {
            new_name <- paste0(REPFile, "_", paste0(strsplit(i, " ")[[1]], collapse = "_"), ".res")
            # Step 1: Run executable (non-blocking)
            system2(
                command = Mage_extraire_executable,
                args = c(REPFile, i),
                wait = TRUE
            )
            res_file <- paste0(REPFile, ".res")
            wait_time <- 0
            while (!file.exists(res_file) && wait_time < 10) {
                Sys.sleep(0.2)
                wait_time <- wait_time + 0.2
            }
            # Rename the output file
            if (file.exists(res_file)) {
                file.remove(new_name)
                file.rename(
                    from = res_file,
                    to = new_name
                )
            } else {
                warning(sprintf("Output file %s.res not found after running with arg '%s'", REPFile, i))
            }
            # Get simulation :
            sim_temp <- read.table(new_name, comment.char = "*")
            # This leds to be sure that we are comparing the same data
            colnames(sim_temp) <- c("X1_obs", column_names_mage_extraire_args[counter])

            # Check
            if (is.na(match(colnames(sim_temp)[-1], colnames(CalData)))) stop("mage_extraire_args are not defined in the same order as during calibration step")
            obs_temp <- CalData[, match(colnames(sim_temp)[-1], colnames(CalData))]
            obs_temp[obs_temp == -9999] <- NA
            res_temp <- obs_temp - sim_temp[, 2]

            sim <- cbind(sim, sim_temp[, -1])
            colnames(sim)[counter + 1] <- paste0("Y", counter, "_sim")
            res <- cbind(res, res_temp)
            colnames(res)[counter + 1] <- paste0("Y", counter, "_res")
            obs <- cbind(obs, obs_temp)
            colnames(obs)[counter + 1] <- paste0("Y", counter, "_obs")
            counter <- counter + 1
        }
        residuals_extraction <- obs %>%
            left_join(sim, by = "X1_obs") %>%
            left_join(res, by = "X1_obs")
        setwd(workspace)
        # Delete
        unlink(temp_path, recursive = TRUE)
    }

    col_names <- paste0("Y", 1:nY)
    position_labels <- paste(legend_wrap_name, sequence_all)

    # Create named vector to use in recoding
    name_map <- setNames(position_labels, col_names)
    # Save names of each wrap plot ordered
    save(name_map, file = file.path(path_post_traitement_data, "names_plot_ordered.RData"))

    # Add uncertainty
    unc_CalData <- CalData_raw[, grep(colnames(CalData_raw), pattern = "Yu_")]
    colnames(unc_CalData) <- paste0("Yu_Y", 1:nY)
    residuals_extraction <- cbind(residuals_extraction, unc_CalData)


    # Step 1: Gather observed and simulated values into long format
    sim_obs_output_variable_long <- residuals_extraction %>%
        # Select columns searched
        select(X1_obs, matches("_obs$"), matches("_sim$")) %>%
        # Pivot observed and simulated values to long format
        pivot_longer(
            cols = matches("(Y\\d+)_(obs|sim)$"),
            names_to = c("variable", "type"),
            names_pattern = "(Y\\d+)_(obs|sim)",
            values_to = "value"
        ) %>%
        # Join with uncertainty values
        left_join(
            residuals_extraction %>%
                select(X1_obs, starts_with("Yu_")) %>%
                pivot_longer(
                    cols = starts_with("Yu_"),
                    names_to = "variable",
                    names_pattern = "Yu_(Y\\d+)",
                    values_to = "value_unc"
                ),
            by = c("X1_obs", "variable")
        ) %>%
        # Keep uncertainty only for observed rows and clean -9999
        mutate(
            value_unc = ifelse(type == "obs", value_unc, NA),
            value_unc = na_if(value_unc, -9999)
        ) %>%
        # Remove rows where both value and value_unc are NA
        filter(!(is.na(value) & is.na(value_unc))) %>%
        # Scale value, rename variables, and order levels
        mutate(
            value = value * 1000,
            value_unc = value_unc * 1000,
            variable = recode(variable, !!!name_map),
            variable = factor(variable, levels = position_labels) # enforce order
        )

    # Step 2: Plot
    conf_level <- 0.95 # Hard-coded
    z_val <- qnorm(1 - (1 - conf_level) / 2)
    sim_obs_output_variable_long <- sim_obs_output_variable_long %>%
        mutate(
            # Set drawing order: obs drawn after sim
            type = factor(type, levels = c("sim", "obs")),
            alpha_val = ifelse(type == "sim", 0.4, 1)
        )

    plot_sim_obs_cal_WS_profiles <-
        ggplot(
            data = sim_obs_output_variable_long,
            aes(
                x = X1_obs, y = value, color = type
            )
        ) +
        geom_line(
            data = filter(
                sim_obs_output_variable_long,
                type == "sim"
            ),
            aes(
                alpha = alpha_val,
                group = variable
            ),
            linewidth = 1
        ) +
        geom_point(
            data = filter(sim_obs_output_variable_long, type == "obs"),
            aes(
                x = X1_obs, y = value, color = type
            )
        ) +
        geom_errorbar(
            data = filter(sim_obs_output_variable_long, type == "obs"),
            aes(
                ymin = value - z_val * value_unc,
                ymax = value + z_val * value_unc,
                alpha = alpha_val,
            ),
            na.rm = TRUE
        ) +
        scale_alpha_identity() +
        facet_wrap(~variable, scales = "free_y") +
        theme_bw() +
        labs(
            x = legend_x_name,
            y = "Stage (mm)",
            title = paste0("Comparison of simulated and observed \nwater surface profiles \n", legend_name_col),
            col = "Data (MAP)"
        ) +
        theme(
            strip.text = element_text(size = 12),
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5)
        )
    # Save plot and data
    ggsave(
        filename = file.path(path_post_traitement, "output_var_comparison.png"),
        plot = plot_sim_obs_cal_WS_profiles,
        dpi = 300, width = dimensions_plot_sim_vs_obs[1], height = dimensions_plot_sim_vs_obs[2],
        units = "px"
    )
    save(plot_sim_obs_cal_WS_profiles,
        file = file.path(path_post_traitement_data, "output_var_comparison_plot.RData")
    )

    save(sim_obs_output_variable_long,
        file = file.path(path_post_traitement_data, "output_var_comparison_data.RData")
    )
    # Residuals
    resdiuals_sim_obs_output_variable_long <- residuals_extraction %>%
        # Select columns searched
        select(X1_obs, matches("_res$")) %>%
        pivot_longer(
            cols = matches("(Y\\d+)_(res)$"),
            names_to = c("variable", "type"),
            names_pattern = "(Y\\d+)_(res)",
            values_to = "value"
        ) %>%
        # Join with uncertainty values
        left_join(
            residuals_extraction %>%
                select(X1_obs, starts_with("Yu_")) %>%
                pivot_longer(
                    cols = starts_with("Yu_"),
                    names_to = "variable",
                    names_pattern = "Yu_(Y\\d+)",
                    values_to = "value_unc"
                ),
            by = c("X1_obs", "variable")
        ) %>%
        # Keep uncertainty only for observed rows and clean -9999
        mutate(
            value_unc = ifelse(type == "res", value_unc, NA),
            value_unc = na_if(value_unc, -9999)
        ) %>%
        # Remove rows where both value and value_unc are NA
        filter(!(is.na(value) & is.na(value_unc))) %>%
        # Scale value, rename variables, and order levels
        mutate(
            value = value * 1000,
            value_unc = value_unc * 1000,
            variable = recode(variable, !!!name_map),
            variable = factor(variable, levels = position_labels) # enforce order
        )

    # Step 2: Plot
    plot_res_cal_WS_profiles <-
        ggplot(
            resdiuals_sim_obs_output_variable_long,
            aes(
                x = X1_obs, y = value, color = type,
                ymin = value - z_val * value_unc,
                ymax = value + z_val * value_unc,
            )
        ) +
        geom_hline(yintercept = 0, linetype = "dashed", col = "black") +
        geom_point() +
        geom_errorbar(na.rm = TRUE) +
        facet_wrap(~variable, scales = "fixed") +
        theme_bw() +
        labs(
            x = legend_x_name,
            y = "Residuals (mm)",
            title = paste0("Residuals water surface profiles \n", legend_name_col),
            col = "Data (MAP)"
        ) +
        theme(
            strip.text = element_text(size = 12),
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5)
        )
    # Save plot and data
    ggsave(
        filename = file.path(path_post_traitement, "residual_output_var_comparison.png"),
        plot = plot_res_cal_WS_profiles,
        dpi = 300, width = dimensions_plot_sim_vs_obs[1], height = dimensions_plot_sim_vs_obs[2],
        units = "px"
    )
    save(plot_res_cal_WS_profiles,
        file = file.path(path_post_traitement_data, "residual_output_var_comparison_plot.RData")
    )

    save(resdiuals_sim_obs_output_variable_long,
        file = file.path(path_post_traitement_data, "residual_output_var_comparison_data.RData")
    )
}
