#### pre-processing to obtain the Calibration data

assign_calibration_and_validation_data <- function(
    Input_Typology,
    Input_Model_Reach,
    CalData,
    All_observations) {
    # 1. Calibration set
    CalData_event <- c()
    for (i in 1:length(Input_Typology)) {
        # Remove the WSE at the nodes
        nodes_extraction <- Input_Model_Reach %>%
            filter(reach %in%
                Input_Typology[[i]]) %>%
            select(-reach)

        nodes_to_remove_cal <- unique(c(nodes_extraction$KP_start, nodes_extraction$KP_end))

        if (nrow(CalData %>%
            filter(id_reach_CAL %in% Input_Typology[[i]])) > 0) {
            if (i == 1) {
                CalData_event <- CalData %>%
                    mutate(Reach_groupped_Cal = case_when(
                        id_reach_CAL %in% Input_Typology[[i]] ~ names(Input_Typology)[i],
                        TRUE ~ Reach_groupped_Cal
                    )) %>%
                    group_by(id_case, Reach_groupped_Cal, id_reach_CAL) %>%
                    # If KP duplicated, take it off
                    filter(!KP %in% nodes_to_remove_cal) %>%
                    ungroup()
            } else {
                CalData_event_temp <- CalData %>%
                    mutate(Reach_groupped_Cal = case_when(
                        id_reach_CAL %in% Input_Typology[[i]] ~ names(Input_Typology)[i],
                        TRUE ~ Reach_groupped_Cal
                    )) %>%
                    group_by(id_case, Reach_groupped_Cal, id_reach_CAL) %>%
                    # If KP duplicated, take it off
                    filter(!KP %in% nodes_to_remove_cal) %>%
                    ungroup()
                CalData_event <- rbind(CalData_event, CalData_event_temp)
            }
        }
    }

    # 2. Validation set
    ValData_event <- All_observations %>%
        anti_join(CalData_event, by = colnames(All_observations)) %>%
        mutate(set = "validation")

    # 3. Combine into a single data frame
    CalValData <- bind_rows(CalData_event, ValData_event) %>% arrange(KP)

    return(CalValData)
}


constructor_CalData <- function(
    observed_data,
    do_manual_uncertainty = FALSE,
    sd_WSE_fixed = 0.05, # in meters
    sd_Q_fixed = 8, # in %
    sd_V_fixed = 3, # in %
    sd_Kmin_fixed = 5, # in m1/3/s
    sd_Kflood_fixed = 5 # in m1/3/s
    ) {
    size_df <- nrow(observed_data)
    # ---------------- WSE ----------------
    if ("WSE" %in% colnames(observed_data)) {
        Y <- data.frame(WSE = observed_data$WSE)

        if (do_manual_uncertainty) {
            Yu <- data.frame(Yu_WSE = rep(sd_WSE_fixed, size_df))
        } else {
            Yu <- data.frame(Yu_WSE = observed_data$Yu_WSE)
        }
    } else {
        Y <- data.frame(WSE = rep(-9999, size_df))
        Yu <- data.frame(Yu_WSE = rep(-9999, size_df))
    }

    # ---------------- Q ----------------
    if ("Q" %in% colnames(observed_data)) {
        Y$Q <- observed_data$Q

        if (do_manual_uncertainty) {
            Yu$Yu_Q <- rep(sd_Q_fixed, size_df)
        } else {
            Yu$Yu_Q <- observed_data$Yu_Q
        }
    } else {
        Y$Q <- rep(-9999, size_df)
        Yu$Yu_Q <- rep(-9999, size_df)
    }

    # ---------------- V ----------------
    if ("V" %in% colnames(observed_data)) {
        Y$V <- observed_data$V

        if (do_manual_uncertainty) {
            Yu$Yu_V <- rep(sd_V_fixed, size_df)
        } else {
            Yu$Yu_V <- observed_data$Yu_V
        }
    } else {
        Y$V <- rep(-9999, size_df)
        Yu$Yu_V <- rep(-9999, size_df)
    }

    # ---------------- Kmin ----------------
    if ("Kmin" %in% colnames(observed_data)) {
        Y$Kmin <- observed_data$Kmin

        if (do_manual_uncertainty) {
            Yu$Yu_Kmin <- rep(sd_Kmin_fixed, size_df)
        } else {
            Yu$Yu_Kmin <- observed_data$Yu_Kmin
        }
    } else {
        Y$Kmin <- rep(-9999, size_df)
        Yu$Yu_Kmin <- rep(-9999, size_df)
    }

    # ---------------- Kflood ----------------
    if ("Kflood" %in% colnames(observed_data)) {
        Y$Kflood <- observed_data$Kflood

        if (do_manual_uncertainty) {
            Yu$Yu_Kflood <- rep(sd_Kflood_fixed, size_df)
        } else {
            Yu$Yu_Kflood <- observed_data$Yu_Kflood
        }
    } else {
        Y$Kflood <- rep(-9999, size_df)
        Yu$Yu_Kflood <- rep(-9999, size_df)
    }

    return(list(
        Y = Y,
        Yu = Yu
    ))
}
