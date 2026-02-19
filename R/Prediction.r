grid_user <- function(
    info_events_reaches) {
    n_prediction <- length(prediction_file)

    grid_user <- c()
    # Loop through each event
    for (event in info_events_reaches) {
        # Initialize an empty list to store all data frames
        all_dfs <- list()

        # Loop through each element in the event (skip the 'type' field)
        for (i in 2:length(event)) {
            all_dfs <- c(all_dfs, list(event[[i]]))
        }
        # Checks
        check_type_pred(event$type)

        # Combine all SU into one event
        df_event <- do.call(rbind, all_dfs)
        rownames(df_event) <- NULL

        if (event$type == "ZdX") {
            check_dX(df_event)

            # Generate discretization for each row and expand into a long-format data frame
            grid_temp <- do.call(rbind, lapply(1:nrow(df_event), function(i) {
                data.frame(
                    event = df_event$event[i],
                    reach = df_event$reach[i],
                    x = seq(from = df_event$xmin[i], to = df_event$xmax[i], length.out = df_event$nb_discretization[i]),
                    t = df_event$tmin[i]
                )
            }))
        } else if (event$type == "QdXT") {
            check_dXT(df_event)

            # Generate discretization for each row and expand into a long-format data frame
            grid_temp <- do.call(rbind, lapply(1:nrow(df_event), function(i) {
                data.frame(
                    event = df_event$event[i],
                    reach = df_event$reach[i],
                    x = df_event$xmin[i],
                    t = df_event$tmin[i]
                )
            }))
        } else if (event$type == "ZdT") {
            check_dT(df_event)

            # Generate discretization for each row and expand into a long-format data frame
            grid_temp <- do.call(rbind, lapply(1:nrow(df_event), function(i) {
                data.frame(
                    event = df_event$event[i],
                    reach = df_event$reach[i],
                    x = df_event$xmin[i],
                    t = seq(from = df_event$tmin[i], to = df_event$tmin[i], length.out = df_event$nb_discretization[i]),
                )
            }))
        } else {
            stop("Never arrive here, if not, a check must be created before")
        }
        grid_user <- rbind(grid_user, grid_temp)
    }
    return(grid_user = grid_user)
}

function(
    paths,
    X_pred,
    nX,
    nY) {
    # checks
    check_calibration_case(paths$path_polynomial)


    # Load data and model used during calibration
    load(file.path(paths$path_post_data, "BaM_objects.RData"))

    Caldata <- data$data
    X <- full_join(X_pred, Caldata[, 1:4], by = c("event", "reach", "x", "t")) %>%
        arrange(event, reach, x)

    names_file_prediction <- colnames(Caldata[, (nX + 1):(nX + nY)])
}
