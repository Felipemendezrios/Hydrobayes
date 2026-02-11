plot_CalData <- function(
    CalData,
    scales_free = "free_y",
    wrap = TRUE) {
    ##########################################
    # Plot WSE (Output number 1)
    ##########################################
    if (any(CalData$WSE != -9999)) {
        CalData$WSE <- convert_9999_to_NA(CalData$WSE)
        plot_WSE <-
            ggplot(data = CalData, aes(
                x = x,
                y = WSE,
                col = factor(reach)
            ))

        if (any(CalData$Yu_WSE != 0)) {
            plot_WSE <- plot_WSE +
                geom_errorbar(aes(
                    ymin = WSE - qnorm(0.975) * Yu_WSE,
                    ymax = WSE + qnorm(0.975) * Yu_WSE
                ))
        }

        plot_WSE <- plot_WSE +
            geom_point() +
            labs(
                x = "Streamwise position (m)",
                y = "Water surface elevation (m)",
                col = "Reaches in the \nhydraulic model"
            ) +
            theme_bw() +
            theme(
                plot.title = element_text(hjust = 0.5),
                legend.title = element_text(hjust = 0.5)
            )

        if (wrap) {
            plot_WSE <- plot_WSE +
                facet_wrap(~event, scales = scales_free, ncol = 1)
        }

        plot_WSE <- plot_WSE +
            theme(
                strip.text = element_text(size = 12, face = "bold"),
                axis.title.x = element_text(size = 13), # x-axis title
                axis.title.y = element_text(size = 13), # y-axis title
                axis.text.x = element_text(size = 12, hjust = 1), # x-axis ticks
                axis.text.y = element_text(size = 12), # y-axis ticks
                legend.title = element_text(size = 14, face = "bold"),
                legend.text = element_text(size = 12)
            )
    } else {
        plot_WSE <- NULL
    }
    ##########################################
    # Plot Q (Output number 2)
    ##########################################
    if (any(CalData$Q != -9999)) {
        CalData$Q <- convert_9999_to_NA(CalData$Q)
        plot_Q <-
            ggplot(data = CalData, aes(
                x = x,
                y = Q,
                col = factor(reach)
            ))

        if (any(CalData$Yu_Q != 0)) {
            plot_Q <- plot_Q +
                geom_errorbar(aes(
                    ymin = Q - qnorm(0.975) * Yu_Q * Q / 100,
                    ymax = Q + qnorm(0.975) * Yu_Q * Q / 100
                ))
        }

        plot_Q <- plot_Q +
            geom_point() +
            labs(
                x = "Streamwise position (m)",
                y = expression("Discharge (m"^
                    {
                        3
                    } * "/s)"),
                col = "Reaches in the \nhydraulic model"
            ) +
            theme_bw() +
            theme(
                plot.title = element_text(hjust = 0.5),
                legend.title = element_text(hjust = 0.5)
            )
        if (wrap) {
            plot_Q <- plot_Q +
                facet_wrap(~event, scales = scales_free, ncol = 1)
        }

        plot_Q <- plot_Q +
            theme(
                strip.text = element_text(size = 12, face = "bold"),
                axis.title.x = element_text(size = 13), # x-axis title
                axis.title.y = element_text(size = 13), # y-axis title
                axis.text.x = element_text(size = 12, hjust = 1), # x-axis ticks
                axis.text.y = element_text(size = 12), # y-axis ticks
                legend.title = element_text(size = 14, face = "bold"),
                legend.text = element_text(size = 12)
            )
    } else {
        plot_Q <- NULL
    }



    ##########################################
    # Plot V (Output number 3)
    ##########################################
    if (any(CalData$V != -9999)) {
        CalData$V <- convert_9999_to_NA(CalData$V)
        plot_V <-
            ggplot(data = CalData, aes(
                x = x,
                y = V,
                col = factor(reach)
            ))

        if (any(CalData$Yu_V != 0)) {
            plot_V <- plot_V +
                geom_errorbar(aes(
                    ymin = Q - qnorm(0.975) * Yu_V * Q / 100,
                    ymax = Q + qnorm(0.975) * Yu_V * Q / 100
                ))
        }

        plot_V <- plot_V +
            geom_point() +
            labs(
                x = "Streamwise position (m)",
                y = "Velocity (m/s)",
                col = "Reaches in the \nhydraulic model"
            ) +
            theme_bw() +
            theme(
                plot.title = element_text(hjust = 0.5),
                legend.title = element_text(hjust = 0.5)
            )
        if (wrap) {
            plot_V <- plot_V +
                facet_wrap(~event, scales = scales_free, ncol = 1)
        }

        plot_V <- plot_V +
            theme(
                strip.text = element_text(size = 12, face = "bold"),
                axis.title.x = element_text(size = 13), # x-axis title
                axis.title.y = element_text(size = 13), # y-axis title
                axis.text.x = element_text(size = 12, hjust = 1), # x-axis ticks
                axis.text.y = element_text(size = 12), # y-axis ticks
                legend.title = element_text(size = 14, face = "bold"),
                legend.text = element_text(size = 12)
            )
    } else {
        plot_V <- NULL
    }



    ##########################################
    # Plot Kmin (Output number 4)
    ##########################################
    if (any(CalData$Kmin != -9999)) {
        CalData$Kmin <- convert_9999_to_NA(CalData$Kmin)
        plot_Kmin <-
            ggplot(data = CalData, aes(
                x = x,
                y = Kmin,
                col = factor(reach)
            ))

        if (any(CalData$Yu_Kmin != 0)) {
            plot_Kmin <- plot_Kmin +
                geom_errorbar(aes(
                    ymin = Kmin - qnorm(0.975) * Yu_Kmin,
                    ymax = Kmin + qnorm(0.975) * Yu_Kmin
                ))
        }

        plot_Kmin <- plot_Kmin +
            geom_point() +
            labs(
                x = "Streamwise position (m)",
                y = expression("Friction coefficient (m"^
                    {
                        1 / 3
                    } * "/s)"),
                col = "Reaches in the \nhydraulic model"
            ) +
            theme_bw() +
            theme(
                plot.title = element_text(hjust = 0.5),
                legend.title = element_text(hjust = 0.5)
            )

        if (wrap) {
            plot_Kmin <- plot_Kmin +
                facet_wrap(~event, scales = scales_free, ncol = 1)
        }

        plot_Kmin <- plot_Kmin +
            theme(
                strip.text = element_text(size = 12, face = "bold"),
                axis.title.x = element_text(size = 13), # x-axis title
                axis.title.y = element_text(size = 13), # y-axis title
                axis.text.x = element_text(size = 12, hjust = 1), # x-axis ticks
                axis.text.y = element_text(size = 12), # y-axis ticks
                legend.title = element_text(size = 14, face = "bold"),
                legend.text = element_text(size = 12)
            )
    } else {
        plot_Kmin <- NULL
    }


    ##########################################
    # Plot Kflood(Output number 5)
    ##########################################
    if (any(CalData$Kflood != -9999)) {
        CalData$Kflood <- convert_9999_to_NA(CalData$Kflood)
        plot_Kflood <-
            ggplot(data = CalData, aes(
                x = x,
                y = Kflood,
                col = factor(reach)
            ))

        if (any(CalData$Yu_Kflood != 0)) {
            plot_Kflood <- plot_Kflood +
                geom_errorbar(aes(
                    ymin = Kflood - qnorm(0.975) * Yu_Kflood,
                    ymax = Kflood + qnorm(0.975) * Yu_Kflood
                ))
        }

        plot_Kflood <- plot_Kflood +
            geom_point() +
            labs(
                x = "Streamwise position (m)",
                y = expression("Friction coefficient (m"^
                    {
                        1 / 3
                    } * "/s)"),
                col = "Reaches in the \nhydraulic model"
            ) +
            theme_bw() +
            theme(
                plot.title = element_text(hjust = 0.5),
                legend.title = element_text(hjust = 0.5)
            )
        if (wrap) {
            plot_Kflood <- plot_Kflood +
                facet_wrap(~event, scales = scales_free, ncol = 1)
        }

        plot_Kflood <- plot_Kflood +
            theme(
                strip.text = element_text(size = 12, face = "bold"),
                axis.title.x = element_text(size = 13), # x-axis title
                axis.title.y = element_text(size = 13), # y-axis title
                axis.text.x = element_text(size = 12, hjust = 1), # x-axis ticks
                axis.text.y = element_text(size = 12), # y-axis ticks
                legend.title = element_text(size = 14, face = "bold"),
                legend.text = element_text(size = 12)
            )
    } else {
        plot_Kflood <- NULL
    }


    return(list(
        plot_WSE = plot_WSE,
        plot_Q = plot_Q,
        plot_V = plot_V,
        plot_Kmin = plot_Kmin,
        plot_Kflood = plot_Kflood
    ))
}
