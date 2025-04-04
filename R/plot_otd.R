plot_otd <- function(ff_gated,
                     df,
                     device_colors = NULL,
                     transformlist,
                     n_mastersample = 1e4,
                     kwargs_loss = list(
                         loss = lossfun_hist,
                         verbose = FALSE,
                         write_intermediate = FALSE,
                         # Do not skip any comparisons
                         should_skip = function(i, j) FALSE,
                         take_time = FALSE,
                         return_as_matrix = TRUE
                     )) {
    # Convert the list of flowFrames into a flowSet
    gated_fs <- flowCore::flowSet(ff_gated)

    # If a single transformation function is provided, apply it to all markers
    if (length(transformlist) == 1) {
        fc_transformlist <- flowCore::transformList(
            flowCore::colnames(gated_fs[[1]]),
            transformlist
        )
    } else {
        fc_transformlist <- flowCore::transformList(
            flowCore::colnames(gated_fs[[1]]),
            transformlist
        )
    }

    # Apply transformation to the flowSet
    gated_fs_transformed <- flowCore::transform(gated_fs, fc_transformlist)

    dt_events <- flowCore::fsApply(
        gated_fs_transformed,
        simplify = FALSE,
        function(x) {
            flowCore::exprs(x) |> data.table::as.data.table()
        }
    )
    dt_events_bound <- data.table::rbindlist(dt_events, idcol = "File")
    mastersample <- dt_events_bound[sample(.N, n_mastersample), ][, File := NULL]
    loss_calculated <- do.call(loss_pairwise, c(
        list(
            datalist_A = list("mastersample" = mastersample),
            datalist_B = dt_events
        ), kwargs_loss
    ))
    loss_calc_mat <- loss_calculated[["dist"]]

    dists_df <- data.frame(t(loss_calculated[["dist"]])) |>
        tibble::rownames_to_column("File") |>
        dplyr::rename("OTD\nto mastersample" = mastersample)
    df_with_dists <- dplyr::left_join(df, dists_df, by = "File")
    p0 <- ggplot2::ggplot(
        df_with_dists, ggplot2::aes(x = Time, y = `OTD\nto mastersample`, color = Device, fill = Device)
    )
    p0 <- p0 + ggplot2::geom_point() + # Scatter plot of data points
        ggpubr::theme_pubr() + # Clean publication-ready theme
        ggplot2::theme(
            # remove y axis lines
            axis.line.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank()
        ) +
        ggplot2::geom_smooth(formula = y ~ x, method = "loess", se = TRUE, alpha = .2) + # Smoothed trend with confidence interval
        ggplot2::xlab("Time") # X-axis label

    # Apply custom device colors if provided
    if (!all(is.null(device_colors))) {
        p0 <- p0 +
            ggplot2::scale_color_manual(values = device_colors) + # Custom colors for lines
            ggplot2::scale_fill_manual(values = device_colors) # Custom colors for fill
    }

    return(p0) # Return the ggplot object
}
