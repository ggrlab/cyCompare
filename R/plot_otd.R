#' Plot OTD-to-Reference Distance Over Time
#'
#' This function visualizes OTD values from each sample to the reference sample over time,
#' optionally colored by device.
#'
#' @param df A `data.table` with metadata. Must contain `"File"`, `"Time"`, and `"Device"` columns.
#' @param otd_df A `data.frame` as returned by `otd_to_reference()`.
#' @param device_colors Optional named vector of colors for each device.
#'
#' @return A `ggplot2` object.
#' @export
plot_otd <- function(df, otd_df, device_colors = NULL) {
    # Merge distance values with metadata
    df_with_dists <- dplyr::left_join(df, otd_df, by = "File")

    # Build plot
    p <- ggplot2::ggplot(
        df_with_dists,
        ggplot2::aes(x = Time, y = OTD_to_referencesample, color = Device, fill = Device)
    ) +
        ggplot2::geom_point() +
        ggplot2::geom_smooth(formula = y ~ x, method = "loess", se = TRUE, alpha = .2) +
        ggpubr::theme_pubr() +
        ggplot2::xlab("Time") +
        ggplot2::theme(
            axis.line.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank()
        )

    # Apply colors if given
    if (!all(is.null(device_colors))) {
        p <- p +
            ggplot2::scale_color_manual(values = device_colors) +
            ggplot2::scale_fill_manual(values = device_colors)
    }

    return(p)
}
