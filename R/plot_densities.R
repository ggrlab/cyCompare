#' Plot Density Distributions for Gated Flow Cytometry Data
#'
#' This function generates density plots for multiple flow cytometry markers across different devices and samples.
#'
#' @param ff_gated A list of `flowFrame` objects containing gated flow cytometry data.
#' @param df A `data.table` containing metadata with at least the columns `"File"`, `"Device"`, and `"Sample"`.
#' @param device_colors An optional named vector of colors for each device to use in the plot.
#' @param transformlist A list of transformation functions for each marker. If a single function is provided, it will be applied to all markers.
#' @param density_n An integer specifying the number of points used to estimate density (default: 500).
#'
#' @return A `ggplot2` object showing the density distributions of markers across samples and devices.
#'
#' @export
plot_densities <- function(ff_gated, df, device_colors, transformlist = NULL, density_n = 500, relevant_columns = NULL) {
    if (all(is.null(relevant_columns))) {
        relevant_columns <- flowCore::colnames(ff_gated[[1]])
    }

    # Transformation functions
    if (length(transformlist) == 1 && !is.null(transformlist)) {
        if (is.function(transformlist)) {
            transformlist <- list(transformlist)
        }
        transformlist <- setNames(
            rep(transformlist, length(relevant_columns)),
            relevant_columns
        )
    }
    for (col_x in flowCore::colnames(ff_gated[[1]])[!flowCore::colnames(ff_gated[[1]]) %in% relevant_columns]) {
        transformlist[[col_x]] <- identity
    }

    gated_dt <- lapply(ff_gated, function(x) {
        flowCore::exprs(x) |> data.table::as.data.table()
    }) |>
        data.table::rbindlist(idcol = "File") |>
        data.table::melt(id.vars = "File")

    gated_dt <- gated_dt[variable %in% relevant_columns]
    compute_density <- function(x, transform_x, y) {
        d <- density(transform_x(x), n = density_n) # 512 points for smoother curves
        data.table::data.table(x = d$x, y = d$y)
    }
    densities <- gated_dt[, compute_density(value, transformlist[[as.character(variable[[1]])]]), by = .(File, variable)]
    densities <- densities[df[, c("File", "Device", "Sample")], on = "File"]
    p0 <- ggplot2::ggplot(
        densities,
        ggplot2::aes(x = x, y = y, color = Device)
    ) +
        ggplot2::geom_line() +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = 0, ymax = y, fill = Device), alpha = 0.2) +
        ggh4x::facet_grid2(Sample ~ variable, scales = "free", independent = "y") +
        ggpubr::theme_pubr() +
        ggplot2::theme(
            # remove axis lines
            axis.line.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank()
        ) +
        ggplot2::ylab("Density") +
        ggplot2::xlab("Transformed MFI") +
        ggplot2::theme(
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank()
        )

    # Apply custom device colors if provided
    if (!all(is.null(device_colors))) {
        p0 <- p0 +
            ggplot2::scale_color_manual(values = device_colors) + # Custom colors for lines
            ggplot2::scale_fill_manual(values = device_colors) # Custom colors for fill
    }

    return(p0) # Return the ggplot object
}
