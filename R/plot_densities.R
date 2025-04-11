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
plot_densities <- function(ff_gated, df, device_colors, transformlist = NULL, density_n = 500, dfcol_grouping_samples = "Device") {
    relevant_columns <- flowCore::colnames(ff_gated[[1]])


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

    if(!all(names(is.null(ff_gated)))){
        stop("ff_gated must be a named list of flowFrame objects")
    }
    gated_dt <- lapply(ff_gated, function(x) {
        flowCore::exprs(x) |> data.table::as.data.table()
    }) |>
        data.table::rbindlist(idcol = "File", fill = TRUE) |>
        data.table::melt(id.vars = "File")

    compute_density <- function(x, transform_x) {
        x_noNA <- x[!is.na(x)]
        if (length(x_noNA) == 0) {
            res <- data.table::data.table(x = NA_real_, y = NA_real_)
        } else {
            d <- density(transform_x(x_noNA), n = density_n) # 512 points for smoother curves
            res <- data.table::data.table(x = d$x, y = d$y)
        }
        return(res)
    }
    densities <- gated_dt[,
        {
            compute_density(value, transformlist[[variable[[1]]]])
        },
        by = .(File, variable)
    ]
    df_part <- data.table::data.table(df)
    df_part <- df_part[, c("File", dfcol_grouping_samples[[1]], "Sample"), with = FALSE]
    densities <- densities[df_part, on = "File"]
    p0 <- ggplot2::ggplot(
        densities,
        ggplot2::aes(x = x, y = y, color = !!rlang::sym(dfcol_grouping_samples[[1]]))
    ) +
        ggplot2::geom_line() +
        # ggplot2::geom_area(ggplot2::aes(fill = Device), alpha = 0.2) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = 0, ymax = y, fill = !!rlang::sym(dfcol_grouping_samples[[1]])), alpha = 0.2) +
        # ggplot2::facet_grid(Sample~variable, scales = "free") +
        ggh4x::facet_grid2(Sample ~ variable, scales = "free_y", independent = "y") +
        # ggh4x::facet_grid2(Sample~variable, scales = "free",independent = "x") +
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
