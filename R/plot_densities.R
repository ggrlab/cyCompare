#' Plot Density Distributions for Gated Flow Cytometry Data
#'
#' This function generates density plots for multiple flow cytometry markers across different devices and samples.
#'
#' @param ff_gated A list of `flowFrame` objects containing gated flow cytometry data.
#' @param df A `data.table` containing metadata with at least the columns `"File"`, `"Device"`, and `"Sample"`.
#' @param device_colors An optional named vector of colors for each device to use in the plot.
#' @param transformlist A list of transformation functions for each marker. If a single function is provided, it will be applied to all markers.
#' @param density_n An integer specifying the number of points used to estimate density (default: 500).
#' @param relevant_columns A character vector of marker names to include in the density plots. If `NULL`, all markers will be included.
#' @param limit_density_quantile
#' A numeric value between 0 and 1 to limit the density to a certain quantile within each variable.
#' If `NA`, no limit is applied and all densities are shown to their full extent. Use
#' if one device shows an extreme density peak that is not representative of the other devices
#' and thus deforms the density plot. 0.95 is a good value to start with.
#'
#' @return A `ggplot2` object showing the density distributions of markers across samples and devices.
#'
#' @export
plot_densities <- function(
    ff_gated,
    df,
    device_colors,
    transformlist = NULL,
    density_n = 500, dfcol_grouping_samples = "Device",
    relevant_columns = NULL,
    limit_density_quantile = NA) {
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

    if (!all(names(is.null(ff_gated)))) {
        stop("ff_gated must be a named list of flowFrame objects")
    }
    gated_dt <- lapply(ff_gated, function(x) {
        flowCore::exprs(x) |> data.table::as.data.table()
    }) |>
        data.table::rbindlist(idcol = "File", fill = TRUE) |>
        data.table::melt(id.vars = "File")

    gated_dt <- gated_dt[variable %in% relevant_columns]
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
            compute_density(value, transformlist[[as.character(variable[[1]])]])
        },
        by = .(File, variable)
    ]
    if (!is.na(limit_density_quantile)) {
        # Limit the density to a certain quantile within each variable
        # calculate the top quantile for each variable
        top_quantile <- densities[, quantile(y, limit_density_quantile), by = variable]
        densities <- densities[top_quantile, on = "variable"][, y := pmin(y, V1)][, V1 := NULL]
    }
    df_part <- data.table::data.table(df)
    df_part <- df_part[, c("File", dfcol_grouping_samples[[1]], "Sample"), with = FALSE]
    densities <- densities[df_part, on = "File"]
    p0 <- ggplot2::ggplot(
        densities,
        ggplot2::aes(x = x, y = y, color = !!rlang::sym(dfcol_grouping_samples[[1]]))
    ) +
        ggplot2::geom_line() +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = 0, ymax = y, fill = !!rlang::sym(dfcol_grouping_samples[[1]])), alpha = 0.2) +
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
