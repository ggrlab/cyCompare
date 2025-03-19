plot_densities <- function(ff_gated, df, device_colors, transformlist = NULL, density_n = 500) {
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

    gated_dt <- lapply(ff_gated, function(x) {
        flowCore::exprs(x) |> data.table::as.data.table()
    }) |>
        data.table::rbindlist(idcol = "File") |>
        data.table::melt(id.vars = "File")

    compute_density <- function(x, transform_x) {
        d <- density(transform_x(x), n = density_n) # 512 points for smoother curves
        data.table::data.table(x = d$x, y = d$y)
    }
    densities <- gated_dt[, compute_density(value, transformlist[[variable[[1]]]]), by = .(File, variable)]
    densities <- densities[df[, c("File", "Device", "Sample")], on = "File"]
    p0 <- ggplot2::ggplot(
        densities,
        ggplot2::aes(x = x, y = y, color = Device)
    ) +
        ggplot2::geom_line() +
        # ggplot2::geom_area(ggplot2::aes(fill = Device), alpha = 0.2) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = 0, ymax = y, fill = Device), alpha = 0.2) +
        # ggplot2::facet_grid(Sample~variable, scales = "free") +
        ggh4x::facet_grid2(Sample ~ variable, scales = "free_y", independent = "y") +
        # ggh4x::facet_grid2(Sample~variable, scales = "free",independent = "x") +
        ggpubr::theme_pubclean() +
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
