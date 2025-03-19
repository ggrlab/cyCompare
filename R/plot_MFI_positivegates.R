plot_MFI_positivegates <- function(dt_count_mfi, marker_to_gate, device_colors, transformlist = NULL) {
    relevant_gates <- unique(unlist(marker_to_gate))
    marker_to_gate_dt <- lapply(marker_to_gate, data.table::as.data.table) |> data.table::rbindlist(idcol = "marker")
    data.table::setnames(marker_to_gate_dt, "V1", "pop")

    dt_count_mfi_relevant <- dt_count_mfi[pop %in% relevant_gates]
    # Merge the 'dt_count_mfi_relevant' data table with 'marker_to_gate_dt' based on the 'pop' column.
    # Join where rows from 'marker_to_gate_dt' are matched with rows from 'dt_count_mfi_relevant'
    # using the 'pop' column as the key.
    # The result will include all columns from 'marker_to_gate_dt' and the matching columns
    # from 'dt_count_mfi_relevant'.
    dt <- marker_to_gate_dt[dt_count_mfi_relevant, on = "pop"]
    dt_long <- data.table::melt(dt, measure.vars = names(marker_to_gate), value.name = "MFI", variable.name = "variable")
    dt_medians <- dt_long[, .(mfi_all_gates = median(MFI)), by = .(File, Device, Time, marker, variable)]
    dt_medians_relevant <- dt_medians[marker == variable]
    for (unique_marker in unique(dt_medians_relevant$marker)) {
        # For each marker, apply the corresponding transformation function to the 'mfi_all_gates' column.
        tryCatch(
            {
                dt_medians_relevant[marker == unique_marker, mfi_all_gates := transformlist[[unique_marker]](mfi_all_gates)]
            },
            error = function(e) {
                warning(paste("Error in transformation for marker", unique_marker, ": ", e$message))
            }
        )
    }
    p0 <- ggplot2::ggplot(
        dt_medians_relevant, ggplot2::aes(x = Time, y = mfi_all_gates, color = Device)
    ) +
        ggplot2::geom_point() + # Scatter plot of data points
        ggpubr::theme_pubclean() + # Clean publication-ready theme
        ggplot2::facet_wrap(~marker, scales = "free_y") + # Facet by population
        # Linear regression with confidence interval
        ggplot2::geom_smooth(formula = y ~ x, method = "loess", se = TRUE, alpha = .2) +
        ggplot2::ylab("Transformed MFI") + # Y-axis label
        ggplot2::xlab("Time") # X-axis label

    if (!all(is.null(device_colors))) {
        p0 <- p0 +
            ggplot2::scale_color_manual(values = device_colors) + # Custom device colors for lines
            ggplot2::scale_fill_manual(values = device_colors) # Custom device colors for fill
    }
    p0
}
