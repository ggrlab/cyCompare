#' Plot Median Fluorescence Intensity (MFI) for Positive Gates
#'
#' This function plots the transformed median fluorescence intensity (MFI) for positive gates
#' over time, grouped by device. It allows for custom transformations and device color mappings.
#'
#' @param dt_count_mfi A `data.table` containing MFI values with columns: `File`, `pop`, `MFI`, `Time`, and `Device`.
#' @param marker_to_gate A named list mapping markers to their respective gating populations.
#' @param device_colors A named vector specifying custom colors for devices (e.g., `c("Device1" = "red", "Device2" = "blue")`).
#' @param transformlist An optional named list of transformation functions for each marker (e.g., `list("CD3" = log10, "CD4" = sqrt)`). If a single value or function is provided, it will be applied to all markers.
#'
#' @param meanratio A logical indicating whether to plot the ratio of MFI to the mean MFI for each marker or the transformed MFI values (default: `FALSE`).
#' @return A `ggplot2` object visualizing the transformed MFI over time for each marker.
#' @export
plot_MFI_positivegates <- function(dt_count_mfi, marker_to_gate, device_colors, transformlist = NULL, meanratio = FALSE) {

    tmp <- marker_to_gate_count(marker_to_gate, dt_count_mfi = dt_count_mfi)
    marker_to_gate_dt <- tmp[["marker_to_gate_dt"]]
    dt_count_mfi_relevant <- tmp[["dt_count_mfi_relevant"]]

    # Merge the filtered MFI data with marker-to-gate mappings using 'pop' as the key
    dt <- marker_to_gate_dt[dt_count_mfi_relevant, on = "pop", allow.cartesian = TRUE]

    # Convert MFI data to long format for plotting
    dt_long <- data.table::melt(
        dt,
        measure.vars = names(marker_to_gate),
        value.name = "MFI",
        variable.name = "variable"
    )
    # Compute median MFI values grouped by File, Device, Time, and marker
    dt_medians <- dt_long[, .(mfi_all_gates = median(MFI)), by = .(File, Device, Time, marker, variable)]

    # Keep only rows where marker matches the variable name
    dt_medians_relevant <- dt_medians[marker == variable]
    if (meanratio) {
        marker_mean <- dt_medians_relevant[, .(mean_mfi = mean(mfi_all_gates)), by = .(marker)]
        dt_medians_relevant_meanratio <- dt_medians_relevant[marker_mean, on = "marker"][, mfi_by_meanMFI := (mfi_all_gates / mean_mfi)]
        p0 <- ggplot2::ggplot(
            dt_medians_relevant_meanratio, ggplot2::aes(x = Time, y = (mfi_by_meanMFI), color = Device, fill = Device)
        ) +
            ggplot2::ylab("Ratio MFI to meanMFI") + # Y-axis label
            ggplot2::facet_wrap(~marker) + # Facet by marker, NO free scales
            ggplot2::geom_hline(yintercept = (1), linetype = "dashed") +
            ggplot2::scale_y_log10()
    } else {
        # Apply transformation functions to MFI values (if provided)
        if (length(transformlist) == 1 && !is.null(transformlist)) {
            if (is.function(transformlist)) {
                transformlist <- list(transformlist)
            }
            transformlist <- setNames(
                rep(transformlist, length(unique(dt_medians_relevant$marker))),
                unique(dt_medians_relevant$marker)
            )
        }
        for (unique_marker in unique(dt_medians_relevant$marker)) {
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
            dt_medians_relevant, ggplot2::aes(x = Time, y = mfi_all_gates, color = Device, fill = Device)
        ) +
            ggplot2::facet_wrap(~marker, scales = "free_y") + # Facet by marker
            ggplot2::ylab("Transformed MFI") # Y-axis label
    }
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
