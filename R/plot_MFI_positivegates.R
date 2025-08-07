#' Plot Median Fluorescence Intensity (MFI) for Positive Gates
#'
#' This function plots the transformed median fluorescence intensity (MFI) values over time
#' for specific positive gates, grouped by device and faceted by marker. It supports optional
#' transformation functions, log-ratio normalization, and device-specific coloring.
#'
#' @param dt_count_mfi A `data.table` containing per-population MFI values with the following columns:
#'   \describe{
#'     \item{File}{Identifier of the FCS file}
#'     \item{pop}{Population name (e.g., "CD4+ T cells")}
#'     \item{MFI}{Median fluorescence intensity for the population}
#'     \item{Time}{Time point or acquisition date}
#'     \item{Device}{Identifier of the instrument used}
#'   }
#' @param marker_to_gate A named list mapping markers to one or more population names used to compute MFIs.
#'   For example: `list("CD4" = c("CD4 T", "CD4 Memory T"))`.
#' @param device_colors An optional named character vector mapping device names to colors,
#'   e.g., `c("Cytometer1" = "steelblue", "Cytometer2" = "firebrick")`.
#' @param transformlist An optional named list of transformation functions for each marker,
#'   e.g., `list("CD3" = log10, "CD4" = sqrt)`. If a single function or list is given,
#'   it is applied uniformly to all markers.
#' @param meanratio Logical. If `TRUE`, plots the ratio of each MFI value to the mean MFI
#'   of its marker across all samples (useful for comparing relative changes). If `FALSE`,
#'   plots the raw (possibly transformed) MFI values. Default is `FALSE`.
#'
#' @return A `ggplot2` object with one facet per marker showing MFI or MFI ratios over time.
#'   Devices are distinguished by color or fill. A LOESS smoothing curve is included to highlight trends.
#'
#' @export
plot_MFI_positivegates <- function(dt_count_mfi, marker_to_gate, device_colors, transformlist = NULL, meanratio = FALSE) {
    tmp <- marker_to_gate_count(marker_to_gate, dt_count_mfi = dt_count_mfi)
    marker_to_gate_dt <- tmp[["marker_to_gate_dt"]]
    dt_count_mfi_relevant <- tmp[["dt_count_mfi_relevant"]]

    # Merge the filtered MFI data with marker-to-gate mappings using 'pop' as the key
    dt <- marker_to_gate_dt[dt_count_mfi_relevant, on = "pop", allow.cartesian = TRUE]

    # Reshape to long format for per-marker plotting
    dt_long <- data.table::melt(
        dt,
        measure.vars = names(marker_to_gate),
        value.name = "MFI",
        variable.name = "variable"
    )

    # Compute median MFI values grouped by File, Device, Time, and marker
    dt_medians <- dt_long[, .(mfi_all_gates = stats::median(MFI)), by = .(File, Device, Time, marker, variable)]

    # Keep only rows where marker matches the variable name
    dt_medians_relevant <- dt_medians[marker == variable]

    if (meanratio) {
        # Normalize by mean MFI across all samples for each marker
        marker_mean <- dt_medians_relevant[, .(mean_mfi = mean(mfi_all_gates)), by = .(marker)]
        dt_medians_relevant_meanratio <- dt_medians_relevant[marker_mean, on = "marker"][
            ,
            mfi_by_meanMFI := (mfi_all_gates / mean_mfi)
        ]

        p0 <- ggplot2::ggplot(
            dt_medians_relevant_meanratio,
            ggplot2::aes(x = Time, y = mfi_by_meanMFI, color = Device, fill = Device)
        ) +
            ggplot2::ylab("Ratio MFI to meanMFI") +
            ggplot2::facet_wrap(~marker) +
            ggplot2::geom_hline(yintercept = 1, linetype = "dashed") +
            ggplot2::scale_y_log10()
    } else {
        transformlist <- transformlist_named(transformlist, unique(dt_medians_relevant$marker))

        # Apply per-marker transformation if available
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
            dt_medians_relevant,
            ggplot2::aes(x = Time, y = mfi_all_gates, color = Device, fill = Device)
        ) +
            ggplot2::facet_wrap(~marker, scales = "free_y") +
            ggplot2::ylab("Transformed MFI")
    }

    # Add scatter, smoothing, and final styling
    p0 <- p0 +
        ggplot2::geom_point() +
        ggpubr::theme_pubr() +
        ggplot2::theme(
            axis.line.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank()
        ) +
        ggplot2::geom_smooth(formula = y ~ x, method = "loess", se = TRUE, alpha = .2) +
        ggplot2::xlab("Time")

    # Apply custom color mappings
    if (!all(is.null(device_colors))) {
        p0 <- p0 +
            ggplot2::scale_color_manual(values = device_colors) +
            ggplot2::scale_fill_manual(values = device_colors)
    }

    return(p0)
}
