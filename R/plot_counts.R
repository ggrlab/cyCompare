#' Plot Population Counts Over Time
#'
#' This function plots the count and percentage of each population relative to the root population
#' over time, grouped by device. It supports optional filtering of specific populations and
#' custom device colors.
#'
#' @param dt_counts A `data.table` containing population counts with columns: `File`, `pop`, `count`, `Time`, and `Device`.
#' @param populations A character vector specifying which populations to plot. If `NULL`, all populations are used.
#' @param device_colors A named vector specifying custom colors for devices (e.g., `c("Device1" = "red", "Device2" = "blue")`).
#'
#' @return A list of two `ggplot2` objects: one for raw counts and another for percentage relative to the root population.
#' @export

plot_counts <- function(dt_counts, populations = NULL, device_colors = NULL) {
    # If populations are not provided, use all unique population names
    if (all(is.null(populations))) {
        populations <- unique(dt_counts$name)
    }

    # Extract root population counts for normalization
    dt_root <- dt_counts[pop == "root"][, c("File", "count"), with = FALSE]
    data.table::setnames(dt_root, "count", "root_count")

    # Subset to selected populations
    dt_pop <- dt_counts[pop %in% populations]

    # Merge population data with root counts based on `File`
    dt_pop_root <- dt_pop[dt_root, on = "File"]

    # Compute the percentage relative to root count
    dt_pop_root[, Percentage_total := count / root_count]
    data.table::setnames(dt_pop_root, "count", "Count") # Rename `count` column for clarity

    # Generate plots for both raw counts and percentage relative to root
    plots <- lapply(
        c("Count", "Percentage_total"),
        function(x) {
            p0 <- ggplot2::ggplot(
                dt_pop_root,
                ggplot2::aes(x = Time, y = !!ggplot2::sym(x), col = Device, fill = Device)
            ) +
                ggplot2::geom_point() + # Scatter plot of data points
                ggpubr::theme_pubr() + # Clean publication-ready theme
                ggplot2::facet_wrap(~pop, scales = "free_y") + # Facet by population
                # Linear regression with confidence interval
                ggplot2::geom_smooth(method = "loess", se = TRUE, alpha = .2)
            if (!all(is.null(device_colors))) {
                p0 <- p0 +
                    ggplot2::scale_color_manual(values = device_colors) + # Custom device colors for lines
                    ggplot2::scale_fill_manual(values = device_colors) # Custom device colors for fill
            }
            p0
        }
    )

    return(plots) # Return the list of plots
}
