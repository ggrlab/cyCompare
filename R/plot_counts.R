#' Plot Population Counts Over Time
#'
#' This function plots the count and percentage of each population relative to the root population
#' over time, grouped by device. It supports optional filtering of specific populations and
#' custom device colors.
#' Data is grouped by device and faceted by population.
#'
#' @param dt_counts A `data.table` with population counts. Must contain columns:
#'   `File`, `pop`, `count`, `Time`, and `Device`.
#' @param populations Optional character vector specifying which populations to plot.
#'   If `NULL`, all populations present in `dt_counts$pop` will be used.
#' @param device_colors Optional named vector defining plot colors per device
#'   (e.g., `c("DeviceA" = "steelblue", "DeviceB" = "firebrick")`).
#'
#' @return A named list of two `ggplot2` plot objects:
#' \describe{
#'   \item{`count`}{Plot of raw event counts over time.}
#'   \item{`percentage`}{Plot of population proportion relative to root population.}
#' }
#'
#' @examples
#' # Simulate example data
#' set.seed(123)
#' dt_counts <- data.table::data.table(
#'     File = rep(paste0("Sample", 1:6), each = 3),
#'     pop = rep(c("root", "CD4+", "CD8+"), times = 6),
#'     count = rpois(18, lambda = 1000),
#'     Time = rep(rep(1:3, each = 3), 2),
#'     Device = rep(c("A", "B"), each = 9)
#' )
#'
#' # Optional: specify populations and custom device colors
#' selected_pops <- c("CD4+", "CD8+")
#' device_cols <- c("A" = "blue", "B" = "green")
#'
#' # Generate plots
#' plots <- plot_counts(dt_counts, populations = selected_pops, device_colors = device_cols)
#'
#' # Display the plots
#' plots$count # Raw count over time
#' plots$percentage # Percentage of root over time

#'
#' @export
plot_counts <- function(dt_counts, populations = NULL, device_colors = NULL) {
    pop <- Percentage_total <- count <- root_count <- Time <- Device <- NULL # For data.table compatibility
    # Step 1: Fallback to all populations if none specified
    if (all(is.null(populations))) {
        # Todo: Before, this was dt_counts$name, why?
        browser()
        populations <- unique(dt_counts$pop)
    }

    # Step 2: Extract root counts for normalization
    dt_root <- dt_counts[pop == "root"][, c("File", "count"), with = FALSE]
    data.table::setnames(dt_root, "count", "root_count")

    # Step 3: Filter for selected populations
    dt_pop <- dt_counts[pop %in% populations]

    # Step 4: Merge root counts into population table by `File`
    dt_pop_root <- dt_pop[dt_root, on = "File"]

    # Step 5: Compute population percentage relative to root
    dt_pop_root[, Percentage_total := count / root_count]
    # Rename `count` column for clarity
    data.table::setnames(dt_pop_root, "count", "Count")

    # Step 6: Generate plots: one for raw count, one for percentage
    metric_names <- c("Count", "Percentage_total")

    plots <- sapply(
        metric_names,
        simplify = FALSE,
        function(metric) {
            p <- ggplot2::ggplot(
                dt_pop_root,
                ggplot2::aes(x = Time, y = !!ggplot2::sym(metric), col = Device, fill = Device)
            ) +
                ggplot2::geom_point() +
                ggpubr::theme_pubclean() +
                ggplot2::facet_wrap(~pop, scales = "free_y") +
                # LOESS with confidence interval
                ggplot2::geom_smooth(formula = y ~ x, method = "loess", se = TRUE, alpha = 0.2)

            # Optional: add device-specific colors
            if (!all(is.null(device_colors))) {
                p <- p +
                    # Custom device colors for lines
                    ggplot2::scale_color_manual(values = device_colors) +
                    # Custom device colors for confidence intervals (fill)
                    ggplot2::scale_fill_manual(values = device_colors)
            }

            return(p)
        }
    )

    return(plots)
}
