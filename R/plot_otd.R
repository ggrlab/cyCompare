#' Plot Optimal Transport Distance (OTD) to a Reference Sample Over Time
#'
#' This function computes the Optimal Transport Distance (OTD) from each flow cytometry sample
#' to a shared reference sample and visualizes the results over time, optionally grouped by device.
#' The distance metric and calculation behavior can be customized through `kwargs_loss`.
#'
#' @param ff_gated A list of `flowFrame` objects representing gated flow cytometry samples.
#' @param df A `data.table` containing metadata with at least `"File"`, `"Device"`, and `"Time"` columns.
#' @param device_colors A named vector of colors for each device. Used to color the plot.
#' @param transformlist A single function or a list of transformation functions for transforming markers.
#' @param n_mastersample Integer. Number of events to sample from the combined data to create the master sample.
#' @param kwargs_loss A named list of arguments passed to `loss_pairwise()`. Should include at least a `loss` function.
#' @param relevant_columns Character vector specifying which columns (markers) to use. If `NULL`, all columns are used.
#'
#' @return A `ggplot2` object showing OTD values to the master sample over time, optionally stratified by device.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' fs <- cytobench::simulate_fs(n_samples = 5, ncells = 300, columns = c("CD4", "CD8"))
#'
#' df_meta <- data.table::data.table(
#'     File = flowCore::sampleNames(fs),
#'     Device = c("A", "A", "B", "B", "A"),
#'     Time = c(1, 2, 3, 4, 5)
#' )
#'
#' device_cols <- c("A" = "steelblue", "B" = "firebrick")
#'
#' p <- plot_otd(
#'     ff_gated = fs,
#'     df = df_meta,
#'     device_colors = device_cols,
#'     transformlist = function(x) asinh(x / 1000),
#'     n_mastersample = 1000
#' )
#'
#' print(p)
#' }
plot_otd <- function(ff_gated,
                     df,
                     device_colors = NULL,
                     transformlist,
                     n_mastersample = 1e4,
                     kwargs_loss = list(
                         loss = lossfun_hist,
                         verbose = FALSE,
                         write_intermediate = FALSE,
                         should_skip = function(i, j) FALSE,
                         take_time = FALSE,
                         return_as_matrix = TRUE
                     ),
                     relevant_columns = NULL) {
    # Default: use all columns from first flowFrame
    if (all(is.null(relevant_columns))) {
        relevant_columns <- flowCore::colnames(ff_gated[[1]])
    }

    # Convert list of flowFrames to a flowSet
    gated_fs <- flowCore::flowSet(ff_gated)

    # Build transformList: apply single function to all columns if needed
    fc_transformlist <- flowCore::transformList(
        relevant_columns,
        if (length(transformlist) == 1) rep(transformlist, length(relevant_columns)) else transformlist
    )

    # Apply transformation
    gated_fs_transformed <- flowCore::transform(gated_fs, fc_transformlist)

    # Convert each flowFrame to a data.table and keep track of sample name
    dt_events <- flowCore::fsApply(
        gated_fs_transformed,
        simplify = FALSE,
        function(x) {
            flowCore::exprs(x) |> data.table::as.data.table()
        }
    )

    # Combine all data into one table and draw a random master sample
    dt_events_bound <- data.table::rbindlist(dt_events, idcol = "File")
    mastersample <- dt_events_bound[sample(.N, n_mastersample), ][, File := NULL]

    # Compute pairwise distances from master sample to each file
    loss_calculated <- do.call(loss_pairwise, c(
        list(
            datalist_A = list("mastersample" = mastersample),
            datalist_B = dt_events
        ), kwargs_loss
    ))

    # Extract distance matrix (1 row: master sample vs all samples)
    loss_calc_mat <- loss_calculated[["dist"]]

    # Convert distances to long-format data.frame
    dists_df <- data.frame(t(loss_calc_mat)) |>
        tibble::rownames_to_column("File") |>
        dplyr::rename("OTD\nto mastersample" = mastersample)

    # Merge with metadata
    df_with_dists <- dplyr::left_join(df, dists_df, by = "File")

    # Create ggplot showing OTD over time
    p0 <- ggplot2::ggplot(
        df_with_dists,
        ggplot2::aes(x = Time, y = `OTD\nto mastersample`, color = Device, fill = Device)
    ) +
        ggplot2::geom_point() + # Raw data points
        ggpubr::theme_pubr() +
        ggplot2::theme(
            axis.line.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank()
        ) +
        ggplot2::geom_smooth(formula = y ~ x, method = "loess", se = TRUE, alpha = .2) +
        ggplot2::xlab("Time")

    # Apply manual device colors if given
    if (!all(is.null(device_colors))) {
        p0 <- p0 +
            ggplot2::scale_color_manual(values = device_colors) +
            ggplot2::scale_fill_manual(values = device_colors)
    }

    return(p0)
}
