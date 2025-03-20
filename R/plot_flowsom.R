#' Generate FlowSOM Clustering and Visualization for Flow Cytometry Data
#'
#' This function performs FlowSOM clustering on gated flow cytometry data and
#' generates PCA and MA plots for visualization.
#'
#' @param ff_gated A list of `flowFrame` objects containing gated flow cytometry data.
#' @param df A `data.table` containing metadata with at least the columns `"File"`, `"Device"`, and `"Sample"`.
#' @param device_colors A named vector of colors for each device to use in the plots.
#' @param transformlist A list of transformation functions for each marker. If a single function is provided, it will be applied to all markers.
#' @param nClus An integer specifying the number of clusters for FlowSOM.
#' @param scale A logical indicating whether to scale the data in FlowSOM clustering .
#' @param xdim An integer specifying the x-dimension of the FlowSOM grid.
#' @param ydim An integer specifying the y-dimension of the FlowSOM grid.
#' @param seed An integer specifying the random seed for FlowSOM clustering .
#' @param ... Additional parameters passed to `FlowSOM::FlowSOM()`.
#'
#' @return A list containing:
#'   \item{plots_pca}{A list of PCA plots for each clustering result, colored by device.}
#'   \item{p_MA}{A list of MA plots comparing cluster proportions between devices.}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' plot_results <- plot_flowsom(
#'     ff_gated = my_ff_list,
#'     df = metadata_dt,
#'     device_colors = my_colors,
#'     transformlist = my_transform_functions
#' )
#'
#' # Access PCA plots
#' plot_results$plots_pca[[1]]
#'
#' # Access MA plots
#' plot_results$p_MA[[1]]
#' }
#'
plot_flowsom <- function(ff_gated,
                         df,
                         device_colors,
                         transformlist,
                         nClus = 5,
                         scale = FALSE,
                         xdim = 3,
                         ydim = 3,
                         seed = 3711283, ...) {
    # Convert the list of flowFrames into a flowSet
    gated_fs <- flowCore::flowSet(ff_gated)

    # If a single transformation function is provided, apply it to all markers
    if (length(transformlist) == 1) {
        fc_transformlist <- flowCore::transformList(
            flowCore::colnames(gated_fs[[1]]),
            transformlist
        )
    } else {
        fc_transformlist <- flowCore::transformList(
            flowCore::colnames(gated_fs[[1]]),
            transformlist
        )
    }

    # Apply transformation to the flowSet
    gated_fs_transformed <- flowCore::transform(gated_fs, fc_transformlist)

    # Perform FlowSOM clustering on the transformed data
    flowsom_all <- FlowSOM::FlowSOM(
        input = gated_fs_transformed,
        transform = FALSE, # Data is already transformed
        transformList = NULL,
        nClus = nClus, # Number of clusters
        scale = scale, # Whether to scale the data
        xdim = xdim, # Grid x-dimension for SOM
        ydim = ydim, # Grid y-dimension for SOM
        seed = seed, # Set random seed for reproducibility
        ...
    )

    # Predict cluster assignments for the input data
    fs_pred <- cytobench::flowSOM_predict(flowsom_all, gated_fs_transformed)

    # Generate PCA plots for FlowSOM clusters
    plots_pca <- lapply(names(fs_pred[["ncells_per_x"]]), function(x) {
        x_data <- fs_pred[["ncells_per_x"]][[x]]

        # Perform PCA on cluster abundance data
        res_pca <- stats::prcomp(x_data |> dplyr::select(-sample), scale = TRUE)

        # Generate PCA plot using ggfortify
        ggfortify:::autoplot.prcomp(
            res_pca,
            data = df,
            colour = "Device"
        ) +
            ggpubr::theme_pubr() +
            ggplot2::theme(legend.position = "top") +
            ggplot2::labs(
                title = "PCA of FlowSOM",
                subtitle = x
            )
    })

    ## Generate MA plots to compare cluster proportions between devices
    p_MA <- lapply(names(fs_pred[["ncells_per_x"]]), function(x) {
        x_data <- fs_pred[["ncells_per_x"]][[x]]
        x_data_numeric <- x_data |> dplyr::select(-sample)

        # Compute cluster proportions per sample
        proportions <- x_data_numeric / rowSums(x_data_numeric)
        x_data[, -1] <- proportions

        # Merge with metadata
        x_data_df <- dplyr::left_join(x_data, df, by = c("sample" = "File"))

        # Convert data to long format for visualization
        x_data_df_long <- x_data_df |>
            tidyr::pivot_longer(
                cols = tidyr::all_of(grep(colnames(x_data_df), pattern = "[cC]luster", value = TRUE)),
                names_to = "cluster_id",
                values_to = "proportion"
            )

        # Reshape data to compare proportions per device
        x_data_df_persample <- x_data_df_long |>
            dplyr::select(Device, SuperSample, Sample, cluster_id, proportion) |>
            tidyr::pivot_wider(
                names_from = c("Device"),
                values_from = "proportion"
            )

        # Generate all pairwise comparisons of devices
        device_combinations <- combn(names(device_colors), 2)
        plotlist <- list()

        for (combination_i in seq_len(ncol(device_combinations))) {
            device_combination <- device_combinations[, combination_i]
            d1 <- device_combination[1]
            d2 <- device_combination[2]

            # Create an MA plot (log-ratio vs mean) for cluster proportions
            p0 <- ggplot2::ggplot(
                x_data_df_persample,
                ggplot2::aes(
                    y = log2(!!ggplot2::sym(d1) / !!ggplot2::sym(d2)), # Log fold-change
                    x = ((!!ggplot2::sym(d1) + !!ggplot2::sym(d2)) / 2) # Mean proportion
                )
            ) +
                ggplot2::ggtitle(paste0(d1, " vs ", d2)) +
                ggplot2::geom_abline(intercept = 0, slope = 0) + # Reference line at log2 fold-change = 0
                ggplot2::theme_bw() +
                ggplot2::scale_x_log10() + # Log-scale for x-axis
                ggplot2::theme(
                    axis.text.x = ggplot2::element_text(size = 6),
                    legend.position = "right",
                    legend.key.size = ggplot2::unit(.5, "cm"),
                    legend.title = element_text(angle = -90)
                )

            # Add fold-change reference lines (x2, x10, x25)
            for (specific_lines in c(2, 10, 25)) {
                p0 <- p0 +
                    ggplot2::geom_hline(
                        yintercept = c(-log2(specific_lines), log2(specific_lines)),
                        linetype = "dashed"
                    ) +
                    ggplot2::annotate(
                        "text",
                        x = -Inf, y = log2(specific_lines),
                        label = paste0("x", specific_lines),
                        hjust = 0, vjust = 0, size = 2.5
                    )
            }

            # Overlay density information using a 2D histogram
            plotlist[[combination_i]] <- p0 +
                ggplot2::stat_bin_2d(
                    ggplot2::aes(fill = log10(ggplot2::after_stat(count))), # Log-transformed count for coloring
                    bins = 100,
                    geom = "raster"
                ) +
                ggplot2::labs(
                    fill = "Log2(n cells in all clusters and samples)"
                ) +
                ggplot2::scale_fill_viridis_c(na.value = NA)
        }
    })

    # Return both PCA and MA plots as a list
    return(
        list(
            plots_pca = plots_pca, # List of PCA plots
            p_MA = p_MA # List of MA plots
        )
    )
}
