#' Plot MA Plots Comparing Cluster Proportions Between Devices
#'
#' @param fs_pred
#' Output from `run_flowsom_clustering()` which comes from `cytobench::flowSOM_predict()` and
#' must contain the `ncells_per_x` list with the number of cells per cluster and sample.
#' @param df Metadata table.
#' @param device_colors Named color vector.
#' @param dfcol_grouping_samples Device column.
#' @param dfcol_train_validation_other Optional column for faceting.
#' @param MA_horizontal_lines_FC Fold-change reference lines.
#' @param MA_bins Number of bins for 2D histogram.
#'
#' @return A list of lists of MA plots.
#' @export
plot_flowsom_ma <- function(fs_pred,
                            df,
                            device_colors,
                            dfcol_grouping_samples = "Device",
                            dfcol_train_validation_other = NULL,
                            MA_horizontal_lines_FC = c(2, 10, 25),
                            MA_bins = 100) {
    lapply(names(fs_pred[["ncells_per_x"]]), function(x) {
        x_data <- fs_pred[["ncells_per_x"]][[x]]
        x_data_numeric <- dplyr::select(x_data, -sample)

        # Compute within-sample proportions
        proportions <- sweep(x_data_numeric, 1, rowSums(x_data_numeric), "/")
        x_data[, -1] <- proportions

        # Join with metadata
        x_data_df <- dplyr::left_join(x_data, df, by = c("sample" = "File"))

        # Pivot to long format
        x_data_df_long <- tidyr::pivot_longer(
            x_data_df,
            cols = tidyr::all_of(grep(colnames(x_data_df), pattern = "[cC]luster", value = TRUE)),
            names_to = "cluster_id",
            values_to = "proportion"
        )

        # Wide format for per-device comparison
        x_data_df_persample <- x_data_df_long |>
            dplyr::select(
                !!rlang::sym(dfcol_grouping_samples[[1]]),
                !!!rlang::syms(dfcol_train_validation_other),
                SuperSample,
                Sample,
                cluster_id,
                proportion
            ) |>
            tidyr::pivot_wider(
                names_from = dfcol_grouping_samples[[1]],
                values_from = "proportion"
            )

        # All pairwise device comparisons
        device_combinations <- combn(names(device_colors), 2)
        plotlist <- list()

        for (combination_i in seq_len(ncol(device_combinations))) {
            device_combination <- device_combinations[, combination_i]
            d1 <- device_combination[1]
            d2 <- device_combination[2]

            # Generate MA plot (log2 fold-change vs average)
            p <- ggplot2::ggplot(
                x_data_df_persample,
                ggplot2::aes(
                    y = log2(!!ggplot2::sym(d1) / !!ggplot2::sym(d2)),
                    x = ((!!ggplot2::sym(d1) + !!ggplot2::sym(d2)) / 2)
                )
            ) +
                ggplot2::geom_abline(intercept = 0, slope = 0) +
                ggplot2::ggtitle(paste(d1, "vs", d2), subtitle = x) +
                ggplot2::scale_x_log10() +
                ggplot2::theme_bw() +
                ggplot2::theme(
                    axis.text.x = ggplot2::element_text(size = 6),
                    legend.position = "right",
                    legend.key.size = ggplot2::unit(.5, "cm"),
                    legend.title = ggplot2::element_text(angle = -90)
                )

            min_x <- min(
                x_data_df_persample[[d1]],
                x_data_df_persample[[d2]],
                na.rm = TRUE
            )
            if (!all(is.null(dfcol_train_validation_other))) {
                p <- p + ggplot2::facet_wrap(dfcol_train_validation_other)
            }

            for (specific_lines in MA_horizontal_lines_FC) {
                p <- p +
                    ggplot2::geom_hline(yintercept = c(-log2(specific_lines), log2(specific_lines)), linetype = "dashed") +
                    ggplot2::annotate(
                        "text",
                        x = min_x, y = log2(specific_lines),
                        label = paste0("x", specific_lines), hjust = 0, vjust = 0, size = 2.5
                    )
            }


            # Add density information via 2D binning
            plotlist[[combination_i]] <- p +
                ggplot2::stat_bin_2d(
                    ggplot2::aes(fill = log10(ggplot2::after_stat(count) + 1)),
                    bins = MA_bins,
                    geom = "tile"
                ) +
                ggplot2::scale_fill_viridis_c(na.value = NA) +
                ggplot2::labs(fill = "Log10(n cells +1)")
        }

        return(plotlist)
    })
}
