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
            p0 <- ggplot2::ggplot(
                x_data_df_persample,
                ggplot2::aes(
                    y = log2(!!ggplot2::sym(d1) / !!ggplot2::sym(d2)),
                    x = ((!!ggplot2::sym(d1) + !!ggplot2::sym(d2)) / 2)
                )
            ) +
                ggplot2::ggtitle(paste0(d1, " vs ", d2), subtitle = x) +
                ggplot2::geom_abline(intercept = 0, slope = 0) +
                ggplot2::theme_bw() +
                ggplot2::scale_x_log10() +
                ggplot2::theme(
                    axis.text.x = ggplot2::element_text(size = 6),
                    legend.position = "right",
                    legend.key.size = ggplot2::unit(.5, "cm"),
                    legend.title = ggplot2::element_text(angle = -90)
                )
            minimum_x_value <- min(
                x_data_df_persample[[d1]],
                x_data_df_persample[[d2]],
                na.rm = TRUE
            )
            # Facet by train/test if applicable
            if (!all(is.null(dfcol_train_validation_other))) {
                p0 <- p0 + ggplot2::facet_wrap(dfcol_train_validation_other)
            }

            # Add horizontal reference lines for fold-changes (e.g., x2, x10, x25)
            for (specific_lines in MA_horizontal_lines_FC) {
                p0 <- p0 +
                    ggplot2::geom_hline(
                        yintercept = c(-log2(specific_lines), log2(specific_lines)),
                        linetype = "dashed"
                    ) +
                    ggplot2::annotate(
                        "text",
                        x = minimum_x_value, y = log2(specific_lines),
                        label = paste0("x", specific_lines),
                        hjust = 0, vjust = 0, size = 5 # 2.5
                    )
            }

            # Add density information via 2D binning
            plotlist[[combination_i]] <- p0 +
                ggplot2::stat_bin_2d(
                    ggplot2::aes(fill = log10(ggplot2::after_stat(count) + 1)),
                    bins = MA_bins,
                    geom = "tile"
                ) +
                ggplot2::labs(fill = "Log10(n cells +1)") +
                ggplot2::scale_fill_viridis_c(na.value = NA)
        }

        return(plotlist)
    })
}
