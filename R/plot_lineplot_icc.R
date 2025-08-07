plot_lineplot_icc <- function(preds_bound,
                              dfcol_grouping_supersamples,
                              dfcol_grouping_samples,
                              dfcol_outcomes,
                              dfcol_train_validation_other,
                              dv_class_positive,
                              device_colors) {
    ## Plot predicted probabilities/scores per sample and device
    dplyr::group_by(
        preds_bound,
        !!!rlang::syms(
            c(
                names_grouping_train,
                "outcome_",
                "model_",
                "clustering_",
                dfcol_grouping_supersamples
                ## Note: Do NOT group by dfcol_grouping_samples here,
                ## we want to compare across devices:
                # dfcol_grouping_samples
            )
        )
    ) |>
        dplyr::group_map(
            ~ {
                current_outcome <- .y[["outcome_"]]
                positive_label <- dv_class_positive[[current_outcome]]
                positive_label_name <- paste0("prob.", positive_label)

                # Filter to predicted positive probabilities only
                preds_positive_outcome <- .x |>
                    dplyr::filter(prediction_name == positive_label_name) |>
                    # Unite the sample grouping columns.
                    tidyr::unite(
                        col = "grouping_color",
                        c(dfcol_grouping_samples),
                        remove = FALSE
                    )

                preds_positive_outcome[["prediction_value"]] <- preds_positive_outcome[["prediction_value"]] * rnorm(
                    nrow(preds_positive_outcome),
                    mean = 5,
                    sd = 1
                )

                # Prepare data for ICC calculation
                per_supersample_preds_wide <- preds_positive_outcome |>
                    dplyr::select(
                        "SuperSample",
                        dfcol_grouping_samples,
                        dfcol_train_validation_other,
                        current_outcome,
                        "prediction_value"
                    ) |>
                    tidyr::pivot_wider(
                        names_from = dfcol_grouping_samples,
                        values_from = prediction_value,
                        id_cols = c("SuperSample", dfcol_train_validation_other),
                        names_prefix = "predXX_"
                    )

                # Calculate ICC (intra-class correlation) across devices
                df_multidevice <- per_supersample_preds_wide |>
                    dplyr::select(dplyr::starts_with("predXX_"))
                if (ncol(df_multidevice) > 1) {
                    res <- psych::ICC(df_multidevice)[["results"]]
                } else {
                    res <- data.frame(
                        type = c("ICC1", "ICC2", "ICC3", "ICC1k", "ICC2k", "ICC3k"),
                        ICC = NA,
                        F = NA,
                        df1 = NA,
                        df2 = NA,
                        p = NA,
                        "lower bound" = NA,
                        "upper bound" = NA,
                        row.names = c(
                            "Single_raters_absolute",
                            "Single_random_raters",
                            "Single_fixed_raters",
                            "Average_raters_absolut",
                            "Average_random_raters",
                            "Average_fixed_raters"
                        )
                    )
                }

                # Flatten and relabel results for summary
                res_oneline <- data.frame(
                    res[1, ], res[2, ], res[3, ],
                    res[4, ], res[5, ], res[6, ]
                )
                names(res_oneline) <- paste0(rep(res[["type"]], each = ncol(res)), "_", colnames(res))
                res_icc <- cbind(.y, res_oneline)

                # Order supersamples by predicted values from training device
                train_pred_column <- paste0("predXX_", current_grouping[["train.Device"]])
                if (train_pred_column %in% colnames(per_supersample_preds_wide)) {
                    per_supersample_preds_wide <- per_supersample_preds_wide |>
                        dplyr::arrange(
                            !!rlang::sym(train_pred_column)
                        )
                    supersamples_order_by_train_device <- per_supersample_preds_wide[["SuperSample"]]
                } else {
                    warning(paste0(
                        "Column ", train_pred_column, " not found in predictions.
                        Did you measure the sample on this device? Ordering by values directly."
                    ))
                    supersamples_order_by_train_device <- preds_positive_outcome |>
                        dplyr::arrange(prediction_value) |>
                        dplyr::pull(SuperSample) |>
                        unique()
                }


                # Plot predictions
                ggplot2::ggplot(
                    preds_positive_outcome,
                    ggplot2::aes(
                        x = factor(SuperSample, levels = supersamples_order_by_train_device),
                        y = prediction_value,
                        col = grouping_color,
                        shape = relevel(factor(!!rlang::sym(current_outcome)), ref = positive_label)
                    )
                ) +
                    ggplot2::geom_point() +
                    ggplot2::scale_color_manual(values = device_colors) +
                    ggpubr::theme_pubclean() +
                    ggplot2::theme(
                        axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
                    ) +
                    ggplot2::labs(
                        shape = current_outcome,
                        color = paste0(c(dfcol_grouping_samples), collapse = "-")
                    )
            }
        )
}
