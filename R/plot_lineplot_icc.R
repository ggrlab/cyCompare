#' Plot Predictions and Compute ICC Across Devices
#'
#' For each grouping combination, this function:
#' - Extracts predicted probabilities for the positive class,
#' - Plots these per sample (faceted by outcome and grouped by device),
#' - Calculates ICC (intra-class correlation) across devices.
#'
#' @param preds_bound Data frame of predictions in long format.
#' @param dfcol_grouping_supersamples Character vector of grouping columns (e.g., "Study").
#' @param dfcol_grouping_samples Character scalar: column representing the device (e.g., "Device").
#' @param dfcol_train_validation_other Character scalar indicating the train/validation/test split.
#' @param dv_class_positive Named vector mapping outcomes to their positive class labels.
#' @param device_colors Named vector assigning colors to devices.
#' @param names_grouping_train
#' Character vector of grouping columns used for training (e.g., "Study", "Device").
#' For each training-grouping, the function computes the ICC and plots predictions.
#'
#' @return A list with:
#'   \item{icc_summary}{ICC statistics per grouping.}
#'   \item{prediction_plots}{List of ggplot objects per grouping.}
#'   \item{plot_icc}{A ggplot object visualizing ICC values across groupings.}
#' @export
plot_lineplot_icc <- function(preds_bound,
                              names_grouping_train,
                              dfcol_grouping_samples,
                              dfcol_grouping_supersamples,
                              dfcol_train_validation_other,
                              dv_class_positive,
                              device_colors) {
    xaxis_current <- ICC2k_ICC <- outcome_ <- NULL # For R CMD check compatibility
    icc_summary_list <- list()
    groupings <- dplyr::group_by(
        preds_bound,
        !!!rlang::syms(c(
            names_grouping_train,
            "outcome_",
            "model_",
            "clustering_",
            dfcol_grouping_supersamples
            ## Note: Do NOT group by dfcol_grouping_samples here,
            ## we want to compare across devices:
            # dfcol_grouping_samples
        ))
    )

    prediction_plot_list <- dplyr::group_map(groupings, ~ {
        df_group <- .x
        metadata <- .y
        current_outcome <- metadata$outcome_
        positive_label <- dv_class_positive[[current_outcome]]
        positive_label_name <- paste0("prob.", positive_label)

        # Subset to positive-class predictions only
        preds_pos <- df_group |>
            dplyr::filter(prediction_name == positive_label_name) |>
            # Unite the sample grouping columns.
            tidyr::unite(
                col = "grouping_color",
                c(dfcol_grouping_samples),
                remove = FALSE
            )

        # Create wide format for ICC
        per_supersample_preds_wide <- preds_pos |>
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
        df_devices <- dplyr::select(per_supersample_preds_wide, dplyr::starts_with("predXX_"))

        if (ncol(df_devices) > 1) {
            res <- psych::ICC(df_devices)[["results"]]
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
        icc_summary_list[[length(icc_summary_list) + 1]] <<- res_icc

        # Order samples by training device prediction
        train_dev <- metadata[[paste0("train.", dfcol_grouping_samples)]]
        train_col <- paste0("predXX_", train_dev)
        if (train_col %in% names(per_supersample_preds_wide)) {
            supersamples_order_by_train_device <- per_supersample_preds_wide |>
                dplyr::arrange(!!rlang::sym(train_col)) |>
                dplyr::pull(SuperSample)
        } else {
            warning(paste0(
                "Column ", train_col, " not found in predictions.
                        Did you measure the sample on this device? Ordering by values directly."
            ))
            supersamples_order_by_train_device <- preds_pos |>
                dplyr::arrange(prediction_value) |>
                dplyr::pull(SuperSample) |>
                unique()
        }

        # Plot predictions per sample
        plot <- ggplot2::ggplot(
            preds_pos,
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
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)) +
            ggplot2::labs(
                shape = current_outcome,
                color = paste(dfcol_grouping_samples, collapse = "-")
            )
    })


    icc_summary <- dplyr::bind_rows(icc_summary_list) |> tibble::remove_rownames()

    p_icc <- ggplot2::ggplot(
        icc_summary |>
            dplyr::rowwise() |>
            dplyr::mutate(
                xaxis_current = do.call(paste, c(dplyr::across(tidyr::all_of(names_grouping_train)), sep = "-"))
            ),
        ggplot2::aes(x = xaxis_current, y = ICC2k_ICC, shape = outcome_)
    ) +
        ggplot2::geom_point() +
        ggpubr::theme_pubclean()

    list(
        plot_icc = p_icc,
        plot_predictions = prediction_plot_list,
        icc_summary = icc_summary
    )
}
