#' Generate Outcome and ROC Plots from cycompare Analysis
#'
#' Produces diagnostic plots from a full `cycompare_outcomes_analyse()` run.
#' This includes outcome distributions, ROC curves, and predicted probabilities grouped by supersample and device.
#' ICC (intra-class correlation) metrics are calculated for each outcome across devices.
#'
#' @inheritParams cycompare_outcomes_analyse
#' @param results_cycompare_analyse
#' A named list returned by `cycompare_outcomes_analyse()`,
#'   containing components `prepared_data`, `clustering`, and `models`.
#'
#' @return Invisibly returns a list of generated plots. Side effects include:
#' \itemize{
#'   \item Outcome visualizations grouped by metadata.
#'   \item ICC plots comparing prediction consistency across devices.
#' }
#'
#' @export
#'
cycompare_outcomes_plot <- function(
    df,
    dfcol_grouping_supersamples = c("Study"),
    dfcol_grouping_samples = "Device",
    dfcol_outcomes = c("outcome_1", "outcome_2"),
    dfcol_train_validation_other = "train_validation_test",
    dv_class_positive = c("outcome_1" = "A", "outcome_2" = 5.1),
    results_cycompare_analyse) {
    device_colors <- results_cycompare_analyse[["prepared_data"]][["device_colors"]]
    pred_applied <- results_cycompare_analyse[["models"]][["applied"]]

    # Combine predictions with corresponding group metadata
    preds_bound_list <- list()
    for (group_i in seq_len(nrow(pred_applied[["groups"]]))) {
        current_grouping <- pred_applied[["groups"]][group_i, ]
        colnames(current_grouping) <- paste0("train.", colnames(current_grouping))
        preds_bound_list[[group_i]] <- cbind(current_grouping, pred_applied[["results"]][[group_i]]) |>
            tibble::as_tibble()
    }
    names_grouping_train <- colnames(current_grouping)
    preds_bound <- dplyr::bind_rows(preds_bound_list)

    #### 1. Plot outcome distribution per study
    p1.1 <- plot_outcome_circles(
        df,
        dfcol_grouping_supersamples = dfcol_grouping_supersamples,
        dfcol_outcomes = dfcol_outcomes
    )

    # Increase memory threshold for future package (if used in downstream async calls)
    options(future.globals.maxSize = 10 * 1024^3)

    #### 2. Calculate ROC grouped by metadata
    calced_roc <- fun_grouped_apply(
        data = NULL,
        result_grouping = results_cycompare_analyse[["models"]][["applied"]],
        make_flowset = FALSE,
        fun = calc_ROC_grouped,
        outdir_base = NULL,
        verbose = FALSE,
        return_results = TRUE,
        dfcol_grouping_supersamples = dfcol_grouping_supersamples,
        dfcol_grouping_samples = dfcol_grouping_samples,
        dfcol_outcomes = dfcol_outcomes,
        dfcol_train_validation_other = dfcol_train_validation_other,
        dv_class_positive = dv_class_positive,
        bygroup = TRUE
    )

    # Plot ROC curves grouped by outcome and study/device
    rocs <- plot_ROC_grouped(
        calced_roc,
        dfcol_grouping_supersamples = dfcol_grouping_supersamples,
        dfcol_grouping_samples = dfcol_grouping_samples,
        dfcol_outcomes = dfcol_outcomes,
        dfcol_train_validation_other = dfcol_train_validation_other,
        dv_class_positive = dv_class_positive,
        device_colors = device_colors
    )

    # # Quick debug: save first ROC plot to PDF (can be removed in production)
    # pdf("removeme.pdf")
    # rocs[["roc.plot"]]
    # dev.off()

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
                res <- psych::ICC(
                    per_supersample_preds_wide |>
                        dplyr::select(dplyr::starts_with("predXX_"))
                )[["results"]]

                # Flatten and relabel results for summary
                res_oneline <- data.frame(
                    res[1, ], res[2, ], res[3, ],
                    res[4, ], res[5, ], res[6, ]
                )
                names(res_oneline) <- paste0(rep(res[["type"]], each = ncol(res)), "_", colnames(res))

                res_icc <- cbind(.y, res_oneline)

                # Order supersamples by predicted values from training device
                supersamples_order_by_train_device <- per_supersample_preds_wide |>
                    dplyr::arrange(
                        !!rlang::sym(paste0("predXX_", current_grouping[["train.Device"]]))
                    ) |>
                    dplyr::pull("SuperSample")

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
                        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
                    ) +
                    ggplot2::labs(
                        shape = current_outcome,
                        color = paste0(c(dfcol_grouping_samples), collapse = "-")
                    )
            }
        )
}
