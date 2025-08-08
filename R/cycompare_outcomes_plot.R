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
    p_lineplot_icc <- plot_lineplot_icc(
        preds_bound = preds_bound,
        dfcol_grouping_supersamples = dfcol_grouping_supersamples,
        dfcol_grouping_samples = dfcol_grouping_samples,
        names_grouping_train = names_grouping_train,
        dfcol_train_validation_other = dfcol_train_validation_other,
        dv_class_positive = dv_class_positive,
        device_colors = device_colors
    )
}
