#' Calculate ROC Metrics Grouped by Metadata
#'
#' Computes ROC (Receiver Operating Characteristic) statistics for predicted outcomes,
#' grouped by relevant metadata such as device, model, and outcome type.
#' The positive class label is dynamically determined for each outcome using the `dv_class_positive` mapping.
#'
#' @inheritParams cycompare_outcomes_analyse
#' @param result_grouping A data frame containing prediction results and metadata columns used for grouping
#' (e.g., `outcome_`, `model_`, `clustering_`, and those specified by `dfcol_*`).
#' @param ... Currently unused; included for extensibility.
#'
#' @return A tibble containing grouped ROC statistics for each group defined by metadata. Includes:
#' \itemize{
#'   \item ROC coordinates at each group's validation threshold.
#'   \item AUC and confidence intervals.
#'   \item Confusion matrix-based metrics (accuracy, kappa, null accuracy, etc.).
#'   \item Log loss and class frequencies.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming `result_df` is pre-processed with the appropriate columns
#' rocs <- calc_ROC_grouped(
#'     df = NULL,
#'     result_grouping = result_df,
#'     dfcol_grouping_supersamples = "Study",
#'     dfcol_grouping_samples = "Device",
#'     dfcol_outcomes = c("outcome_1", "outcome_2"),
#'     dfcol_train_validation_other = "train_validation_test",
#'     dv_class_positive = c("outcome_1" = "A", "outcome_2" = 5.1)
#' )
#' }
calc_ROC_grouped <- function(
    df = NULL, # Placeholder argument for compatibility
    result_grouping,
    dfcol_grouping_supersamples = c("Study"),
    dfcol_grouping_samples = "Device",
    dfcol_outcomes = c("outcome_1", "outcome_2"),
    dfcol_train_validation_other = "train_validation_test",
    dv_class_positive = c("outcome_1" = "A", "outcome_2" = 5.1),
    ...) {
    # Group the results by metadata columns to compute ROC per group
    grouped_rocs <- result_grouping |>
        dplyr::group_by(
            !!!rlang::syms(c(
                "outcome_",
                "model_",
                "clustering_",
                dfcol_grouping_supersamples,
                dfcol_grouping_samples,
                dfcol_train_validation_other
            ))
        ) |>
        dplyr::group_map(
            ~ {
                data <- .x
                current_group <- .y

                # Determine the correct positive class for this outcome
                pos_label <- dv_class_positive[[current_group[["outcome_"]]]]

                # Define other class levels for binary ROC
                other_labels <- setdiff(
                    unique(data[[current_group[["outcome_"]]]]),
                    pos_label
                )

                # Compute ROC curve
                roc <- pROC::roc(
                    response = data[[current_group[["outcome_"]]]],
                    predictor = data[["prediction_value"]],
                    direction = "<",
                    levels = c(pos_label, other_labels)
                )

                # Bind current grouping info with ROC coordinates at the selected threshold
                cbind(
                    current_group,
                    coords_helper(roc, threshold = unique(data[["threshold_proc_validationset"]]))
                )
            }
        ) |>
        do.call(what = rbind) |> # Combine results from all groups into one data frame
        tibble::as_tibble()

    return(grouped_rocs)
}
