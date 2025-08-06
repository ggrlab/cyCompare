#' Calculate ROC Metrics Grouped by Metadata
#'
#' Computes ROC (Receiver Operating Characteristic) statistics for prediction results, grouped by specified metadata columns.
#' The positive class label is dynamically determined per outcome using the `dv_class_positive` mapping.
#'
#' @param df Data frame input (not used but kept for compatibility).
#' @param result_grouping A data frame containing prediction results and relevant metadata columns for grouping.
#' @param dfcol_grouping_supersamples Character vector of column names used for higher-level grouping (e.g., "Study").
#' @param dfcol_grouping_samples Character, name of the column used for sample-level grouping (e.g., "Device").
#' @param dfcol_outcomes Character vector of column names representing different outcomes (e.g., `c("outcome_1", "outcome_2")`).
#' @param dfcol_train_validation_other Character, column name that differentiates between train/validation/test sets.
#' @param dv_class_positive Named vector specifying the positive class for each outcome (e.g., `c("outcome_1" = "A", "outcome_2" = 5.1)`).
#' @param ... Additional arguments passed to internal functions (currently unused).
#'
#' @return A tibble containing grouped ROC metrics, including coordinates extracted at specified thresholds.
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
