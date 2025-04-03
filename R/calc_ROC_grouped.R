calc_ROC_grouped <- function(
    df = NULL, # df is NOT used here, but is required for the function signature
    result_grouping,
    dfcol_grouping_supersamples = c("Study"),
    dfcol_grouping_samples = "Device",
    dfcol_outcomes = c("outcome_1", "outcome_2"),
    dfcol_train_validation_other = "train_validation_test",
    dv_class_positive = c("outcome_1" = "A", "outcome_2" = 5.1),
    ...) {
    # grouped_rocs <- df_in_group |>
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
                pos_label <- dv_class_positive[[current_group[["outcome_"]]]]
                other_labels <- setdiff(
                    unique(data[[current_group[["outcome_"]]]]),
                    pos_label
                )
                roc <- pROC::roc(
                    response = data[[current_group[["outcome_"]]]],
                    predictor = data[["prediction_value"]],
                    direction = "<",
                    levels = c(pos_label, other_labels)
                )
                return(
                    cbind(current_group, coords_helper(roc, threshold = unique(data[["threshold_proc_validationset"]])))
                )
            }
        ) |>
        do.call(what = rbind) |>
        tibble::as_tibble()

    return(grouped_rocs)
}
