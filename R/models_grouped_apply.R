#' Train Grouped Models Using a Classification Pipeline
#'
#' This function applies models to grouped flow cytometry data.
#' It merges phenotype proportions (from FlowSOM clustering) with metadata and applies
#' pre-trained models (from `cytobench::wrapper_count_models`) to predict outcome variables.
#'
#' Grouping is controlled by metadata columns. If `bygroup = TRUE`, the classification is
#' restricted to the subset defined by `grouping`; otherwise, all samples are used.
#'
#' @inheritParams cycompare_outcomes_analyse
#' @param outdir A character string specifying the output directory to save model results. Not used here.
#' @param result_grouping
#' A list returned from the grouping pipeline. Must include `clustering` (with proportions) and
#' `models` (trained model list).
#' @param grouping A `data.frame` with one row per sample,  that defines grouping variables for merging with `df`.
#' @param bygroup Logical;
#'   If FALSE, the function will apply each model to the entire dataset
#' (`result_grouping[["clustering"]][[counts_proportions]][[finalmodel_cluster_selection]]`).
#'   If TRUE, the `FALSE` dataset is limited to the rows mathing the `grouping` variable.
#' @param counts_proportions A character string indicating the type of clustering proportions to use
#'   (e.g., "proportions_per_x"). This should match the names in `result_grouping[["clustering"]]`.
#'   This is used to retrieve the clustering proportions for each sample.
#'   Default is "proportions_per_x".
#'   This should match the names in `result_grouping[["clustering"]]`.
#' @param ... Not used.
#'
#' @return A `data.table` in long format with predictions for each sample, outcome, model, and clustering:
#' \describe{
#'   \item{outcome_}{The modeled outcome.}
#'   \item{model_}{The model variant name.}
#'   \item{clustering_}{The FlowSOM clustering ID used.}
#'   \item{prediction_name}{The name of the predicted class/label.}
#'   \item{prediction_value}{The predicted probability or score.}
#'   \item{threshold_proc_validationset}{The classification threshold used.}
#' }
#'
#' @export
models_grouped_apply <- function(
    df,
    outdir,
    result_grouping,
    grouping,
    dfcol_train_validation_other,
    dfcol_outcomes,
    counts_proportions = "proportions_per_x",
    bygroup = FALSE,
    ...) {
    # If bygroup is TRUE, only keep rows matching grouping
    if (bygroup) {
        grouping_noTVT <- grouping |>
            dplyr::select(-tidyr::all_of(dfcol_train_validation_other))

        # Join grouping metadata with feature data
        df_in_group <- dplyr::left_join(
            grouping_noTVT,
            df,
            by = colnames(grouping_noTVT)
        )
    } else {
        df_in_group <- df
    }

    # Retrieve trained models
    models <- result_grouping[["models"]][["final_models"]]


    # Apply all models across all outcomes and clusterings
    predictions_long <- sapply(names(models), simplify = FALSE, function(outcome_x) {
        models_x <- models[[outcome_x]]
        lapply(models_x, function(model_variant) {
            sapply(names(model_variant), simplify = FALSE, function(clustering_x) {
                single_model <- model_variant[[clustering_x]]
                clustering_df <- result_grouping[["clustering"]][[counts_proportions]][[clustering_x]]

                # Join phenotype cluster proportions with input data
                df_in_group_pheno_values <- dplyr::left_join(df_in_group, clustering_df, by = "File")

                # Apply the model to the input data
                preds <- single_model$predict_newdata(df_in_group_pheno_values)

                # Keep only prediction columns (drop metadata)
                preds_dt <- data.table::as.data.table(preds)
                tmp <- preds_dt[, -which(colnames(preds_dt) %in% c("row_ids", "truth", "response")), with = FALSE]

                # Merge predictions with input
                df_preds <- cbind(df_in_group_pheno_values, tmp)

                # Convert wide predictions to long format
                df_preds_long <- tidyr::pivot_longer(
                    df_preds,
                    cols = dplyr::all_of(colnames(tmp)),
                    names_to = "prediction_name",
                    values_to = "prediction_value"
                )

                # Add threshold used for classification
                df_preds_long[["threshold_proc_validationset"]] <-
                    single_model[["model"]][["threshold_proc_closest.topleft"]]

                return(df_preds_long)
            }) |>
                data.table::rbindlist(idcol = "clustering_")
        }) |>
            data.table::rbindlist(idcol = "model_")
    }) |>
        data.table::rbindlist(idcol = "outcome_")

    return(predictions_long)
}
