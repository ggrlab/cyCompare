#' Train Grouped Models Using a Classification Pipeline
#'
#' This function applies classification models to flow cytometry data grouped by metadata.
#' It handles preprocessing, joins with grouping information, and fits models using `cytobench::wrapper_count_models`.
#'
#' @param df A data frame containing feature columns used for classification.
#' @param outdir A character string specifying the output directory to save model results. Not used here.
#' @param result_grouping A list object containing grouping results, including phenotype proportions.
#' @param grouping A data frame that defines grouping variables for merging with `df`.
#' @param result_grouping_models A list of models to be applied to the data.
#' @param bygroup Logical;
#'   If FALSE, the function will apply each model to the entire dataset (result_grouping[["clustering"]][["proportions_per_x"]][["metaCluster"]]).
#'   If TRUE, the `FALSE` dataset is limited to the rows mathing the `grouping` variable.
#' @param dfcol_train_validation_other Character vector of column names indicating data split (e.g., train/validation/test).
#' @param dfcol_outcomes Character vector of outcome column names to be modeled.
#' @param ... Not used.
#'
#' @return A list of trained models and evaluation results.
#' @export
models_grouped_apply <- function(
    df,
    outdir,
    result_grouping,
    grouping,
    result_grouping_models,
    dfcol_train_validation_other,
    dfcol_outcomes,
    bygroup = FALSE,
    ...) {
    if (bygroup) {
        # Remove the train/validation/test columns to form grouping keys only
        grouping_noTVT <- grouping |> dplyr::select(-tidyr::all_of(dfcol_train_validation_other))

        # Merge the grouping metadata with the feature data
        df_in_group <- dplyr::left_join(
            grouping_noTVT,
            df,
            by = colnames(grouping_noTVT)
        )
    } else {
        df_in_group <- df
    }

    # Get phenotype proportions per file from result_grouping
    clustering_df <- result_grouping[["clustering"]][["proportions_per_x"]][["metaCluster"]]
    models <- result_grouping[["models"]][["final_models"]]

    # Join the phenotype data to the merged data
    df_in_group_pheno_values <- dplyr::left_join(df_in_group, clustering_df, by = "File")

    # Fit models using wrapper_count_models from cytobench
    predictions_long <- sapply(names(models), simplify = FALSE, function(outcome_x) {
        models_x <- models[[outcome_x]]
        lapply(models_x, function(model_variant) {
            sapply(names(model_variant), simplify = FALSE, function(clustering_x) {
                single_model <- model_variant[[clustering_x]]
                preds <- single_model$predict_newdata(df_in_group_pheno_values)
                tmp <- data.table::as.data.table(preds)[, -c(1:3)]

                df_preds <- cbind(
                    df_in_group_pheno_values,
                    tmp
                )
                df_preds_long <- tidyr::pivot_longer(
                    df_preds,
                    cols = dplyr::all_of(colnames(tmp)),
                    names_to = "prediction_name",
                    values_to = "prediction_value"
                )
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
