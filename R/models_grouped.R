#' Train Grouped Models Using a Classification Pipeline
#'
#' This function applies classification models to flow cytometry data grouped by metadata.
#' It handles preprocessing, joins with grouping information, and fits models using `cytobench::wrapper_count_models`.
#'
#' @param df A data frame containing feature columns used for classification.
#' @param outdir A character string specifying the output directory to save model results.
#' @param result_grouping A list object containing grouping results, including phenotype proportions.
#' @param grouping A data frame that defines grouping variables for merging with `df`.
#' @param dfcol_train_validation_other
#'  Character vector of column names indicating data split (e.g., train/validation/test).
#' @param dfcol_outcomes Character vector of outcome column names to be modeled.
#' @param hparam_n_evaluations Integer; number of hyperparameter evaluations during model tuning. Default is 3.
#' @param seed Integer; seed for reproducibility. Default is 42.
#' @param learners_classification
#'  A list of mlr3 learners with tuning parameters set. Default includes a tuned `classif.ranger`.
#' @param dv_class_positive A named vector specifying the positive class for each outcome variable.
#' @param loss_measure An `mlr3` measure object for model evaluation. Default is classification log-loss.
#' @param ... Additional arguments passed to `cytobench::wrapper_count_models`.
#'
#' @return A list of trained models and evaluation results.
#' @export
models_grouped <- function(
    df,
    outdir,
    result_grouping,
    grouping,
    dfcol_train_validation_other,
    dfcol_outcomes,
    hparam_n_evaluations = 3,
    seed = 42,
    learners_classification = list(
        mlr3::lrn(
            "classif.ranger",
            predict_type = "prob",
            predict_sets = c("train", "test"),
            max.depth = paradox::to_tune(2, 20), # range for max depth of trees
            num.trees = paradox::to_tune(c(500, 1000, 1500, 2000)), # trees to evaluate
            importance = "impurity"
        )
    ),
    dv_class_positive = c("outcome_1" = "A", "outcome_2" = 5.1),
    loss_measure = mlr3::msr("classif.logloss"),
    ...) {
    # Remove the train/validation/test columns to form grouping keys only
    grouping_noTVT <- grouping |> dplyr::select(-tidyr::all_of(dfcol_train_validation_other))

    # Merge the grouping metadata with the feature data
    df_in_group <- dplyr::left_join(
        grouping_noTVT,
        df,
        by = colnames(grouping_noTVT)
    )

    # Get phenotype proportions per file from result_grouping
    tmp_df <- result_grouping[["proportions_per_x"]][["metaCluster"]]

    # Join the phenotype data to the merged data
    df_in_group_pheno_values <- dplyr::left_join(df_in_group, tmp_df, by = "File")

    # Fit models using wrapper_count_models from cytobench
    modelx <- cytobench::wrapper_count_models(
        df_list = list("metaCluster" = df_in_group_pheno_values),
        tvt_col = dfcol_train_validation_other,
        outdir = outdir,
        dvs_potential = dfcol_outcomes,
        dvs_multiclass = c(), # currently only binary classification is supported
        ivs_regex = "[cC]luster_", # regex to identify input variables
        hparam_n_evaluations = hparam_n_evaluations,
        seed = seed,
        learners_classification = learners_classification,
        dv_class_positive = dv_class_positive,
        measures = loss_measure,
        hpoptimized_final_trainsets = c("train"),
        ...
    )

    return(modelx)
}
