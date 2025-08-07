#' Train Grouped Classification Models Using Clustered Flow Cytometry Data
#'
#' This function trains classification models on clustered flow cytometry data grouped by metadata.
#' The grouping information is used to join phenotypic proportions (from FlowSOM results) to
#' the feature data, which is then passed to `cytobench::wrapper_count_models` for model training.
#'
#' @inheritParams cycompare_outcomes_analyse
#' @param outdir A character string specifying the output directory for model results.
#' @param result_grouping A list containing precomputed grouping results (e.g., from FlowSOM clustering),
#' typically the `applied` result from `fun_grouped_apply()`. Must include `"proportions_per_x"` slot.
#' @param grouping A data frame that defines grouping variables for merging with `df`.
#' @param hparam_n_evaluations Integer; number of hyperparameter evaluations during model tuning. Default is 3.
#' @param seed Integer; random seed for reproducibility. Default is 42.
#' @param learners_classification A named list of [mlr3::Learner] objects with hyperparameter ranges defined.
#' Defaults to a tuned random forest (`classif.ranger`).
#' @param dv_class_positive Named vector specifying the positive class for each outcome.
#' This ensures consistent performance evaluation for binary classification tasks.
#' @param loss_measure An `mlr3` measure object to evaluate model performance. Default is classification log-loss.
#' @param ... Additional arguments passed to `cytobench::wrapper_count_models`.
#'
#' @return A list of trained models and associated evaluation metadata, as returned by `cytobench::wrapper_count_models`.
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
            max.depth = paradox::to_tune(2, 20),
            num.trees = paradox::to_tune(c(500, 1000, 1500, 2000)),
            importance = "impurity"
        )
    ),
    dv_class_positive = c("outcome_1" = "A", "outcome_2" = 5.1),
    loss_measure = mlr3::msr("classif.logloss"),
    ...) {
    # Remove the TVT (train/validation/test) assignment columns to form the unique grouping key
    grouping_noTVT <- grouping |> dplyr::select(-tidyr::all_of(dfcol_train_validation_other))

    # Merge the grouping metadata with the feature data
    df_in_group <- dplyr::left_join(
        grouping_noTVT,
        df,
        by = colnames(grouping_noTVT)
    )

    # Todo: Why do I exclusively use "metaCluster" here?
    browser()
    # Extract phenotype proportions (e.g., cluster frequencies) from FlowSOM result
    tmp_df <- result_grouping[["proportions_per_x"]][["metaCluster"]]

    # Join feature data with phenotype proportions based on sample identifier
    df_in_group_pheno_values <- dplyr::left_join(df_in_group, tmp_df, by = "File")

    # Train classification models using outcome-specific learners and tuning settings
    # Todo: Why do I exclusively use "metaCluster" here?
    browser()
    modelx <- cytobench::wrapper_count_models(
        df_list = list("metaCluster" = df_in_group_pheno_values),
        tvt_col = dfcol_train_validation_other,
        outdir = outdir,
        dvs_potential = dfcol_outcomes,
        dvs_multiclass = c(), # Only binary outcomes supported for now
        ivs_regex = "[cC]luster_", # Feature columns must match this pattern
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
