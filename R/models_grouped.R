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
            predict_type = "prob", predict_sets = c("train", "test"),
            max.depth = paradox::to_tune(2, 20), # minimum and maximum depth
            num.trees = paradox::to_tune(c(500, 1000, 1500, 2000)),
            importance = "impurity"
        )
    ),
    dv_class_positive = c("outcome_1" = "A", "outcome_2" = 5.1),
    loss_measure = mlr3::msr("classif.logloss"),
    ...) {
    grouping_noTVT <- grouping |> dplyr::select(-tidyr::all_of(dfcol_train_validation_other))
    df_in_group <- dplyr::left_join(
        grouping_noTVT,
        df,
        by = colnames(grouping_noTVT)
    )
    tmp_df <- result_grouping[["proportions_per_x"]][["metaCluster"]]
    df_in_group_pheno_values <- dplyr::left_join(df_in_group, tmp_df, by = "File")
    modelx <- cytobench::wrapper_count_models(
        # modelx <- wrapper_count_models(
        df_list = list("metaCluster" = df_in_group_pheno_values),
        tvt_col = dfcol_train_validation_other,
        outdir = outdir,
        dvs_potential = dfcol_outcomes,
        dvs_multiclass = c(),
        ivs_regex = "[cC]luster_",
        hparam_n_evaluations = hparam_n_evaluations,
        seed = seed,
        learners_classification = learners_classification,
        dv_class_positive = dv_class_positive,
        measures = loss_measure,
        hpoptimized_final_trainsets = c("train"),
    )
    return(modelx)
}
