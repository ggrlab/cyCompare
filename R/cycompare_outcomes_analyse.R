#' Core Analysis Pipeline for Flow Cytometry Outcomes Across Devices
#'
#' Performs a full pipeline for outcome prediction across multiple flow cytometry devices.
#' This includes data preprocessing (gating, transformation), FlowSOM clustering, and outcome modeling.
#' The function returns results that can be directly used for downstream statistical analysis or plotting.
#'
#' @param flowframes
#' Named list of `flowFrame` objects with raw flow cytometry data. Names must match entries in `df$File`.
#' @param df
#' A `data.table` containing sample metadata. Must include columns specified by `dfcol_*`
#' and a unique `File` column linking to `names(flowframes)`..
#' @param ff_columns_relevant
#' Character vector of marker/feature columns to use for clustering and modeling.
#' @param dfcol_grouping_supersamples
#' Column(s) in `df` used to define supersample groupings (e.g., `"Study"`).
#' A supersample collects samples of same origin: A study collects samples from
#' multiple donors/patients, a patient might be in multiple times.
#' @param dfcol_grouping_samples
#' Column in `df` that defines sample-level grouping (e.g., `"Device"`).
#'  Each final biologically different sample could be measured with a different technology,
#' cytometer, or lab but SHOULD have the same values.
#' @param dfcol_train_validation_other
#' Column indicating sample usage role: `"train"`, `"test"`, etc.
#' @param dfcol_outcomes
#'  Character vector of outcome column names to model (e.g., `c("outcome_1", "outcome_2")`).
#' @param outcome_models
#' Named list of modeling functions. Defaults to `cv.glmnet`.
#' @param outdir_base
#' Base directory for saving intermediate results (default: `tempdir()`).
#' @param transformlist
#' Either a function (applied as `flowCore::transformList`) or a `transformlist` object.
#' Default: `asinh(x / 1e3)`.
#' @param gatingsets
#' Named list of gating sets used to extract primary populations. One gating set per sample.
#' @param gatename_primary
#' String specifying the gate used for downstream data (e.g., `"CD3+"`).
#' @param n_events_postgate
#' Maximum number of events to retain per sample after gating (default: 10,000).
#' @param marker_to_gate
#' Named vector mapping markers to their respective gate names.
#' @param device_colors
#' Either a named vector or a function returning colors for each device.
#' @param clustering_n_subsampling
#' Number of FlowSOM subsampling repeats (default: 1).
#' @param clustering_n_subsampled_cells
#' Number of cells per sample to use per subsample (default: 10,000).
#' @param clustering_subsampling_seed_first
#' Seed for the first clustering repeat (default: 42).
#' @param kwargs_flowsom
#' Named list of FlowSOM parameters (e.g., `nClus`, `xdim`, `ydim`, `seed`, etc.).
#' @param kwargs_modelling
#' Named list of modeling parameters (e.g., learners, tuning, loss function).
#' @param dv_class_positive
#'  Named vector specifying the positive class for each outcome (e.g., `c("outcome_1" = "A", "outcome_2" = 5.1)`).
#'
#' @return A named list with the following structure:
#' \describe{
#'   \item{prepared_data}{Transformed flowFrames, cell counts, and device colors.}
#'   \item{clustering}{FlowSOM models trained on training data and applied to all samples.}
#'   \item{models}{
#'      Trained models according to `dfcol_train_validation_other` splits
#'      and their predictions across all samples.
#'  }
#' }
#'
#' @export
cycompare_outcomes_analyse <- function(
    flowframes,
    df,
    ff_columns_relevant,
    dfcol_grouping_supersamples = c("Study"),
    dfcol_grouping_samples = "Device",
    dfcol_train_validation_other = "tvt",
    dfcol_outcomes = c("outcome_1", "outcome_2"),
    outcome_models = list("glmnet" = glmnet::cv.glmnet),
    outdir_base = tempdir(),
    transformlist = function(x) asinh(x / 1e3),
    gatingsets,
    gatename_primary,
    n_events_postgate = 10e3,
    marker_to_gate,
    device_colors = function(n) {
        RColorBrewer::brewer.pal(n, "Dark2")
    },
    clustering_n_subsampling = 1,
    clustering_n_subsampled_cells = 1e4,
    clustering_subsampling_seed_first = 42,
    kwargs_flowsom = list(
        nClus = 5,
        scale = FALSE,
        xdim = 3,
        ydim = 3,
        seed = 3711283
    ),
    kwargs_modelling = list(
        hparam_n_evaluations = 3,
        seed = 42,
        learners_classification = list(
            mlr3::lrn(
                "classif.ranger",
                predict_type = "prob", predict_sets = c("train", "test"),
                max.depth = paradox::to_tune(2, 20),
                num.trees = paradox::to_tune(c(500, 1000, 1500, 2000)),
                importance = "impurity"
            )
        ),
        loss_measure = mlr3::msr("classif.logloss")
    ),
    dv_class_positive = c("outcome_1" = "A", "outcome_2" = 5.1)) {
    ### 1. Gating and preparation
    prepared <- cycompare_preparation(
        flowframes = flowframes,
        df = df,
        ff_columns_relevant = ff_columns_relevant,
        device_colors = device_colors,
        gatingsets = gatingsets,
        gatename_primary = gatename_primary,
        n_events_postgate = n_events_postgate,
        dfcol_grouping_samples = dfcol_grouping_samples
    )
    gated_ff <- prepared[["gated_ff"]]
    counts_joint <- prepared[["counts_joint"]]
    device_colors <- prepared[["device_colors"]]

    ### 2. Transformation (e.g. asinh)
    if (all(is.null(transformlist))) {
        transformlist <- NULL
    } else if (is.function(transformlist)) {
        transformlist <- flowCore::transformList(
            from = ff_columns_relevant,
            tfun = transformlist
        )
    } else if (!"transformlist" %in% class(transformlist)) {
        stop("transformlist should be NULL, a function, or a transformlist object")
    }

    if (!is.null(transformlist)) {
        gated_transformed_ff <- lapply(gated_ff, function(ff_x) {
            flowCore::transform(ff_x, transformlist)
        })
    } else {
        gated_transformed_ff <- gated_ff
    }
    gated_ff <- NULL # Free memory

    # nolint start
    ### 3. Train FlowSOM clustering models on training samples
    ## fun_grouped applies "fun" to each unique combination of grouping columns in df
    # E.g. given that we restrict to "train" samples, with
    #   dfcol_grouping_supersamples = "Study"
    #   dfcol_grouping_samples = "Device"
    # The resulting two groups will be:
    # # A tibble: 2 x 3
    #   Study Device   train_validation_test
    #   <chr> <chr>    <chr>
    # 1 BAD   aurora   train
    # 2 BAD   fortessa train
    #
    # Apart from this "groups" return element, the "result" element is a list corresponding
    # to the groups in order of the [["groups"]] tibble.
    # nolint end
    clusterings_ontrain <- do.call(
        fun_grouped,
        c(
            list(
                data = gated_transformed_ff,
                fun = flowsom_repeatsubsampling,
                df = df |> dplyr::filter(
                    !!rlang::sym(dfcol_train_validation_other) == "train"
                ),
                dfcol_grouping_supersamples = dfcol_grouping_supersamples,
                dfcol_grouping_samples = dfcol_grouping_samples,
                dfcol_train_validation_other = dfcol_train_validation_other,
                outdir_base = file.path(outdir_base, "clustering"),
                verbose = TRUE,
                columns_clustering = ff_columns_relevant,
                n_subsampling = clustering_n_subsampling,
                n_subsampled_cells = clustering_n_subsampled_cells,
                subsampling_seed_first = clustering_subsampling_seed_first,
                transform = FALSE # Already transformed
            ),
            kwargs_flowsom
        )
    )

    ### 4. Apply clustering to all samples
    # The results given in "results_grouping" are applied to _all_ the `data` in fun_grouped_apply.
    # Effectively, all trained FlowSOM models are going to be applied to all samples in gated_transformed_ff.
    applied_fs <- fun_grouped_apply(
        data = gated_transformed_ff,
        fun = flowsom_repeatsubsampling_apply,
        result_grouping = clusterings_ontrain,
        outdir_base = file.path(outdir_base, "clustering_applied"),
        verbose = FALSE,
        return_results = TRUE,
        remove_results_keywords = c("flowsom_newdata", "cells_clusters_from_train")
    )

    # nolint start
    ### 5. Train outcome models on clustered training data
    # The results given in "results_grouping" are applied to _all_ the `data` in fun_grouped_apply.
    # models_grouped()
    #   - restricts df to the matching groups.
    #   - merges the df with the clustering results by group
    #   - trains models using the with cytobench::wrapper_count_models. This also ensures proper training/validation/test USAGE, not SPLITTING!
    # nolint end
    model_fs <- do.call(
        fun_grouped_apply,
        c(
            list(
                data = df,
                result_grouping = applied_fs,
                make_flowset = FALSE,
                fun = models_grouped,
                outdir_base = file.path(outdir_base, "models"),
                verbose = FALSE,
                return_results = TRUE,
                dfcol_train_validation_other = dfcol_train_validation_other,
                dfcol_outcomes = dfcol_outcomes,
                dv_class_positive = dv_class_positive
            ),
            kwargs_modelling
        )
    )

    ### 6. Apply models to all data
    # The results given in "results_grouping" are applied to _all_ the `data` in fun_grouped_apply.
    # If result_grouping is a list of elements, `fun` receives the grouping-specific data
    # from EACH element of the list.
    # Because "bygroup" is set to FALSE, the function will be applied to the entire `df`.
    model_fs_applied <- do.call(
        fun_grouped_apply,
        c(
            list(
                data = df,
                result_grouping = list(
                    "clustering" = applied_fs,
                    "models" = model_fs
                ),
                make_flowset = FALSE,
                fun = models_grouped_apply,
                outdir_base = file.path(outdir_base, "models_applied"),
                verbose = FALSE,
                return_results = TRUE,
                dfcol_train_validation_other = dfcol_train_validation_other,
                dfcol_outcomes = dfcol_outcomes,
                bygroup = FALSE
            ),
            kwargs_modelling
        )
    )

    ### 7. Return structured results for downstream use
    return(
        list(
            "prepared_data" = list(
                "gated_transformed_ff" = gated_transformed_ff,
                "counts_joint" = counts_joint,
                "device_colors" = device_colors
            ),
            "clustering" = list(
                "trained" = clusterings_ontrain,
                "applied" = applied_fs
            ),
            "models" = list(
                "trained" = model_fs,
                "applied" = model_fs_applied
            )
        )
    )
}
