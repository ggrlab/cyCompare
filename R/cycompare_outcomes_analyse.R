#' Core Analysis of Flow Cytometry Data Across Devices with outcomes
#'
#' This function performs preprocessing, transformation, FlowSOM clustering, and predictive modeling
#' on flow cytometry data to compare devices and generate outcome-based insights.
#' It returns intermediate structured data ready for downstream plotting or reporting.
#'
#' @param flowframes A named list of `flowFrame` objects containing flow cytometry data.
#' @param df A `data.table` with metadata. Must contain columns specified in `dfcol_*` parameters and `File` as unique column linking to `names(flowframes)`.
#' @param ff_columns_relevant Character vector of markers to use for analysis and clustering.
#' @param dfcol_grouping_supersamples Metadata column used to define supersample groups (e.g., "Study").
#' @param dfcol_grouping_samples Metadata column used to group samples (e.g., "Device").
#' @param dfcol_train_validation_other Column defining sample roles ("train", "test", etc.).
#' @param dfcol_outcomes Character vector of column names containing outcome variables.
#' @param outcome_models Named list of modeling functions, defaulting to `cv.glmnet`.
#' @param outdir_base Base directory for saving intermediate results (default: `tempdir()`).
#' @param transformList A transformation function or a `flowCore::transformList` object (default: `asinh(x / 1e3)`).
#' @param gatingsets A named list of gating sets, one per dataset.
#' @param gatename_primary Character string indicating the population used for downstream analysis.
#' @param n_events_postgate Maximum number of events per sample to retain after gating.
#' @param marker_to_gate Named vector mapping marker names to gate names.
#' @param device_colors Either a named vector or a function that generates colors for each device.
#' @param clustering_n_subsampling Number of repeated subsamplings for FlowSOM (default: 1).
#' @param clustering_n_subsampled_cells Number of cells per sample for FlowSOM (default: 10,000).
#' @param clustering_subsampling_seed_first Seed for the first subsampling (default: 42).
#' @param kwargs_flowsom Named list of additional FlowSOM parameters (e.g. `nClus`, `scale`, `xdim`, `ydim`, `seed`).
#' @param kwargs_modelling Named list of model training parameters (e.g., learners, tuning, seed).
#'
#' @return A named list with three main components:
#' \describe{
#'   \item{prepared_data}{List with transformed flowFrames, counts, and device colors.}
#'   \item{clustering}{List with FlowSOM models trained on training data and applied to all samples.}
#'   \item{models}{List with outcome models trained and applied to clustering results.}
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
    transformList = function(x) asinh(x / 1e3),
    gatingsets,
    gatename_primary,
    n_events_postgate = 10e3,
    marker_to_gate,
    device_colors = function(n) {
        RColorBrewer::brewer.pal(n, "Dark2")
    },
    # FlowSOM parameters
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
    ### 1. Gating and basic preparation
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

    ### 2.Transformation of flow data
    if (all(is.null(transformList))) {
        transformList <- NULL
    } else if (is.function(transformList)) {
        transformList <- flowCore::transformList(
            from = ff_columns_relevant,
            tfun = transformList
        )
    } else if (!"transformList" %in% class(transformList)) {
        stop("transformList should be NULL, a function, or a transformList object")
    }

    if (!is.null(transformList)) {
        gated_transformed_ff <- lapply(gated_ff, function(ff_x) {
            flowCore::transform(ff_x, transformList)
        })
    } else {
        gated_transformed_ff <- gated_ff
    }
    gated_ff <- NULL # Clean up memory

    # nolint start
    ### 3.Clustering using FlowSOM on training samples only
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

    ### 4.Apply trained clustering model to all samples
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
    ### 5.Train outcome models on clustered training data
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

    ### 6. Apply each trained models to all samples
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

    ### 7.Return structured results for downstream use
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
