#' Compare Flow Cytometry Data Across Devices
#'
#' This function performs a comparative analysis of flow cytometry data across multiple devices.
#' It includes basic sample statistics, gating, density plots, FlowSOM clustering, and marker intensity comparisons.
#'
#' @param flowframes A named list of `flowFrame` objects containing flow cytometry data.
#' @param df A `data.table` containing metadata with at least the columns `"File"`, `"Device"`, and `"Sample"`.
#' @param ff_columns_relevant A character vector specifying the relevant markers for analysis.
#' @param transformlist A transformation function or a named list of functions for transforming marker intensities
#'        (default: `asinh(x / 1e3)`).
#' @param gatingsets A named list of gating sets for each dataset.
#' @param gatename_primary A character string specifying the primary gating population.
#' @param n_events_postgate An integer specifying the maximum number of events to keep post-gating.
#' @param marker_to_gate A named vector mapping marker names to their corresponding gates.
#' @param device_colors A named vector or function that provides colors for each device.
#'        If a function is provided, it should take the number of devices as input and return a vector of colors.
#' @param nClus An integer specifying the number of clusters for FlowSOM clustering (default: 5).
#' @param scale A logical indicating whether to scale the data in FlowSOM clustering (default: `FALSE`).
#' @param xdim An integer specifying the x-dimension of the FlowSOM grid (default: 3).
#' @param ydim An integer specifying the y-dimension of the FlowSOM grid (default: 3).
#' @param seed An integer specifying the random seed for FlowSOM clustering (default: `3711283`).
#' @param ... Additional parameters passed to `plot_flowsom()`.
#'
#' @return A named list of ggplot2 objects containing:
#'   \item{"Samples over time per device"}{A plot showing the number of samples collected over time per device.}
#'   \item{"Counts and percentages"}{Plots of gated cell counts and percentages per sample.}
#'   \item{"Positive population MFI"}{Plots showing the median fluorescence intensity (MFI) of positive gated populations.}
#'   \item{"Density plots"}{Density distributions of marker intensities across devices and samples.}
#'   \item{"Flowsom_PCA"}{Principal Component Analysis (PCA) plots of FlowSOM clustering results.}
#'   \item{"Flowsom_MA"}{MA plots comparing cluster proportions between devices.}
#'
#' @export

cycompare_outcomes <- function(
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
    clustering_n_subsampling = 1, # Number of subsamplings
    clustering_n_subsampled_cells = 1e4, # Number of cells to be subsampled (up or downsampling automatically) for FlowSOM
    clustering_subsampling_seed_first = 42, # Seed for the first subsampling
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
                max.depth = paradox::to_tune(2, 20), # minimum and maximum depth
                num.trees = paradox::to_tune(c(500, 1000, 1500, 2000)),
                importance = "impurity"
            )
        ),
        dv_class_positive = c("outcome_1" = "A", "outcome_2" = 5.1),
        loss_measure = mlr3::msr("classif.logloss")
    ),
    ...) {
    prepared <- cycompare_preparation(
        flowframes = flowframes,
        df = df,
        ff_columns_relevant = ff_columns_relevant,
        device_colors = device_colors,
        gatingsets = gatingsets,
        gatename_primary = gatename_primary,
        n_events_postgate = n_events_postgate
    )
    gated_ff <- prepared[["gated_ff"]]
    counts_joint <- prepared[["counts_joint"]]
    device_colors <- prepared[["device_colors"]]

    #### Transform the data
    if (all(is.null(transformList))) {
        transformList <- NULL
    } else if (is.function(transformList)) {
        transformList <- flowCore::transformList(
            from = ff_columns_relevant,
            tfun = transformList
        )
    } else {
        # I expect that this is then a proper transformList for flowCore
        if (!"transformList" %in% class(transformList)) {
            stop("transformList should be NULL, a function or a transformList object")
        }
    }

    if (!all(is.null(transformList))) {
        gated_transformed_ff <- lapply(gated_ff, function(ff_x) {
            flowCore::transform(ff_x, transformList)
        })
    } else {
        gated_transformed_ff <- gated_ff
    }
    gated_ff <- NULL

    possible_groupings <- df |>
        dplyr::select(dfcol_grouping_supersamples, dfcol_grouping_samples, dfcol_train_validation_other) |>
        dplyr::distinct()
    possible_groupings_training <- possible_groupings |>
        dplyr::filter(!!rlang::sym(dfcol_train_validation_other) == "train")



    # #### 1. Basic plots
    # ## 1.1 Outcomes per study
    # p1.1 <- plot_outcome_circles(
    #     df,
    #     dfcol_grouping_supersamples = dfcol_grouping_supersamples,
    #     dfcol_outcomes = dfcol_outcomes
    # )



    ### Clustering
    clusterings_ontrain <- do.call(
        fun_grouped,
        c(
            list(
                data = gated_transformed_ff,
                fun = flowsom_repeatsubsampling,
                df = df |> dplyr::filter(
                    # Clustering should only be trained on the TRAINING set samples
                    !!rlang::sym(dfcol_train_validation_other) == "train"
                ),
                dfcol_grouping_supersamples = dfcol_grouping_supersamples,
                dfcol_grouping_samples = dfcol_grouping_samples,
                dfcol_train_validation_other = dfcol_train_validation_other,
                outdir_base = file.path(outdir_base, "clustering"),
                verbose = TRUE,
                # FlowSOM parameters
                columns_clustering = ff_columns_relevant,
                n_subsampling = clustering_n_subsampling, # Number of subsamplings
                n_subsampled_cells = clustering_n_subsampled_cells, # Number of cells to be subsampled (up or downsampling automatically)
                subsampling_seed_first = clustering_subsampling_seed_first,
                # Do not transform the data within FlowSOM, it has been transformed already!
                transform = FALSE
            ),
            kwargs_flowsom
        )
    )

    ### Apply clustering to all samples
    applied_fs <- fun_grouped_apply(
        data = gated_transformed_ff,
        fun = flowsom_repeatsubsampling_apply,
        result_grouping = clusterings_ontrain,
        outdir_base = file.path(outdir_base, "clustering_applied"),
        verbose = FALSE,
        return_results = TRUE,
        remove_results_keywords = c("flowsom_newdata", "cells_clusters_from_train")

        # if n_metacluster is not provided, it will be taken from the FlowSOM result
        # n_metacluster = kwargs_flowsom[["nClus"]]
    )

    browser()

    # options(
    #     future.globals.maxSize = 10 * 1024^3
    # )
    pacman::p_load("mlr3learners")

    models_fs <- fun_grouped_apply(
        data = df,
        result_grouping = applied_fs,
        make_flowset = FALSE,
        fun = modelling,
        outdir_base = file.path(outdir_base, "models"),
        verbose = FALSE,
        return_results = TRUE,
        dfcol_train_validation_other = dfcol_train_validation_other,
        dfcol_outcomes = dfcol_outcomes,
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
        loss_measure = mlr3::msr("classif.logloss")
    )

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
                dfcol_outcomes = dfcol_outcomes
            ),
            kwargs_modelling
        )
    )
    return()
}

