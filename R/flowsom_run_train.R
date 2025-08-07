#' Train and Predict FlowSOM Clustering Based on Training Samples
#'
#' This function applies marker transformation, runs FlowSOM clustering on a subset of training samples,
#' and predicts cluster assignments for all samples. Clustering is performed using `flowsom_repeatsubsampling()`.
#'
#' @param ff_gated A list of `flowFrame` objects or a `flowSet`, representing gated cytometry data.
#' @param df
#' A `data.table` containing metadata with at least a `File` column. Only used
#' if `dfcol_train_validation_other` is specified, then `dfcol_train_validation_other` column
#' is used to select training samples for training the clustering.
#' @param transformlist
#' Either a single transformation function or a named list of functions per marker.
#' @param seed Integer seed for reproducibility of subsampling.
#' @param dfcol_train_validation_other
#' Optional column in `df` used to select training samples (e.g., `"TrainTest"`).
#' @param relevant_columns
#' Character vector specifying markers to use for clustering. Defaults to all columns.
#' @param ... Additional arguments passed to `flowsom_repeatsubsampling()`.
#'
#' @return A list containing cluster prediction results, as returned by `cytobench::flowSOM_predict()`.
#'
#' @export
flowsom_run_train <- function(ff_gated,
                              df,
                              transformlist,
                              seed = 1,
                              dfcol_train_validation_other = NULL,
                              relevant_columns,
                              ...) {
    File <- NULL
    if (!inherits(ff_gated, "flowSet")) {
        ff_gated <- flowCore::flowSet(ff_gated)
    }
    if (is.null(relevant_columns)) {
        relevant_columns <- flowCore::colnames(ff_gated[[1]])
    }

    # Build transformlist
    fc_transformlist <- flowCore::transformList(
        flowCore::colnames(ff_gated[[1]]),
        transformlist_named(transformlist, flowCore::colnames(ff_gated[[1]]))
    )

    # Apply transformation
    fs_transformed <- flowCore::transform(ff_gated, fc_transformlist)

    # Filter training samples if applicable
    if (!is.null(dfcol_train_validation_other)) {
        train_ids <- df[df[[dfcol_train_validation_other]] == "train", File]
        fs_train <- fs_transformed[train_ids]
    } else {
        fs_train <- fs_transformed
    }
    # Run FlowSOM clustering only on training samples - if specified, otherwise on all
    flowsom_model <- flowsom_repeatsubsampling(
        fs_train,
        outdir = NULL,
        columns_clustering = relevant_columns,
        subsampling_seed_first = seed,
        ...
    )
    # Predict clusters for all samples
    cytobench::flowSOM_predict(flowsom_model[["fs_res_train"]], fs_transformed)
}
