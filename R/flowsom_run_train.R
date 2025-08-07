#' Run FlowSOM Clustering and Predict Cluster Assignments
#'
#' Transforms the data, selects training samples, performs FlowSOM clustering,
#' and predicts cluster assignments for all samples.
#'
#' @param ff_gated A list of flowFrame objects.
#' @param df Metadata table with a `File` column.
#' @param transformlist Transformation function(s).
#' @param nClus Number of clusters.
#' @param scale Logical; if TRUE, scale data.
#' @param xdim, ydim Dimensions for FlowSOM grid.
#' @param seed Random seed.
#' @param dfcol_train_validation_other Optional column name for filtering training data (e.g., "TrainTest").
#'
#' @return Output of `cytobench::flowSOM_predict()`.
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
