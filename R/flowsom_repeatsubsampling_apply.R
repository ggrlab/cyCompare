#' Apply FlowSOM prediction on repeated subsamplings
#'
#' This function performs repeated subsampling of a list of flowFrames, applies a trained FlowSOM model
#' to the subsampled data, and returns predicted cluster results with normalized cluster proportions.
#'
#' @param ff_list A list of `flowFrame` objects to be subsampled and clustered.
#' @param outdir A character string indicating the output directory where results will be saved.
#'   If `outdir` exists, the results will be written using `flowsom_wrapper_saving()`.
#' @param flowsom_result A trained FlowSOM object (typically returned by `flowsom_repeatsubsampling()` which called `cytobench::flowSOM_optimal()`).
#' @param n_metacluster Integer. Number of metaclusters to use when predicting clusters.
#'   If NULL, the number of metaclusters in `flowsom_result` will be used.
#' @param n_subsampling Integer. Number of repeated subsampling iterations to perform. Default is 1.
#' @param n_subsampled_cells Integer. Number of cells per subsample. If the original sample is smaller,
#'   upsampling is performed; otherwise, downsampling. If you want to use the original sample size,
#'  set `n_subsampled_cells` to `Inf`. See `cytobench::subsample_ff()` for details.
#' @param subsampling_seed_first Integer. Base seed for reproducible subsampling. Each subsample iteration
#'   uses `subsampling_seed_first + i - 1`.
#' @param remove_results_keywords Character vector. Names of result slots in `flowsom_result` to be removed from output.
#'   This is useful to avoid returning large objects that are not needed for further analysis.
#' @param ... Not used here, but potentially relevant if the calling function is used in a
#' pipeline with additional parameters.
#'
#' @return A nested list of FlowSOM prediction results. For each data slot, proportions of clusters per sample
#'   are returned (after normalization). Slots defined in `remove_results_keywords` are removed. Each list entry
#'   contains subsampling-specific results.
#'
#' @export
flowsom_repeatsubsampling_apply <- function(
    ff_list,
    outdir,
    flowsom_result,
    n_metacluster = NULL,
    n_subsampling = 1, # Number of subsamplings
    n_subsampled_cells = 10000, # Number of cells to be subsampled (up or downsampling automatically)
    subsampling_seed_first = 427764,
    remove_results_keywords = c("flowsom_newdata", "cells_clusters_from_train"),
    ...) {
    # Subsample the data
    subsample_multiple <- unlist(
        lapply(1:n_subsampling, function(i) {
            fs_train_subsampled <- flowCore::fsApply(ff_list, function(x) {
                cytobench::subsample_ff(x, n_cells = n_subsampled_cells, seed = subsampling_seed_first + i - 1)
            })
            flowCore::sampleNames(fs_train_subsampled) <- paste0(
                flowCore::sampleNames(fs_train_subsampled),
                "_subsampled",
                i
            )
            return(flowCore::flowSet_to_list(fs_train_subsampled))
        })
    )
    fs_multiple <- flowCore::flowSet(subsample_multiple)
    fs_full_predictedFS <- cytobench::flowSOM_predict(
        flowsom_result = flowsom_result,
        flowset = fs_multiple,
        n_metacluster = min(
            n_metacluster,
            # Maximum number of metaclusters from FlowSOM itself:
            flowsom_result[["fs_res_train"]][["map"]][["nMetaclusters"]]
        )
    )
    fs_full_predictedFS[["proportions_per_x"]] <- lapply(fs_full_predictedFS[["ncells_per_x"]], function(x) {
        cluster_columns <- colnames(x)[grepl("[cC]luster", colnames(x))]
        ncells_per_row <- rowSums(x[, cluster_columns])
        ncells_per_row_mat <- matrix(
            ncells_per_row,
            ncol = length(cluster_columns), nrow = nrow(x),
            byrow = FALSE
        )
        # column1, column2, column3
        # ncells_row1, ncells_row1, ncells_row1
        # ncells_row2, ncells_row2, ncells_row2
        # ...
        x[, cluster_columns] <- x[, cluster_columns] / ncells_per_row_mat
        if (!all(rowSums(x[, cluster_columns]) - 1 < .Machine$double.eps)) {
            # I know this check is slow, but these matrices are going to be small anyways
            # Especially irrelevant regarding the previous clustering application.
            stop("Proportions do not sum to 1 per sample.")
        }
        x
    })
    for (kw_x in remove_results_keywords) {
        fs_full_predictedFS[[kw_x]] <- NULL
    }
    if (exists(outdir)) {
        flowsom_wrapper_saving(
            fs_applied = fs_full_predictedFS,
            dir_out = dirname(outdir),
            modelname = basename(outdir)
        )
    }
    fs_full_predictedFS <- lapply(fs_full_predictedFS, function(x) {
        lapply(x, function(y) {
            if ("sample" %in% colnames(y)) {
                y <- dplyr::mutate(
                    y,
                    File = sub("_subsampled[0-9]+$", "", sample),
                    .before = "sample"
                )
            }
            return(y)
        })
    })

    return(fs_full_predictedFS)
}
