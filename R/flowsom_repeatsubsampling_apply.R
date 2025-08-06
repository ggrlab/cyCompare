#' Apply FlowSOM Prediction After Repeated Subsampling
#'
#' Performs repeated subsampling of  a list of flowFrames and applies a trained FlowSOM model
#' to each subsample. Cluster proportions are normalized per sample. Results can be saved to disk
#' or returned in memory.
#'
#' @param ff_list A list of `flowFrame` objects to be subsampled and clustered.
#' @param outdir A character string indicating the output directory where results will be saved.
#'   If the path exists, results will be written using `flowsom_wrapper_saving()`.
#' @param flowsom_result
#' A trained FlowSOM result (usually from `flowsom_repeatsubsampling()` or
#' `cytobench::flowSOM_optimal()`).
#' @param n_metacluster
#' Integer. Number of metaclusters to use when predicting.
#'   If NULL, the number of metaclusters in `flowsom_result` will be used.
#' @param n_subsampling Integer. Number of subsampling iterations (default: 1).
#' @param n_subsampled_cells Integer. Number of cells per subsample. If the original sample is smaller,
#'  upsampling is performed; otherwise, downsampling. See `cytobench::subsample_ff()` for details.
#'  Use `Inf` to avoid subsampling.
#' @param subsampling_seed_first
#' Integer. Base seed used for deterministic subsampling.
#' Each subsample iteration uses `subsampling_seed_first + i - 1`.
#' @param remove_results_keywords Character vector of component names to strip from the result.
#' Names of result slots in `flowsom_result` to be removed from output.
#' This is useful to avoid returning large objects that are not needed for further analysis.
#' @param ...
#' Not used here, but potentially relevant if the calling function is used in a
#' pipeline with additional parameters.
#'
#' @return A nested list containing FlowSOM prediction results for each subsample, with cluster proportions normalized per sample.
#' Large intermediate objects can be removed based on `remove_results_keywords`. Each list entry
#'   contains subsampling-specific results.
#'
#' @export
flowsom_repeatsubsampling_apply <- function(
    ff_list,
    outdir,
    flowsom_result,
    n_metacluster = NULL,
    n_subsampling = 1,
    n_subsampled_cells = 10000,
    subsampling_seed_first = 427764,
    remove_results_keywords = c("flowsom_newdata", "cells_clusters_from_train"),
    ...) {
    # --- 1. Subsample the input flowFrames multiple times ---
    multiple_subsamples <- subsample_multiple(
        ff_list = ff_list,
        n_subsampling = n_subsampling,
        n_subsampled_cells = n_subsampled_cells,
        subsampling_seed_first = subsampling_seed_first
    )

    # Combine all subsampled flowFrames into a single flowSet
    fs_multiple <- flowCore::flowSet(multiple_subsamples)

    # --- 2. Apply the trained FlowSOM model to the new subsampled data ---
    fs_full_predictedFS <- cytobench::flowSOM_predict(
        flowsom_result = flowsom_result,
        flowset = fs_multiple,
        n_metacluster = min(
            n_metacluster %||% Inf,
            # Ensure we do not exceed the number of metaclusters used in training
            flowsom_result[["fs_res_train"]][["map"]][["nMetaclusters"]]
        )
    )

    # --- 3. Normalize cluster counts to proportions for each sample ---
    fs_full_predictedFS[["proportions_per_x"]] <- lapply(
        fs_full_predictedFS[["ncells_per_x"]],
        function(x) {
            cluster_cols <- grep("[cC]luster", colnames(x), value = TRUE)
            ncells_per_row <- rowSums(x[, cluster_cols])
            ncells_per_row_mat <- matrix(
                ncells_per_row,
                nrow = nrow(x),
                ncol = length(cluster_cols),
                byrow = FALSE
            )
            x[, cluster_cols] <- x[, cluster_cols] / ncells_per_row_mat
            if (!all(abs(rowSums(x[, cluster_cols]) - 1) < .Machine$double.eps)) {
                stop("Proportions do not sum to 1 per sample.")
            }
            x
        }
    )

    # --- 4. Remove large or unnecessary intermediate slots ---
    for (kw in remove_results_keywords) {
        fs_full_predictedFS[[kw]] <- NULL
    }

    # --- 5. Save results to disk (optional) ---
    if (!is.null(outdir)) {
        flowsom_wrapper_saving(
            fs_applied = fs_full_predictedFS,
            dir_out = dirname(outdir),
            modelname = basename(outdir)
        )
    }

    # --- 6. Clean up sample names in each data.frame ---
    fs_full_predictedFS <- lapply(fs_full_predictedFS, function(lst) {
        lapply(lst, function(tbl) {
            if ("sample" %in% colnames(tbl)) {
                tbl <- dplyr::mutate(tbl, File = sub("_subsampled[0-9]+$", "", sample), .before = "sample")
            }
            return(tbl)
        })
    })

    return(fs_full_predictedFS)
}
