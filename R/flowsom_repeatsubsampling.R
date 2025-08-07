#' Train FlowSOM Using Repeated Subsampling
#'
#' This function performs FlowSOM clustering after repeated subsampling of a flow cytometry training set.
#' Each subsample is taken independently using a reproducible seed and concatenated into one training set
#' for FlowSOM model fitting.
#'
#' @param ff_list A list or `flowSet` of flowFrames containing the training data.
#' @param outdir Optional output directory to save FlowSOM results. If `NULL`, results are returned only in memory.
#' @param columns_clustering Character vector of marker channels to use for clustering.
#' @param n_subsampling Integer; number of repeated subsampling iterations. Default: 1.
#' @param n_subsampled_cells Integer; number of cells per subsample. Will up/downsample automatically.
#' @param subsampling_seed_first Integer; base seed used to ensure reproducibility of each subsampling.
#' @param ... Additional arguments passed to `cytobench::flowSOM_optimal()`.
#'
#' @return
#' The output of `cytobench::flowSOM_optimal()`, typically a list containing the FlowSOM model and cluster assignments.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- flowsom_repeatsubsampling(
#'     ff_list = some_flowset,
#'     outdir = tempdir(),
#'     columns_clustering = c("CD3", "CD4", "CD8"),
#'     n_subsampling = 3,
#'     n_subsampled_cells = 5000
#' )
#' }
flowsom_repeatsubsampling <- function(
    ff_list,
    outdir,
    columns_clustering = c(
        "FITC-A", # CD45RA
        "PE-A", # CCR7
        "ECD-A", # CD28
        "PC5.5-A", # PD1
        "PC7-A", # CD27
        "APC-A", # CD4
        "AF700-A", # CD8
        # "AA750-A",  # CD3   - has been gated for
        "PB-A" # CD57
        # "KrO-A"     # CD45  - has been gated for
    ),
    n_subsampling = 1,
    n_subsampled_cells = 10000,
    subsampling_seed_first = 427764,
    ...) {
    # --- 1. Repeatedly subsample each flowFrame ---
    train_multiple <- subsample_multiple(
        ff_list = ff_list,
        n_subsampling = n_subsampling,
        n_subsampled_cells = n_subsampled_cells,
        subsampling_seed_first = subsampling_seed_first
    )

    # --- 2. Combine all subsamples into one training set ---
    fs_train_multiple <- flowCore::flowSet(train_multiple)

    # --- 3. Run FlowSOM ---
    fs_results <- cytobench::flowSOM_optimal(
        fs_train = fs_train_multiple,
        outdir = outdir,
        colsToUse = columns_clustering,
        ...
    )

    # --- 4. Save results if directory is specified ---
    if (!is.null(outdir)) {
        dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
        qs::qsave(fs_results, file.path(outdir, "r2-FlowSOM_result_train_predictions.qs"))
    }

    return(fs_results)
}
