# Wrapper to
# 1) Load processed data
# 1.2) Subsample the data (multiple times) to the given number of cells
# 2) load the phenodata
# 3) FlowSOM on the training data
flowsom_repeatsubsampling <- function(
    gated_ff,
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
    n_subsampling = 1, # Number of subsamplings
    n_subsampled_cells = 10000, # Number of cells to be subsampled (up or downsampling automatically)
    subsampling_seed_first = 427764,
    ...) {


    # Subsample the data
    train_multiple <- unlist(
        lapply(1:n_subsampling, function(i) {
            fs_train_subsampled <- flowCore::fsApply(gated_ff, function(x) {
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
    fs_train_multiple <- flowCore::flowSet(train_multiple)
    fs_results <- cytobench::flowSOM_optimal(
        fs_train = fs_train_multiple,
        outdir = outdir,
        colsToUse = columns_clustering,
        ...
    )
    if (!is.null(outdir)) {
        dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
        qs::qsave(fs_results, file.path(outdir, "r2-FlowSOM_result_train_predictions.qs"))
    }
    return(fs_results)
}
