# Wrapper to
# 1) Load processed data
# 1.2) Subsample the data (multiple times) to the given number of cells
# 2) load the phenodata
# 3) FlowSOM on the training data
flowsom_repeatsubsampling_apply <- function(
    ff_list,
    outdir,
    flowsom_result,
    n_metacluster = 10,
    n_subsampling = 1, # Number of subsamplings
    n_subsampled_cells = 10000, # Number of cells to be subsampled (up or downsampling automatically)
    subsampling_seed_first = 427764,
    keep_clustered_cells = FALSE,
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
    if (!keep_clustered_cells) {
        fs_full_predictedFS[["cells_clusters_from_train"]] <- NULL
    }
    if (exists(outdir)) {
        flowsom_wrapper_saving(
            fs_applied = fs_full_predictedFS,
            dir_out = dirname(outdir),
            modelname = basename(outdir)
        )
    }

    return(fs_full_predictedFS)
}
