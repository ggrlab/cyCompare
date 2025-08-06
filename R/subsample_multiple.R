#' Perform Repeated Subsampling of Flow Cytometry Data
#'
#' Applies repeated subsampling to a list or `flowSet` of flowFrames. Each subsample iteration
#' uses a reproducible seed and appends a unique `_subsampledX` suffix to the sample names.
#'
#' @param ff_list A list or `flowSet` of `flowFrame` objects to be subsampled.
#' @param n_subsampling Integer. Number of repeated subsampling rounds to perform (default: 1).
#' @param n_subsampled_cells Integer. Number of cells per subsample. Automatically downsamples or upsamples.
#' @param subsampling_seed_first Integer. Base seed used for reproducibility across iterations.
#'
#' @return A list of subsampled `flowFrame` objects, one per sample per iteration. Sample names are suffixed accordingly.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' subsampled_frames <- subsample_multiple(
#'     ff_list = my_flowset,
#'     n_subsampling = 3,
#'     n_subsampled_cells = 5000,
#'     subsampling_seed_first = 123
#' )
#' }
subsample_multiple <- function(
    ff_list,
    n_subsampling = 1,
    n_subsampled_cells = 10000,
    subsampling_seed_first = 427764) {
    # Perform subsampling independently per iteration
    unlist(
        lapply(seq_len(n_subsampling), function(i) {
            fs_subsampled <- flowCore::fsApply(ff_list, function(x) {
                cytobench::subsample_ff(
                    x,
                    n_cells = n_subsampled_cells,
                    seed = subsampling_seed_first + i - 1
                )
            })

            # Label each sample with subsampling round
            flowCore::sampleNames(fs_subsampled) <- paste0(
                flowCore::sampleNames(fs_subsampled),
                "_subsampled", i
            )

            flowCore::flowSet_to_list(fs_subsampled)
        })
    )
}
