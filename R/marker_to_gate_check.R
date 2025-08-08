#' Check Marker-to-Gate Mapping Against GatingSets
#'
#' Validates whether all gates specified in a `marker_to_gate` mapping are present in the provided `GatingSet` objects.
#' The `marker_to_gate` mapping is used to define positive populations for each marker, which are
#' expected to be pre-gated. These positive populations serve as reference points to assess whether the
#' median fluorescence intensity (MFI) of markers shifts (e.g. over time or between devices).
#'
#' @param marker_to_gate
#' A named character vector mapping each marker to one or more population gate names.
#' These gates are assumed to represent positively stained reference populations.
#' @param gatingsets
#' A named list of `GatingSet` objects, one per sample or FCS file.
#' Names must match those used elsewhere in the analysis.
#'
#' @return
#' No return value. Warnings are issued if any expected gates are missing from the corresponding GatingSet.
#'
#' @details
#' This function is a diagnostic tool that ensures the integrity of reference populations used for downstream
#' MFI-based comparisons. Missing gates will be reported with sample-specific warnings.
#'
#' @export
marker_to_gate_check <- function(marker_to_gate, gatingsets) {
    # Extract all population paths from each GatingSet
    all_paths <- lapply(gatingsets, flowWorkspace::gs_get_pop_paths)

    # Get unique gates referenced in the marker_to_gate mapping
    relevant_gates <- unique(unlist(marker_to_gate))

    # Check for missing gates in each GatingSet
    relevant_gates_in_gatingsets <- sapply(names(all_paths), simplify = FALSE, function(x) {
        missing_gates <- setdiff(relevant_gates, all_paths[[x]])
        return(missing_gates)
    })

    # Warn if any gates are missing
    if (any(sapply(relevant_gates_in_gatingsets, length) > 0)) {
        warning("Not all relevant gates are present in the GatingSet:\n")
        for (x in names(relevant_gates_in_gatingsets)) {
            if (length(relevant_gates_in_gatingsets[[x]]) > 0) {
                warning(paste0(
                    "GatingSet of sample ", x, " is missing the following gates:\n",
                    paste0(relevant_gates_in_gatingsets[[x]], collapse = "\n")
                ))
            }
        }
    }
}
