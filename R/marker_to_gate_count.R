#' Extract Positive Populations for Marker-Level MFI Comparisons
#'
#' This function prepares a mapping between markers and their corresponding pre-gated
#' positive populations (as defined by `marker_to_gate`). It filters the provided
#' `dt_count_mfi` data to only include these relevant populations.
#'
#' The conceptual idea is to define, for each marker, one or more positive populations
#' that were gated prior to analysis. These are used to assess changes in marker signal intensity -
#' for example, to evaluate whether median fluorescence intensity (MFI) differs between devices
#' or over time. This helps monitor signal drift or systematic bias in the measurement process.
#'
#' @param marker_to_gate
#' A named list or vector mapping marker names to gating population(s)
#'        (e.g., `list("CD4" = "/Singlets/CD4+", "CD8" = "CD8+")`).
#' @param dt_count_mfi
#' A `data.table` containing MFI and count statistics, including columns
#'        `"pop"` and `"File"` (typically from gated populations).
#'
#' @return A named list with:
#' \describe{
#'   \item{marker_to_gate_dt}{A `data.table` with two columns: `"marker"` and `"pop"` representing the mapping.}
#'   \item{dt_count_mfi_relevant}{A filtered `data.table` with rows corresponding to relevant gated populations.}
#' }
#'
#' @export
marker_to_gate_count <- function(marker_to_gate, dt_count_mfi) {
    browser()
    # Extract unique gating populations from marker_to_gate
    relevant_gates <- unique(unlist(marker_to_gate))

    # Convert marker_to_gate list into a data.table with "marker" and "pop" columns
    marker_to_gate_dt <- lapply(marker_to_gate, data.table::as.data.table) |>
        data.table::rbindlist(idcol = "marker")
    data.table::setnames(marker_to_gate_dt, "V1", "pop") # Rename column

    # Filter dt_count_mfi to keep only relevant gating populations
    dt_count_mfi_relevant <- dt_count_mfi[pop %in% relevant_gates]

    # Check if any gates are missing
    expected <- length(relevant_gates) * nrow(unique(dt_count_mfi[, "File"]))
    if (!nrow(dt_count_mfi_relevant) == expected) {
        missing_gates <- setdiff(relevant_gates, unique(dt_count_mfi_relevant[["pop"]]))
        warning(
            "Not all relevant gates are present in the MFI data: \n  ",
            paste0(missing_gates, collapse = "  \n")
        )
    }

    return(list(
        marker_to_gate_dt = marker_to_gate_dt,
        dt_count_mfi_relevant = dt_count_mfi_relevant
    ))
}
