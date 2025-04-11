marker_to_gate_count <- function(marker_to_gate, dt_count_mfi) {
    # Extract unique gating populations from marker_to_gate
    relevant_gates <- unique(unlist(marker_to_gate))

    # Convert marker_to_gate list into a data.table with "marker" and "pop" columns
    marker_to_gate_dt <- lapply(marker_to_gate, data.table::as.data.table) |>
        data.table::rbindlist(idcol = "marker")
    data.table::setnames(marker_to_gate_dt, "V1", "pop") # Rename column

    # Filter dt_count_mfi to keep only relevant gating populations
    dt_count_mfi_relevant <- dt_count_mfi[pop %in% relevant_gates]
    if (!nrow(dt_count_mfi_relevant) == length(relevant_gates) * nrow(unique(dt_count_mfi[, "File"]))) {
        missing_gates <- setdiff(relevant_gates, unique(dt_count_mfi_relevant[["pop"]]))
        warning("Not all relevant gates are present in the MFI data: \n  ", paste0(missing_gates, collapse = "  \n"))
    }
    return(list(
        marker_to_gate_dt = marker_to_gate_dt,
        dt_count_mfi_relevant = dt_count_mfi_relevant
    ))
}
