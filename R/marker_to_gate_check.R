
marker_to_gate_check <- function(marker_to_gate, gatingsets) {
    all_paths <- lapply(gatingsets, flowWorkspace::gs_get_pop_paths)
    relevant_gates <- unique(unlist(marker_to_gate))
    relevant_gates_in_gatingsets <- sapply(names(all_paths), simplify = FALSE, function(x) {
        missing_gates <- setdiff(relevant_gates, all_paths[[x]])

        return(missing_gates)
    })
    relevant_gates_in_gatingsets[[1]] <- c("a", "bv")
    if (any(sapply(relevant_gates_in_gatingsets, length) > 0)) {
        warning("Not all relevant gates are present in the GatingSet: \n  ")
        for (x in names(relevant_gates_in_gatingsets)) {
            if (length(relevant_gates_in_gatingsets[[x]]) > 0) {
                warning(paste0(
                    "GatingSet of sample ", x, " is missing the following gates: \n",
                    paste0(
                        relevant_gates_in_gatingsets[[x]],
                        collapse = "\n"
                    )
                ))
            }
        }
    }
}
