#' Apply Gating Hierarchy to a FlowFrame
#'
#' This function applies a predefined gating hierarchy to a `flowFrame` object
#' and extracts the gated population based on a specified gate name.
#'
#' @param ff A `flowFrame` object representing the cytometry data.
#' @param gating A `GatingSet` object containing the gating hierarchy.
#' @param gatename A character string specifying the population to extract.
#'
#' @return A list containing:
#'   \item{ff_gated}{A `flowFrame` object with the gated population.}
#'   \item{counts}{A matrix with population counts from the `GatingSet`.}
#'
#' @details
#' The function first writes the `flowFrame` to a temporary FCS file.
#' It then applies the gating hierarchy from the `GatingSet` to the new FCS file.
#' The function extracts the gated population corresponding to `gatename`
#' and retrieves population counts.
#' @export
gate_flowframe <- function(ff, gating, gatename) {
    # Write the flowFrame to a temporary FCS file in memory
    fcs_filepath <- cytobench::write_memory_FCS(ff)

    # Apply gating hierarchy to the new FCS file
    gating_applied <- suppressMessages(
        flowWorkspace::gh_apply_to_new_fcs(
            gating,
            fcs_filepath
        )
    )

    # Extract the gated population for the specified gate name
    gating_applied_realized <- gating_applied |>
        flowWorkspace::gs_pop_get_data(gatename) |>
        flowWorkspace::realize_view()
    ff_realized <- gating_applied_realized[[1]]
    flowCore::keyword(ff_realized)[["GUID"]] <- flowCore::keyword(ff)[["GUID"]]

    # Retrieve the count matrix from the GatingSet
    counts <- flowWorkspace::gs_pop_get_count_fast(gating_applied)
    counts[, name := flowCore::keyword(ff)[["GUID"]]]

    # Return the gated flowFrame and population counts
    return(
        list(
            ff_gated = ff_realized, # Extract the first gated flowFrame
            counts = counts
        )
    )
}
