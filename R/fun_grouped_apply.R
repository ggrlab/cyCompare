#' Apply a Function to Grouped Flow Cytometry Data
#'
#' This function applies a specified analysis function (e.g., FlowSOM training)
#' to subsets of flow cytometry data defined by unique combinations of df columns.
#' Each grouping results in a separate function call, which can be useful
#' for training models per sample, supersample, or condition.
#'
#' @param gated_ff A named list or vector of flowFrames, indexed by the "File" column of df.
#' @param make_flowset Logical; if TRUE, convert subset to a flowSet before applying the function.
#' @param fun Function to apply on each group. Default: `flowsom_repeatsubsampling`.
#' The function should accept a flowSet and an output directory as first two arguments.
#' @param df Data frame containing all following `dfcol_XXX` and `File` column.
#' @param dfcol_grouping_supersamples Character; column name(s) for supersample grouping.
#' @param dfcol_grouping_samples Character; column name(s) for sample grouping.
#' @param dfcol_train_validation_other Character; column name to distinguish training/validation/other.
#' @param dfcol_othergroups Character vector (optional); other metadata columns for grouping.
#' @param outdir_base Character (optional); base directory for writing outputs. The output directory for each group will be
#' created by combining `outdir_base` with the grouping names. E.g.
#' `outdir_base/columnA.valueA___columnB.valueB`.
#' @param verbose Logical; print progress messages.
#' @param ... Additional arguments passed to `fun`.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{results}{List of outputs returned by `fun` for each group.}
#'   \item{groups}{Data frame of the grouping combinations used.}
#' }
#' @export
fun_grouped_apply <- function(
    ff_list,
    result_grouping,
    make_flowset = TRUE,
    fun = flowsom_repeatsubsampling,
    outdir_base = NULL,
    verbose = FALSE,
    return_results = TRUE,
    ...) {
        
    if (!return_results && is.null(outdir_base)) {
        stop("If return_results is FALSE, outdir_base must be specified and the function should save the results there.")
    }

    possible_groupings <- result_grouping[["groups"]]
    if (make_flowset) {
        # Convert gated_ff to a flowSet
        ff_list <- flowCore::flowSet(ff_list)
    }

    # Iterate through each group and apply the analysis function
    # all_results <- future.apply::future_lapply(
    #     future.seed = TRUE,
    #     future.conditions = "message",
    all_results <- lapply(
        seq_len(nrow(possible_groupings)),
        function(grouping_i) {
            grouping_x <- possible_groupings[grouping_i, ]

            if (verbose) {
                message(
                    "Processing grouping: ",
                    paste0(names(grouping_x), ".", grouping_x, collapse = "___"),
                    "\n"
                )
            }

            # Create output directory for this grouping, if specified
            current_outdir <- NULL
            if (!is.null(outdir_base)) {
                current_outdir <- file.path(
                    outdir_base,
                    paste0(names(grouping_x), ".", grouping_x, collapse = "___")
                )
            }

            # Apply the user-specified function
            res <- fun(
                ff_list = ff_list,
                outdir = current_outdir,
                result_grouping[["results"]][[grouping_i]],
                ...
            )
            if (!return_results) {
                # if return_results is FALSE, set res to NULL
                # to avoid returning the results of the function
                # This is useful if the function returns a large object which is
                # SAVED within the function!
                res <- NULL
            }
            return(res)
        }
    )

    # Return results and grouping combinations
    return(
        list(
            "results" = all_results,
            "groups" = possible_groupings
        )
    )
}
