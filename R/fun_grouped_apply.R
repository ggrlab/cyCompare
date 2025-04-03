#' Apply a Function to Predefined Groupings of Flow Cytometry Data
#'
#' This function applies a user-defined function (e.g., FlowSOM training) to data
#' across multiple groupings, typically the result of `fun_grouped()`. It allows optional conversion
#' to a `flowSet`, supports parallel execution, and can either return the results or rely on
#' side effects such as file output.
#'
#' @param data
#' The data which should be subset. E.g.
#'  - A named list of flowFrames, indexed by the "File" column of df.
#' @param result_grouping A list returned from `fun_grouped()` that contains `groups` and optionally `results`.
#' @param make_flowset Logical; if TRUE, convert `ff_list` into a `flowSet`. Default is TRUE.
#' @param fun A function to apply to each group. It should accept at least arguments `ff_list`, `outdir`, and a grouping-specific input (e.g., `result_grouping[["results"]][[i]]`).
#' @param outdir_base Optional string path to save results. Required if `return_results` is FALSE.
#' @param verbose Logical; if TRUE, prints progress messages.
#' @param return_results Logical; if FALSE, the function does not return results and assumes side effects (e.g., saving to disk).
#' @param ... Additional arguments passed to `fun`.
#'
#' @return A list containing:
#' \describe{
#'   \item{results}{List of results for each grouping, or NULL if `return_results` is FALSE.}
#'   \item{groups}{Data frame of group combinations processed.}
#' }
#'
#' @export
fun_grouped_apply <- function(
    data,
    result_grouping,
    make_flowset = TRUE,
    fun = flowsom_repeatsubsampling_apply,
    outdir_base = NULL,
    verbose = FALSE,
    return_results = TRUE,
    ...) {
    # Ensure output directory is set when no results are returned

    if (!return_results && is.null(outdir_base)) {
        stop(
            "If return_results is FALSE, outdir_base must be specified and the function should save the results there."
        )
    }

    # Extract group definitions
    is_result_list <- FALSE
    if ("groups" %in% names(result_grouping)) {
        possible_groupings <- result_grouping[["groups"]]
        result_grouping_results <- result_grouping[["results"]]
    } else if (all(sapply(result_grouping, function(x) "groups" %in% names(x)))) {
        # If result_grouping is a list of lists, extract the first element
        all_groupings <- lapply(result_grouping, function(x) x[["groups"]])
        if (!all(sapply(all_groupings, function(x) all(x == all_groupings[[1]])))) {
            stop("All groupings must be identical if a list of grouping results is given.")
        }
        possible_groupings <- all_groupings[[1]]
        result_grouping_results <- lapply(result_grouping, function(x) x[["results"]])
        is_result_list <- TRUE
    } else {
        stop("result_grouping must contain a 'groups' element or be a list of lists with 'groups' elements.")
    }
    # Optionally convert input list to a flowSet
    if (make_flowset) {
        data <- flowCore::flowSet(data)
    }

    # Apply the function to each group in parallel
    all_results <- future.apply::future_lapply(
        seq_len(nrow(possible_groupings)),
        future.seed = TRUE,
        future.conditions = "message",
        function(grouping_i) {
            grouping_x <- possible_groupings[grouping_i, ]

            # Print current grouping if verbose is enabled
            if (verbose) {
                message(
                    "Processing grouping: ",
                    paste0(names(grouping_x), ".", grouping_x, collapse = "___"),
                    "\n"
                )
            }

            # Set output directory for this grouping
            current_outdir <- NULL
            if (!is.null(outdir_base)) {
                current_outdir <- file.path(
                    outdir_base,
                    paste0(names(grouping_x), ".", grouping_x, collapse = "___")
                )
            }

            if (is_result_list) {
                # Get the result for this grouping
                result_grouping_results_x <- lapply(result_grouping_results, function(x) {
                    x[[grouping_i]]
                })
            } else {
                # Get the result for this grouping
                result_grouping_results_x <- result_grouping_results[[grouping_i]]
            }
            # Apply the user-defined function to this group
            res <- fun(
                data,
                outdir = current_outdir,
                grouping = grouping_x,
                result_grouping_results_x,
                ...
            )

            # Drop result if configured to only save via side effect
            if (!return_results) {
                res <- NULL
            }

            return(res)
        }
    )
    # Return all results and the group combinations
    return(
        list(
            "results" = all_results,
            "groups" = possible_groupings
        )
    )
}
