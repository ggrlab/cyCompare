#' Apply a Function to Predefined Groupings of Flow Cytometry Data
#'
#' Applies a user-defined function (e.g., FlowSOM application or modeling) to each group from `fun_grouped()` output.
#' This supports grouped application of downstream tasks across supersample/device combinations.
#'
#' @inheritParams cycompare_outcomes_analyse
#' @inheritParams fun_grouped
#' @param data
#' Named list of flowFrames or similar objects. Names must match `File` column from metadata.
#' @param result_grouping
#' A list returned from `fun_grouped()`, containing `groups` and optionally `results`.
#' Can also be a list of such lists (e.g., clustering + modeling).
#'
#' @return A named list with:
#' \describe{
#'   \item{results}{List of outputs from `fun()` per group (or `NULL` if `return_results = FALSE`).}
#'   \item{groups}{Data frame of group combinations used.}
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
    # Require outdir_base if results are not returned (assumes side effects)
    if (!return_results && is.null(outdir_base)) {
        stop("If return_results is FALSE, outdir_base must be specified.")
    }

    # --- Parse result_grouping input ---
    is_result_list <- FALSE
    if ("groups" %in% names(result_grouping)) {
        # Standard case: one set of groupings
        possible_groupings <- result_grouping[["groups"]]
        result_grouping_results <- result_grouping[["results"]]
    } else if (all(sapply(result_grouping, function(x) "groups" %in% names(x)))) {
        # If result_grouping is a list of lists, extract the first element
        all_groupings <- lapply(result_grouping, function(x) x[["groups"]])
        if (!all(sapply(all_groupings, function(x) all(x == all_groupings[[1]])))) {
            stop("All groupings must be identical if a list of result_groupings is given.")
        }
        possible_groupings <- all_groupings[[1]]
        result_grouping_results <- lapply(result_grouping, function(x) x[["results"]])
        is_result_list <- TRUE
    } else {
        stop("result_grouping must contain a 'groups' element or be a list of lists with 'groups' elements.")
    }

    # --- Optional conversion to flowSet ---
    if (make_flowset) {
        data <- flowCore::flowSet(data)
    }

    # --- Apply user function per group in parallel ---
    all_results <- future.apply::future_lapply(
        seq_len(nrow(possible_groupings)),
        future.seed = TRUE,
        future.conditions = "message",
        function(grouping_i) {
            grouping_x <- possible_groupings[grouping_i, ]

            if (verbose) {
                # Print current grouping if verbose is enabled
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

            # Extract group-specific result(s)
            if (is_result_list) {
                result_grouping_results_x <- lapply(result_grouping_results, function(x) {
                    x[[grouping_i]]
                })
            } else {
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

    # --- Return all results and the group combinations ---
    return(
        list(
            results = all_results,
            groups = possible_groupings
        )
    )
}
