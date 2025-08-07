#' Build a Named List of Transformation Functions for Each Marker
#'
#' This utility constructs a named list of transformation functions to be used with `flowCore::transformList()`.
#' If a single transformation function is provided, it is replicated and named for each marker in `relevant_columns`.
#' If `flowcore = TRUE`, it returns a `flowCore::transformList()` object instead.
#'
#' @param transformlist Either a single function, a list of functions, or a named list of functions to apply.
#' @param relevant_columns Character vector of marker names to apply transformations to.
#' @param flowcore Logical. If `TRUE`, returns a `flowCore::transformList` object. Default is `FALSE`.
#'
#' @return A named list of transformation functions, or a `flowCore::transformList` if `flowcore = TRUE`.
#'
#' @export
#'
#' @examples
#' # Example 1: Single function applied to multiple markers
#' transformlist_named(function(x) asinh(x / 1000), c("CD4", "CD8"))
#'
#' # Example 2: FlowCore-compatible object
#' transformlist_named(function(x) asinh(x / 1000), c("CD4", "CD8"), flowcore = TRUE)
transformlist_named <- function(transformlist, relevant_columns, flowcore = FALSE) {
    # If a single unnamed function is provided, replicate and name it for all columns
    if (length(transformlist) == 1 && !is.null(transformlist)) {
        if (is.function(transformlist)) {
            transformlist <- list(transformlist)
        }
        transformlist <- setNames(rep(transformlist, length(relevant_columns)), relevant_columns)
    }

    # Return either a flowCore::transformList or the named list
    if (flowcore) {
        return(
            flowCore::transformList(
                relevant_columns,
                transformlist
            )
        )
    } else {
        return(transformlist)
    }
}
