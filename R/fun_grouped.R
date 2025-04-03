#' Apply a Function to Grouped Flow Cytometry Data
#'
#' This function applies a specified analysis function (e.g., FlowSOM training)
#' to subsets of flow cytometry data defined by unique combinations of df columns.
#' Each grouping results in a separate function call, which can be useful
#' for training models per sample, supersample, or condition.
#'
#' @param data
#' The data which should be subset. E.g.
#'  - A named list of flowFrames, indexed by the "File" column of df.
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
fun_grouped <- function(
    data,
    make_flowset = TRUE,
    fun = flowsom_repeatsubsampling,
    df,
    dfcol_grouping_supersamples,
    dfcol_grouping_samples,
    dfcol_train_validation_other,
    dfcol_othergroups = NULL,
    outdir_base = NULL,
    verbose = FALSE,
    return_results = TRUE,
    subset_fun = function(data, elements) {
        data[elements]
    },
    ...) {
    if (!return_results && is.null(outdir_base)) {
        stop("If return_results is FALSE, outdir_base must be specified and the function should save the results there.")
    }
    # Get all unique combinations of the specified grouping columns
    possible_groupings <- df |>
        dplyr::select(
            dfcol_grouping_supersamples,
            dfcol_grouping_samples,
            dfcol_train_validation_other,
            dfcol_othergroups
        ) |>
        dplyr::distinct()

    # Iterate through each group and apply the analysis function
    all_results <- future.apply::future_lapply(
        seq_len(nrow(possible_groupings)),
        future.seed = TRUE,
        future.conditions = "message",
        function(grouping_i) {
            grouping_x <- possible_groupings[grouping_i, ]

            if (verbose) {
                message(
                    "Processing grouping: ",
                    paste0(names(grouping_x), ".", grouping_x, collapse = "___"),
                    "\n"
                )
            }

            # Subset metadata to only include rows matching this grouping
            df_subset <- dplyr::left_join(
                grouping_x,
                df,
                by = c(
                    dfcol_grouping_supersamples,
                    dfcol_grouping_samples,
                    dfcol_train_validation_other,
                    dfcol_othergroups
                )
            )

            # Create output directory for this grouping, if specified
            if (!is.null(outdir_base)) {
                current_outdir <- file.path(
                    outdir_base,
                    paste0(names(grouping_x), ".", grouping_x, collapse = "___")
                )
            }

            # Stop if no matching rows found
            if (nrow(df_subset) == 0) {
                stop("No data for this grouping")
            }

            # Extract flowFrames by file names
            data_subset <- subset_fun(data, df_subset[["File"]])

            # Optionally convert to a flowSet
            if (make_flowset) {
                data_subset <- flowCore::flowSet(data_subset)
            }

            # Sanity check: make sure the number of files matches
            if (length(data_subset) != nrow(df_subset)) {
                stop(
                    "Some files are missing from the data_subset set. ",
                    "I wanted to load ", nrow(df_subset),
                    " files, but only found ", length(data_subset),
                    ". Files are loaded by their names ('File' column in df)."
                )
            }

            # Apply the user-specified function
            res <- fun(
                data_subset,
                outdir = current_outdir,
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
