#' Apply a Function to Grouped Flow Cytometry Data
#'
#' Applies a user-specified function to grouped subsets of flow cytometry data.
#' Groupings are defined by combinations of metadata columns in `df`, such as study or device.
#' Each group's subset is passed to `fun()` (e.g., FlowSOM training) independently,
#' allowing per-device or per-supersample processing.
#'
#' @inheritParams cycompare_outcomes_analyse
#' @param data
#' The data which should be subset.
#' A named list of flowFrames (or similar), with names matching entries in `df$File`.
#'
#' @param make_flowset Logical; if `TRUE`, convert the subset to a `flowSet` before calling `fun`.
#' @param fun
#' Function to apply to each group. It must accept a subset of the data and an `outdir` as its first two arguments.
#' @param dfcol_othergroups
#' Optional character vector of additional grouping columns.
#' @param outdir_base
#' Optional character path; used as the root directory for saving results.
#' Each group will have a subdirectory named like `columnA.valueA___columnB.valueB/`.
#' @param verbose Logical; print progress messages (default: `FALSE`).
#' @param return_results
#' Logical; whether to return the result of `fun()` for each group. If `FALSE`,
#' assumes `fun()` is expected to save results to disk into `outdir_base`.
#' @param subset_fun
#' Function to extract a subset from `data` given a vector of `File` names. Defaults to indexing `data[files]`.
#' @param ... Additional arguments passed to `fun`.

#'
#' @return A named list with:
#' \describe{
#'   \item{results}{List of outputs from `fun()` per group (or `NULL` if `return_results = FALSE`).}
#'   \item{groups}{Data frame containing one row per grouping combination used.}
#' }
#'
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
        stop("If return_results is FALSE, you must specify outdir_base to write results.")
    }

    # --- 1. Define all unique group combinations ---
    possible_groupings <- df |>
        dplyr::select(
            dfcol_grouping_supersamples,
            dfcol_grouping_samples,
            dfcol_train_validation_other,
            dfcol_othergroups
        ) |>
        dplyr::distinct()

    # --- 2. Iterate through each group and apply the analysis function ---
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

            # Subset metadata to this group's rows
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

            # Stop if no rows found
            if (nrow(df_subset) == 0) {
                stop("No data for this grouping")
            }

            # --- 3. Define group output path (optional) ---
            if (!is.null(outdir_base)) {
                current_outdir <- file.path(
                    outdir_base,
                    paste0(names(grouping_x), ".", grouping_x, collapse = "___")
                )
            }

            # --- 4. Subset the data based on df_subset[["File"]] ---
            data_subset <- subset_fun(data, df_subset[["File"]])

            # Convert to flowSet if requested
            if (make_flowset) {
                data_subset <- flowCore::flowSet(data_subset)
            }

            # --- 5. Sanity check on subset size ---
            if (length(data_subset) != nrow(df_subset)) {
                stop(
                    "Mismatch: expected ", nrow(df_subset), " files but found ",
                    length(data_subset), " in subset."
                )
            }

            # --- 6. Apply the user-defined function ---
            res <- fun(
                data_subset,
                outdir = current_outdir,
                ...
            )

            # Discard result if not requested (e.g., result was saved to disk)
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

    # --- 7. Return all results with corresponding group metadata ---
    return(
        list(
            "results" = all_results,
            "groups" = possible_groupings
        )
    )
}
