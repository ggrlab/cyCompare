fun_grouped <- function(
    gated_ff,
    make_flowset = TRUE,
    fun = flowsom_repeatsubsampling,
    df,
    dfcol_grouping_supersamples,
    dfcol_grouping_samples,
    dfcol_train_validation_other,
    dfcol_othergroups = NULL,
    outdir_base = NULL,
    verbose = FALSE,
    ...) {
    possible_groupings <- df |>
        dplyr::select(
            dfcol_grouping_supersamples,
            dfcol_grouping_samples,
            dfcol_train_validation_other,
            dfcol_othergroups
        ) |>
        dplyr::distinct()
    all_results <- list()
    # For each grouping (supersample and sample), train a separate FlowSOM
    for (grouping_i in seq_len(nrow(possible_groupings))) {
        grouping_x <- possible_groupings[grouping_i, ]
        if (verbose) {
            cat(
                "Processing grouping: ",
                paste0(
                    names(grouping_x),
                    ".",
                    grouping_x,
                    collapse = "___"
                ),
                "\n"
            )
        }
        # Effectively select all rows from df that match the grouping
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
        if (!is.null(outdir_base)) {
            current_outdir <- file.path(
                outdir_base,
                paste0(
                    names(grouping_x),
                    ".",
                    grouping_x,
                    collapse = "___"
                )
            )
        }

        if (nrow(df_subset) == 0) {
            stop("No data for this grouping")
        }
        fs_subset <- gated_ff[df_subset[["File"]]]
        if (make_flowset) {
            fs_subset <- flowCore::flowSet(fs_subset)
        }
        if (length(fs_subset) != nrow(df_subset)) {
            stop(
                "Some files are missing from the fs_subset set. I wanted to load ",
                nrow(df_subset), " files, but only found ", length(fs_subset),
                ". Files are loaded by their names ('File' column in df)."
            )
        }
        all_results[[grouping_i]] <- fun(
            fs_subset,
            current_outdir,
            ...
        )
    }
    return(
        list(
            "results" = all_results,
            "groups" = possible_groupings
        )
    )
}
