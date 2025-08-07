#' Calculate Optimal Transport Distance (OTD) to a Reference Sample
#'
#' This function transforms gated flow cytometry data, constructs a reference sample,
#' and computes the Optimal Transport Distance (OTD) from each sample to the reference.
#'
#' @param ff_gated A list of `flowFrame` objects representing gated samples.
#' @param transformlist A single function or a named list of functions used for transformation.
#' @param n_referencesample Integer. Number of events to subsample for the reference sample.
#' @param kwargs_loss Named list of additional arguments to pass to `loss_pairwise()`.
#' @param relevant_columns Character vector of column names to use. If `NULL`, all markers are used.
#'
#' @return A `data.frame` with one row per sample and a column `OTD_to_referencesample`.
#'         The first column is `"File"`, matching the input metadata table.
#'
#' @export
otd_to_reference <- function(ff_gated,
                             transformlist,
                             n_referencesample = 1e4,
                             kwargs_loss = list(
                                 loss = lossfun_hist,
                                 verbose = FALSE,
                                 write_intermediate = FALSE,
                                 should_skip = function(i, j) FALSE,
                                 take_time = FALSE,
                                 return_as_matrix = TRUE
                             ),
                             relevant_columns = NULL) {
    # Default: use all marker columns
    if (is.null(relevant_columns)) {
        relevant_columns <- flowCore::colnames(ff_gated[[1]])
    }

    # Convert list to flowSet
    gated_fs <- flowCore::flowSet(ff_gated)

    # Build transformList
    fc_transformlist <- flowCore::transformList(
        relevant_columns,
        if (length(transformlist) == 1) rep(transformlist, length(relevant_columns)) else transformlist
    )

    # Apply transformation
    gated_fs_transformed <- flowCore::transform(gated_fs, fc_transformlist)

    # Convert flowFrames to data.tables
    dt_events <- flowCore::fsApply(
        gated_fs_transformed,
        simplify = FALSE,
        function(x) flowCore::exprs(x) |> data.table::as.data.table()
    )

    # Create reference sample from all events
    dt_events_bound <- data.table::rbindlist(dt_events, idcol = "File")
    referencesample <- dt_events_bound[sample(.N, n_referencesample), ][, File := NULL]

    # Compute pairwise distance from reference sample to each sample
    loss_calculated <- do.call(loss_pairwise, c(
        list(
            datalist_A = list("referencesample" = referencesample),
            datalist_B = dt_events
        ), kwargs_loss
    ))

    # Convert to tidy data.frame
    dists_df <- data.frame(t(loss_calculated[["dist"]])) |>
        tibble::rownames_to_column("File") |>
        dplyr::rename("OTD_to_referencesample" = referencesample)

    return(dists_df)
}
