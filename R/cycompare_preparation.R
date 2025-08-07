#' Prepare Gated FlowFrames and Metadata for Cross-Device Comparison
#'
#' Applies a primary gating step to flow cytometry samples, assigns colors to devices,
#' and returns gated data and metadata for downstream clustering, modeling, or plotting.
#'
#' @inheritParams cycompare_outcomes_analyse
#' @param seed Integer.
#' Random seed for reproducibility during subsampling. Default is 42.
#' @param marker_to_gate
#' Named vector mapping markers to their associated gate name.
#' Later used for plotting.
#' @param transformlist
#' A named list or function defining transformations for `ff_columns_relevant`.
#' All markers listed must be present. If provided as a function or unnamed list, it is recycled.
#'
#' @return A named list with the following elements:
#' \describe{
#'   \item{gated_ff}{Named list of gated and downsampled `flowFrame`s, reduced to `ff_columns_relevant`.}
#'   \item{counts_joint}{Merged `data.table` combining per-sample cell counts (from gating) and metadata (`df`).}
#'   \item{device_colors}{Named vector of device color assignments. Auto-generated if not provided explicitly.}
#'   \item{marker_to_gate}{The (possibly adjusted) mapping from markers to gates used during validation.}
#'   \item{gatename_primary}{Name of the gate used to subset events (typically `"Live"` or `"CD3+"`).}
#' }
#'
#' @export
cycompare_preparation <- function(
    flowframes,
    df,
    ff_columns_relevant,
    device_colors = function(n) {
        RColorBrewer::brewer.pal(n, "Dark2")
    },
    gatingsets,
    gatename_primary,
    n_events_postgate = 10e3,
    marker_to_gate = NULL,
    seed = 42,
    dfcol_grouping_supersamples = c("Study"),
    dfcol_grouping_samples = "Device",
    dfcol_train_validation_other = NULL,
    transformlist = NULL) {
    # --- 1. Check required columns in metadata ---
    mandatory_df_cols <- c(
        "File", "SuperSample", "Sample", dfcol_grouping_samples[[1]],
        dfcol_train_validation_other, dfcol_grouping_supersamples, "Time"
    )
    if (!all(mandatory_df_cols %in% colnames(df))) {
        stop(
            "Missing required columns in `df`: ",
            paste(setdiff(mandatory_df_cols, colnames(df)), collapse = ", ")
        )
    }

    # Warn if multiple grouping columns for samples are provided
    if (length(dfcol_grouping_samples) > 1) {
        warning(
            "Only the first sample grouping column is used for device color assignment: ",
            dfcol_grouping_samples[1]
        )
    }

    # --- 2. Identify all unique devices in metadata ---
    unique_devices <- unique(df[[dfcol_grouping_samples[[1]]]])
    if (all(is.null(names(flowframes)))) {
        stop("flowframes must be a named list")
    }

    # Assign colors if a function is provided instead of a named vector
    if (!all(unique_devices %in% names(device_colors))) {
        if (is.function(device_colors)) {
            # Auto-assign colors (ensure at least 3 colors for RColorBrewer)
            device_colors <- setNames(
                device_colors(length(unique_devices))[seq_along(unique_devices)],
                unique_devices
            )
            warning(
                "Device colors were auto-assigned:\n",
                paste(names(device_colors), device_colors, sep = " = ", collapse = ", ")
            )
        } else {
            stop("`device_colors` must be a named vector or a function(number_of_devices)")
        }
    }

    # --- 3. Set up fallback gating if none is provided ---
    if (all(is.null(gatingsets))) {
        gatingsets <- sapply(
            names(flowframes),
            simplify = FALSE,
            function(x) {
                empty_gs <- flowWorkspace::GatingSet(flowCore::flowSet(flowframes[[x]][1, ]))
                return(empty_gs)
            }
        )
        gatename_primary <- "root"
        if (!all(is.null(marker_to_gate))) {
            # Now that gates exist, so marker_to_gate can only be "root" for all markers
            marker_to_gate[TRUE] <- "root"
        }
    }

    # --- 4. Validate marker-to-gate mapping ---
    marker_to_gate_check(marker_to_gate, gatingsets)
    missing_markers <- setdiff(names(marker_to_gate), flowCore::colnames(flowframes[[1]]))
    if (length(missing_markers) > 0) {
        stop(
            "Markers in `marker_to_gate` missing from flowframes. ",
            "Please use the colnames, not markernames.\nMissing: ",
            paste(missing_markers, collapse = ", "),
            "\nPresent: ",
            paste0(flowCore::colnames(flowframes[[1]]), collapse = ", ")
        )
    }

    # --- 5. Apply gating ---
    # Apply the primary gate to each sample using cytobench::gate_cells
    gated_ff <- sapply(
        names(flowframes),
        simplify = FALSE,
        function(x) {
            gated <- cytobench::gate_cells(
                flowset = flowCore::flowSet(flowframes[[x]]),
                gatingset = gatingsets[[x]],
                gatename = gatename_primary,
                verbose = FALSE
            )
            # Track which sample was processed
            gated[["counts"]][["sample"]] <- x

            # Use colnames instead of markernames in count table (ensures transformation works)
            name_map <- setNames(
                names(flowCore::markernames(gated[["flowset_gated"]])),
                flowCore::markernames(gated[["flowset_gated"]])
            )

            colnames(gated[["counts"]])[colnames(gated[["counts"]]) %in% names(name_map)] <- na.omit(name_map[colnames(gated[["counts"]])])
            return(gated)
        }
    )

    # --- 6. Merge gating counts with metadata ---
    counts_ff <- lapply(gated_ff, function(x) x[["counts"]]) |> data.table::rbindlist(fill = TRUE)
    data.table::setnames(counts_ff, "sample", "File")
    if (!all(df[["File"]] %in% counts_ff[["File"]])) {
        stop("Mismatch between `df$File` and gated counts")
    }
    counts_joint <- data.table::data.table(df)[counts_ff, on = "File"]

    # Sanity check: were any gates nearly empty?
    if (quantile(counts_joint[pop == gatename_primary][["count"]], 0.9) < 100) {
        stop("Primary gate yields <100 events in 90th percentile - check `gatename_primary`.")
    }

    # --- 7. Extract gated flowFrames ---
    gated_ff <- lapply(gated_ff, function(x) x[["flowset_gated"]][[1]])

    # --- 8. Downsample and subset markers ---
    gated_ff <- lapply(gated_ff, cytobench::subsample_ff,
        n_cells = n_events_postgate, seed = seed
    )
    # Retain only relevant marker columns
    gated_ff <- lapply(gated_ff, function(x) x[, ff_columns_relevant])

    # --- 9. process transformlist if given ---
    # transformlist_named() is going to be called where needed again!
    # This here is just for testing if the transformlist is valid
    transformlist_list <- transformlist_named(
        transformlist,
        relevant_columns = ff_columns_relevant,
        flowcore = FALSE
    )
    # Ensure every marker has a transform
    if (!all(ff_columns_relevant %in% names(transformlist_list))) {
        stop(
            "All `ff_columns_relevant` must be covered by `transformlist`. Missing: ",
            paste(setdiff(ff_columns_relevant, names(transformlist_list)), collapse = ", ")
        )
    }

    # --- 10. Return prepared object ---
    return(
        list(
            gated_ff = gated_ff, # List of gated, downsampled, transformed flowFrames
            counts_joint = counts_joint, # Combined metadata and per-sample gating statistics
            device_colors = device_colors, # Named vector of colors per device
            marker_to_gate = marker_to_gate, # Final validated marker-to-gate mapping
            gatename_primary = gatename_primary # The primary gating population used
        )
    )
}
