#' Prepare Gated FlowFrames and Metadata for Cross-Device Comparison
#'
#' This function applies a primary gating step to flow cytometry data,
#' assigns colors to devices, and prepares both gated events and metadata
#' for downstream analyses such as visualization or normalization.
#'
#' @param flowframes A named list of `flowFrame` objects. Names must match sample identifiers in `df`.
#' @param df A data.frame containing metadata for each sample, including a "Device" column.
#' @param ff_columns_relevant A character vector of column names to retain after gating (e.g., marker channels).
#' @param device_colors Either a named character vector mapping device names to colors, or a function `f(n)` returning `n` colors.
#'                      Defaults to using `RColorBrewer::brewer.pal(n, "Dark2")`.
#' @param gatingsets A named list of `GatingSet` objects used for gating, one per sample.
#' @param gatename_primary The name of the primary gate used to subset events (e.g., "Live").
#' @param n_events_postgate Exact number of events to retain per sample after gating. Might upsample cells. Default is 10,000. See `cytoBench::subsample_ff`.
#' @param seed Random seed for reproducibility when downsampling. Default is 42.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{gated_ff}{A list of gated and downsampled `flowFrame` objects, keeping only relevant columns.}
#'   \item{counts_joint}{A data.table with per-sample gated cell counts merged with metadata.}
#'   \item{device_colors}{Named vector of device colors used, auto-assigned if a function was provided.}
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
    dfcol_grouping_samples = "Device",
    transformlist = NULL) {
    mandatory_df_cols <- c("File", "SuperSample", "Sample", dfcol_grouping_samples[[1]], "Time")
    if (!all(mandatory_df_cols %in% colnames(df))) {
        tmp <- paste0(
            "df must contain the following columns: ",
            paste(mandatory_df_cols, collapse = ", "), ". Missing:",
            paste(mandatory_df_cols[!mandatory_df_cols %in% colnames(df)], collapse = ", ")
        )
        stop(tmp)
    }
    if (length(dfcol_grouping_samples) > 1) {
        warning(
            "dfcol_grouping_samples is more than one element, colors will be assigned according to the FIRST element only: ",
            dfcol_grouping_samples[1]
        )
    }
    # Identify all unique devices in metadata
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
                "The following device colors were automatically assigned:",
                paste0(
                    names(device_colors),
                    " = ",
                    device_colors,
                    collapse = ", "
                )
            )
        } else {
            stop("device_colors must be a named vector or a function(number_of_devices)")
        }
    }
    if (all(is.null(gatingsets))) {
        gatingsets <- sapply(
            names(flowframes),
            simplify = FALSE,
            function(x) {
                tmpfile <- cytobench::write_memory_FCS(
                    flowframes[[x]][1, ]
                )

                empty_gs <- flowWorkspace::GatingSet(flowCore::flowSet(flowframes[[x]][1, ]))
                return(empty_gs)
            }
        )
        gatename_primary <- "root"
        if (!all(is.null(marker_to_gate))) {
            # Now now gates exist, so marker_to_gate can only be "root" for all markers
            marker_to_gate[TRUE] <- "root"
        }
    }

    marker_to_gate_check(
        marker_to_gate = marker_to_gate,
        gatingsets = gatingsets
    )

    # Apply the primary gate to each sample using cytobench::gate_cells
    gated_ff <- sapply(
        names(flowframes),
        simplify = FALSE,
        function(x) {
            tmp <- cytobench::gate_cells(
                flowset = flowCore::flowSet(flowframes[[x]]),
                gatingset = gatingsets[[x]],
                gatename = gatename_primary,
                verbose = FALSE
            )
            tmp[["counts"]][["sample"]] <- x
            # markernames to colnames in gated_ff
            # Such that transformlist later works by COLUMN name
            names_dict <- setNames(
                names(flowCore::markernames(tmp[["flowset_gated"]])),
                flowCore::markernames(tmp[["flowset_gated"]])
            )
            colnames(tmp[["counts"]])[colnames(tmp[["counts"]]) %in% names(names_dict)] <- na.omit(names_dict[colnames(tmp[["counts"]])])
            tmp
        }
    )
    gated_ff <- lapply(gated_ff, function(x) x[, ff_columns_relevant])

    if (length(transformlist) == 1 && !is.null(transformlist)) {
        if (is.function(transformlist)) {
            transformlist <- list(transformlist)
        }
        transformlist <- setNames(
            rep(transformlist, length(ff_columns_relevant)),
            ff_columns_relevant
        )
    }

    if (!all(ff_columns_relevant %in% names(transformlist))) {
        stop("All ff_columns_relevant must be present in transformlist")
    }

    # Collect and join cell counts with metadata
    counts_ff <- lapply(gated_ff, function(x) x[["counts"]]) |> data.table::rbindlist(fill = TRUE)
    data.table::setnames(counts_ff, "sample", "File")
    if (!all(df[["File"]] %in% counts_ff[["File"]])) {
        stop("Not all files in df are present in counts_ff based on 'File' column.")
    }
    counts_joint <- data.table::data.table(df)[counts_ff, on = "File"]

    # Sanity check: ensure gate captures sufficient events
    if (quantile(counts_joint[pop == gatename_primary][["count"]], .9) < 100) {
        stop(
            "The primary gate has less than 100 cells in the 90th percentile of samples. ",
            "Did you select the right gate for these samples?"
        )
    }

    # Extract gated flowFrames
    gated_ff <- lapply(gated_ff, function(x) x[["flowset_gated"]][[1]])

    # Downsample and select only relevant columns
    gated_ff <- lapply(gated_ff,
        cytobench::subsample_ff,
        n_cells = n_events_postgate,
        seed = seed
    )
    gated_ff <- lapply(gated_ff, function(x) x[, ff_columns_relevant])


    if (length(transformlist) == 1 && !is.null(transformlist)) {
        if (is.function(transformlist)) {
            transformlist <- list(transformlist)
        }
        transformlist <- setNames(
            rep(transformlist, length(ff_columns_relevant)),
            ff_columns_relevant
        )
    }

    if (!all(ff_columns_relevant %in% names(transformlist))) {
        stop("All ff_columns_relevant must be present in transformlist")
    }

    # Return result
    return(
        list(
            gated_ff = gated_ff,
            counts_joint = counts_joint,
            device_colors = device_colors,
            marker_to_gate = marker_to_gate,
            gatename_primary = gatename_primary
        )
    )
}
