cycompare_preparation <- function(
    flowframes,
    df,
    ff_columns_relevant,
    device_colors = function(n) {
        RColorBrewer::brewer.pal(n, "Dark2")
    },
    gatingsets,
    gatename_primary,
    max_events_postgate = 10e3) {
    unique_devices <- unique(df[["Device"]])
    if (!all(unique_devices %in% names(device_colors))) {
        if (is.function(device_colors)) {
            # Automatically assign colors
            device_colors <- setNames(
                # RColorBrewer takes a minimum of 3 colors, thus the seq_along
                device_colors(length(unique_devices))[seq_along(unique_devices)],
                unique_devices
            )
            warning("The following device colors were automatically assigned:")
            warning(paste0(device_colors, collapse = ", "))
        } else {
            stop("device_colors must be a named vector or a function(number_of_devices)")
        }
    }


    #### Gate each sample to the primary gate
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
            tmp
        }
    )
    counts_ff <- lapply(gated_ff, function(x) x[["counts"]]) |> data.table::rbindlist()
    data.table::setnames(counts_ff, "sample", "File")
    counts_joint <- counts_ff[df, on = "File"]

    if (quantile(counts_joint[pop == gatename_primary][["count"]], .9) < 100) {
        stop(
            "The primary gate has less than 100 cells in the 90th percentile of samples. Did you select the right gate for these samples?"
        )
    }


    gated_ff <- lapply(gated_ff, function(x) x[["flowset_gated"]][[1]])
    # Limit the number of events post-gating and select only the relevant columns
    gated_ff <- lapply(gated_ff, function(x) {
        nevents <- min(flowCore::nrow(x), max_events_postgate)
        x[sample(flowCore::nrow(x), nevents, replace = FALSE), ff_columns_relevant]
    })

    return(
        list(
            gated_ff = gated_ff,
            counts_joint = counts_joint,
            device_colors = device_colors
        )
    )
}
