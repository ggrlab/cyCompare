cycompare <- function(
    flowframes,
    df,
    ff_columns_relevant,
    gatingsets,
    gatename_primary,
    marker_to_gate,
    outcome_columns_df,
    outcome_models,
    device_colors = function(n) {
        RColorBrewer::brewer.pal(n, "Dark2")
    }) {
    unique_devices <- unique(df[["Device"]])
    if (!all(unique_devices %in% names(device_colors))) {
        if (is.function(device_colors)) {
            # Automatically assign colors
            device_colors <- setNames(
                # RColorBrewer takes a minimum of 3 colors, thus the seq_along
                device_colors(length(unique_devices))[seq_along(unique_devices)],
                unique_devices
            )
            print("The following device colors were automatically assigned:")
            print(device_colors)
        } else {
            stop("device_colors must be a named vector or a function(number_of_devices)")
        }
    }

    #### 1. Basic plots
    ## 1.1 Samples over time per device
    p1.1 <- plot_samples_by_time(df)


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

    # Plots of counts and percentages
    p1.2_3 <- plot_counts(counts_joint, gatename_primary, device_colors)

    #### Figure 2 plots
    # 1. Positive population MFI
    p2.1 <- plot_MFI_positivegates(
        dt_count_mfi = counts_joint,
        marker_to_gate = marker_to_gate,
        device_colors = device_colors,
        transformlist = function(x) {
            asinh(x / 1e3)
        }
    )
    browser()
    gated_ff <- lapply(gated_ff, function(x) x[["ff_gated"]])
}
