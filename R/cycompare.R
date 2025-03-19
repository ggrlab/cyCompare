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
            device_colors <- setNames(
                device_colors(length(unique_devices)),
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
            gate_flowframe(
                ff = flowframes[[x]],
                gating = gatingsets[[x]],
                gatename = gatename_primary
            )
        }
    )
    browser()
    counts_ff <- lapply(gated_ff, function(x) x[["counts"]]) |> data.table::rbindlist(idcol = "File")
    counts_joint <- counts_ff[df, on = "File"]
    plot_counts(counts_joint, gatename_primary)

    gated_ff <- lapply(gated_ff, function(x) x[["ff_gated"]])
}
