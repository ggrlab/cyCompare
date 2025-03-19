cycompare <- function(
    flowframes,
    df,
    ff_columns_relevant,
    gatingsets,
    gatename_primary,
    marker_to_gate,
    outcome_columns_df,
    outcome_models) {
    # browser()

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
}
