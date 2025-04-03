cycompare_outcomes_plot <- function(
    df,
    dfcol_grouping_supersamples = c("Study"),
    dfcol_grouping_samples = "Device",
    dfcol_outcomes = c("outcome_1", "outcome_2"),
    dfcol_train_validation_other = "train_validation_test",
    dv_class_positive = c("outcome_1" = "A", "outcome_2" = 5.1),
    results_cycompare_analyse) {
    device_colors <- results_cycompare_analyse[["prepared_data"]][["device_colors"]]
    browser()
    #### 1. Basic plots
    ## 1.1 Outcomes per study
    p1.1 <- plot_outcome_circles(
        df,
        dfcol_grouping_supersamples = dfcol_grouping_supersamples,
        dfcol_outcomes = dfcol_outcomes
    )

    options(
        future.globals.maxSize = 10 * 1024^3
    )
    results_cycompare_analyse[["models"]][["applied"]][["groups"]]
    rocs <- plot_ROC_grouped(
        results_cycompare_analyse,
        dfcol_grouping_supersamples = dfcol_grouping_supersamples,
        dfcol_grouping_samples = dfcol_grouping_samples,
        dfcol_outcomes = dfcol_outcomes,
        dfcol_train_validation_other = dfcol_train_validation_other,
        dv_class_positive = dv_class_positive,
        device_colors = device_colors
    )

    pdf("removeme.pdf")
    rocs[["roc.plot"]]
    dev.off()
    
}
