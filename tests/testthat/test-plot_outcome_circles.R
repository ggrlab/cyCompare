test_that("Test plotting", {
    df1 <- data.frame(
        batch = rep(c("A", "B"), each = 50),
        device = rep(c("X", "Y"), times = 50),
        outcome1 = sample(c("pos", "neg"), 100, replace = TRUE),
        outcome2 = sample(c("high", "low"), 100, replace = TRUE)
    )

    # usual usage:
    plots <- plot_outcome_circles(
        df1,
        dfcol_grouping_supersamples = c("batch"),
        dfcol_outcomes = c("outcome1", "outcome2")
    )
    # "device" would usually be a SAMPLE grouping variable, not a SUPERSAMPLE
    # but here its just for testing
    plots <- plot_outcome_circles(
        df1,
        dfcol_grouping_supersamples = c("batch", "device"),
        dfcol_outcomes = c("outcome1", "outcome2")
    )
    print(plots)
    testthat::expect_true(TRUE)
})
