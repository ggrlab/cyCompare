test_that("CyCompare coords_helper", {
    data(aSAH, package = "pROC")
    r <- pROC::roc(aSAH$outcome, aSAH$s100b,
        levels = c("Good", "Poor")
    )
    tmp <- coords_helper(proc_result = r)
    tmp2 <- coords_helper(proc_result = r, threshold = 0.5)
    testthat::expect_equal(tmp$proc, tmp2$proc)
})
