test_that("CyCompare", {
    n_events_postgate <- 499
    tmp <- prepare_testdata()
    tmp_cy <- cycompare(
        flowframes = flowCore::flowSet_to_list(tmp[["fs"]]),
        ff_columns_relevant = names(tmp[["relevant_mn"]]),
        df = tmp[["df"]],
        n_events_postgate = n_events_postgate,
        gatingsets = tmp[["gatingsets"]],
        gatename_primary = "/FITC_PE_gate",
        marker_to_gate = tmp[["marker_to_gate"]]
    )
    tmpdir <- cytobench::local_tempdir_time()
    pdf(file.path(tmpdir, "removeme.pdf"), width = 20)
    print(tmp_cy)
    dev.off()
    # The events post gating should be n_events_postgate.
    # Therefore the total number of cells should be
    # number of samples * n_events_postgate
    testthat::expect_equal(
        nrow(tmp_cy$results_flowsom$cells_clusters_from_train),
        length(tmp$fs) * n_events_postgate
    )
    df_result1 <- data.frame(
        "File" = c(
            "StudyA_d01_s01_CyA",
            "StudyA_d02_s02_CyA",
            "StudyA_d03_s03_CyA",
            "StudyA_d01_s01_CyB",
            "StudyA_d02_s02_CyB",
            "StudyA_d03_s03_CyB"
        ),
        "OTD_to_referencesample" = c(
            0.002165620,
            0.002185727,
            0.002199899,
            0.002504403,
            0.002519673,
            0.002542104
        )
    )
    testthat::expect_equal(
        df_result1,
        tmp_cy$calculated_otd,
        tolerance = 1e-7
    )
})
test_that("CyCompare, no gatingset", {
    n_events_postgate <- 499
    tmp <- prepare_testdata()
    tmp_cy <- cycompare(
        flowframes = flowCore::flowSet_to_list(tmp[["fs"]]),
        ff_columns_relevant = names(tmp[["relevant_mn"]]),
        df = tmp[["df"]],
        n_events_postgate = n_events_postgate,
        gatingsets = NULL,
        # If the gatingsets are NULL; gatename_primary is replaced by "root"
        gatename_primary = "/FITC_PE_gate",
        marker_to_gate = tmp[["marker_to_gate"]]
    )


    tmpdir <- cytobench::local_tempdir_time()
    pdf(file.path(tmpdir, "removeme.pdf"), width = 20)
    print(tmp_cy)
    dev.off()
    # The events post gating should be n_events_postgate.
    # Therefore the total number of cells should be
    # number of samples * n_events_postgate
    testthat::expect_equal(
        nrow(tmp_cy$results_flowsom$cells_clusters_from_train),
        length(tmp$fs) * n_events_postgate
    )
    df_result1 <- data.frame(
        "File" = c(
            "StudyA_d01_s01_CyA",
            "StudyA_d02_s02_CyA",
            "StudyA_d03_s03_CyA",
            "StudyA_d01_s01_CyB",
            "StudyA_d02_s02_CyB",
            "StudyA_d03_s03_CyB"
        ),
        "OTD_to_referencesample" = c(
            0.002165620,
            0.002185727,
            0.002199899,
            0.002504403,
            0.002519673,
            0.002542104
        )
    )
    # The following is - even though we set the seed equally -
    # not equal because the gatingsets are not identical, therefore
    # NOT the exact same cells have been selected to calculate the OTD on.
    testthat::expect_false(
        all(
            abs(
                df_result1$OTD_to_referencesample -
                    tmp_cy$calculated_otd$OTD_to_referencesample
            ) < 1e-7
        )
    )
})
