test_that("CyCompare Badalona samples", {
    tmp <- prepare_testdata()
    ff_list_downsampled <- tmp[["ff_list_downsampled"]]
    df <- tmp[["df"]]
    relevant_mn <- tmp[["relevant_mn"]]
    marker_to_gate <- tmp[["marker_to_gate"]]
    gatingsets <- tmp[["gatingsets"]]

    tmp <- cycompare(
        flowframes = ff_list_downsampled,
        ff_columns_relevant = names(relevant_mn),
        df = df,
        n_events_postgate = 1e4,
        gatingsets = gslist,
        gatename_primary = "/Singlets/CD45+/CD3+",
        marker_to_gate = marker_to_gate
    )
    pdf("removeme.pdf", width = 20)
    print(tmp)
    dev.off()
    testthat::expect_true(TRUE)
})
test_that("CyCompare Badalona samples, no gatingset", {
    tmp <- prepare_testdata()
    ff_list_downsampled <- tmp[["ff_list_downsampled"]]
    df <- tmp[["df"]]
    relevant_mn <- tmp[["relevant_mn"]]
    marker_to_gate <- tmp[["marker_to_gate"]]
    tmp <- cycompare(
        flowframes = ff_list_downsampled,
        ff_columns_relevant = names(relevant_mn),
        df = df,
        n_events_postgate = 1e4,
        gatingsets = NULL,
        # If the gatingsets are NULL; gatename_primary is replaced by "root"
        gatename_primary = "/Singlets/CD45+/CD3+",
        # If the gatingsets are NULL; All elements of marker_to_gate are replaced by "root"
        # marker_to_gate[TRUE] <- "root"
        marker_to_gate = marker_to_gate
    )
    pdf("removeme.pdf", width = 20)
    print(tmp)
    dev.off()
    testthat::expect_true(TRUE)
})
