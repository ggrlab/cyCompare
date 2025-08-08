test_that("CyCompare, outcomes", {
    tempdir <- cytobench::local_tempdir_time()
    prepdata <- prepare_testdata_fcs(tempdir = tempdir)

    # pacman::p_load("mlr3learners")

    tempdir2 <- cytobench::local_tempdir_time()
    result_analysis <- cycompare_outcomes_analyse(
        flowframes = flowCore::flowSet_to_list(prepdata$fs),
        ff_columns_relevant = names(prepdata$relevant_mn),
        df = prepdata$df,
        dfcol_grouping_supersamples = c("Study"),
        dfcol_grouping_samples = "Device",
        dfcol_outcomes = c("outcome_1", "outcome_2"),
        outcome_models = list("glmnet" = glmnet::cv.glmnet),
        n_events_postgate = 1e4,
        gatingsets = prepdata[["gatingsets"]],
        gatename_primary = "/FITC_PE_gate",
        marker_to_gate = prepdata$marker_to_gate,
        dfcol_train_validation_other = "train_validation_test",
        outdir_base = tempdir2
    )
    generated_filelist <- list.files(tempdir2, recursive = TRUE)
    expected_filelist <- c(
        "clustering/Study.EX___Device.CyA___train_validation_test.train/p1-FlowSOM.pdf",
        "clustering/Study.EX___Device.CyA___train_validation_test.train/r1-FlowSOM_result_train.qs",
        "clustering/Study.EX___Device.CyA___train_validation_test.train/r2-FlowSOM_result_train_predictions.qs",
        "clustering/Study.EX___Device.CyB___train_validation_test.train/p1-FlowSOM.pdf",
        "clustering/Study.EX___Device.CyB___train_validation_test.train/r1-FlowSOM_result_train.qs",
        "clustering/Study.EX___Device.CyB___train_validation_test.train/r2-FlowSOM_result_train_predictions.qs",
        "clustering_applied/Study.EX___Device.CyA___train_validation_test.train.qs",
        "clustering_applied/Study.EX___Device.CyB___train_validation_test.train.qs",
        "clustering_applied/ncells_per_cluster.csv",
        "clustering_applied/ncells_per_metaCluster.csv",
        "clustering_applied/proportions_per_cluster.csv",
        "clustering_applied/proportions_per_metaCluster.csv",
        "models/Study.EX___Device.CyA___train_validation_test.train/final_models.qs",
        "models/Study.EX___Device.CyA___train_validation_test.train/learners_tuned.qs",
        "models/Study.EX___Device.CyA___train_validation_test.train/predictions.csv",
        "models/Study.EX___Device.CyB___train_validation_test.train/final_models.qs",
        "models/Study.EX___Device.CyB___train_validation_test.train/learners_tuned.qs",
        "models/Study.EX___Device.CyB___train_validation_test.train/predictions.csv"
    )
    testthat::expect_equal(
        sort(generated_filelist),
        sort(expected_filelist)
    )

    qs::qsave(result_analysis, file.path(tempdir2, "result_analysis.qs"))
    result_analysis <- qs::qread(file.path(tempdir2, "result_analysis.qs"))
    plotted_results <- cycompare_outcomes_plot(
        df = prepdata$df,
        dfcol_grouping_supersamples = c("Study"),
        dfcol_grouping_samples = "Device",
        dfcol_outcomes = c("outcome_1", "outcome_2"),
        results_cycompare_analyse = result_analysis
    )
    pdf(file.path(tempdir2, "removeme.pdf"), width = 20)
    print(plotted_results)
    dev.off()
})
