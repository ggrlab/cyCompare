#' Compare Flow Cytometry Data Across Devices
#'
#' Performs a standardized cross-device comparison pipeline for flow cytometry data.
#' This includes gating, marker transformation, FlowSOM clustering, density and MFI plots,
#' and optionally optimal transport distance (OTD) analysis. Can save or load preprocessed data for efficiency.
#'
#' @inheritParams cycompare_outcomes_analyse
#' @param postgate_sample_seed Integer seed used for reproducible post-gating subsampling (default: 42).
#' @param prepared_saveload
#' Optional character path. If the file exists, will load preprocessed data with `qs::qread()`;
#' if not, will save it after processing using `qs::qsave()`.
#' @param nClus Number of clusters for FlowSOM (default: 5).
#' @param scale Logical, whether to scale markers for FlowSOM (default: `FALSE`).
#' @param xdim An integer specifying the x-dimension of the FlowSOM grid (default: 3).
#' @param ydim An integer specifying the y-dimension of the FlowSOM grid (default: 3).
#' @param do_otd Logical. Whether to perform optimal transport distance (OTD) comparison (default: `TRUE`).
#' @param OTD_kwargs_loss Named list of arguments for the OTD loss function.
#' @param do_flowsom Logical. Whether to perform FlowSOM clustering and plot comparisons (default: `TRUE`).
#' @param flowsom_seed An integer specifying the random seed for FlowSOM clustering (default: `3711283`).
#'
#' @return A named list of `ggplot2` plots:
#' \describe{
#'   \item{Samples over time per device}{Number of samples per device over time.}
#'   \item{Counts and percentages}{Per-sample cell counts and gated proportions.}
#'   \item{Positive population MFI}{MFI values of gated markers across devices.}
#'   \item{Positive population MFI ratio}{Ratio of MFI values across devices.}
#'   \item{Density plots}{Marker intensity distributions for each device.}
#'   \item{OTD to mastersample}{Pairwise distances between devices or conditions.}
#'   \item{Flowsom_PCA}{Cluster-level PCA across samples.}
#'   \item{Flowsom_MA}{Mean-vs-difference cluster abundance plots.}
#' }
#'
#' @export
#'
cycompare <- function(
    flowframes,
    df,
    ff_columns_relevant,
    dfcol_grouping_supersamples = NULL,
    dfcol_grouping_samples = "Device",
    dfcol_train_validation_other = "train_validation_test",
    transformlist = function(x) asinh(x / 1e3),
    gatingsets,
    gatename_primary,
    n_events_postgate = 10e3,
    postgate_sample_seed = 42,
    marker_to_gate,
    prepared_saveload = FALSE,
    device_colors = function(n) {
        RColorBrewer::brewer.pal(n, "Dark2")
    },
    # OTD parameters
    do_otd = TRUE,
    OTD_kwargs_loss = list(
        loss = lossfun_hist,
        verbose = FALSE,
        write_intermediate = FALSE,
        # Do not skip any comparisons
        should_skip = function(i, j) FALSE,
        take_time = FALSE,
        return_as_matrix = TRUE
    ),
    # FlowSOM parameters
    do_flowsom = TRUE,
    nClus = 5,
    scale = FALSE,
    xdim = 3,
    ydim = 3,
    flowsom_seed = 3711283) {
    # --- Load or prepare data ---
    if (is.character(prepared_saveload) && file.exists(prepared_saveload)) {
        message("Loading prepared data from ", prepared_saveload)
        prepared <- qs::qread(prepared_saveload)
    } else {
        prepared <- cycompare_preparation(
            flowframes = flowframes,
            df = df,
            ff_columns_relevant = ff_columns_relevant,
            device_colors = device_colors,
            gatingsets = gatingsets,
            gatename_primary = gatename_primary,
            n_events_postgate = n_events_postgate,
            seed = postgate_sample_seed,
            transformlist = transformlist,
            dfcol_grouping_samples = dfcol_grouping_samples,
            dfcol_train_validation_other = dfcol_train_validation_other,
            dfcol_grouping_supersamples = dfcol_grouping_supersamples,
            marker_to_gate = marker_to_gate
        )
        if (is.character(prepared_saveload)) {
            qs::qsave(prepared, prepared_saveload)
        }
    }

    # Extract components from prepared list
    gated_ff <- prepared[["gated_ff"]]
    counts_joint <- prepared[["counts_joint"]]
    device_colors <- prepared[["device_colors"]]
    marker_to_gate <- prepared[["marker_to_gate"]]
    gatename_primary <- prepared[["gatename_primary"]]

    # --- 1. Basic plots ---
    ## 1.1: Samples over time per device
    p1.1 <- plot_samples_by_time(
        df,
        dfcol_grouping_supersamples = dfcol_grouping_supersamples
    )

    ## 1.2 & 1.3: Gated cell counts and percentages
    p1.2_3 <- plot_counts(
        dt_counts = counts_joint,
        device_colors = device_colors
    )

    # --- 2. MFI and density plots ---
    ## 2.1: Median fluorescence intensity for positive populations
    p2.1 <- plot_MFI_positivegates(
        dt_count_mfi = counts_joint,
        marker_to_gate = marker_to_gate,
        device_colors = device_colors,
        transformlist = transformlist
    )

    ## 2.1: Ratio version of above
    p2.1_ratio <- plot_MFI_positivegates(
        dt_count_mfi = counts_joint,
        marker_to_gate = marker_to_gate,
        device_colors = device_colors,
        transformlist = transformlist,
        meanratio = TRUE
    )

    ## 2.2: Density distributions per marker/device
    p2.2 <- plot_densities(
        ff_gated = gated_ff,
        df = df,
        device_colors = device_colors,
        transformlist = transformlist,
        relevant_columns = ff_columns_relevant
    )

    # --- 3. Optimal Transport Distance (OTD) ---
    if (do_otd) {
        calculated_otd <- otd_to_reference(
            ff_gated = gated_ff,
            transformlist = transformlist,
            n_referencesample = n_events_postgate,
            kwargs_loss = OTD_kwargs_loss,
            relevant_columns = ff_columns_relevant
        )
        p3.1 <- plot_otd(
            df = df,
            otd_df = calculated_otd,
            device_colors = device_colors
        )
    } else {
        calculated_otd <- NULL
        p3.1 <- NULL
    }

    # --- 4. FlowSOM clustering and plots ---
    if (do_flowsom) {
        res_fs <- flowsom_run_train(
            ff_gated = gated_ff,
            relevant_columns = ff_columns_relevant,
            transformlist = transformlist,
            nClus = nClus,
            scale = scale,
            xdim = xdim,
            ydim = ydim,
            seed = flowsom_seed
        )
        p4.1 <- plot_flowsom_pca(
            fs_pred = res_fs,
            df = df,
            device_colors = device_colors,
            dfcol_grouping_samples = dfcol_grouping_samples,
            dfcol_train_validation_other = dfcol_train_validation_other
        )
        p4.2 <- plot_flowsom_ma(
            fs_pred = res_fs,
            df = df,
            device_colors = device_colors,
            dfcol_grouping_samples = dfcol_grouping_samples,
            dfcol_train_validation_other = dfcol_train_validation_other,
            MA_horizontal_lines_FC = c(2, 10, 25),
            MA_bins = 100
        )
    } else {
        res_fs <- NULL
        p4.1 <- NULL
        p4.2 <- NULL
    }

    # --- Return all collected plots ---
    return(
        list(
            "Samples over time per device" = p1.1,
            "Counts and percentages" = p1.2_3,
            "Positive population MFI" = p2.1,
            "Positive population MFI ratio" = p2.1_ratio,
            "Density plots" = p2.2,
            "OTD to mastersample" = p3.1,
            "Flowsom_PCA" = p4.1,
            "Flowsom_MA" = p4.2,
            "results_flowsom" = res_fs,
            "calculated_otd" = calculated_otd
        )
    )
}
