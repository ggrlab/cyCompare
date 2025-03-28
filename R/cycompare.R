#' Compare Flow Cytometry Data Across Devices
#'
#' This function performs a comparative analysis of flow cytometry data across multiple devices.
#' It includes basic sample statistics, gating, density plots, FlowSOM clustering, and marker intensity comparisons.
#'
#' @param flowframes A named list of `flowFrame` objects containing flow cytometry data.
#' @param df A `data.table` containing metadata with at least the columns `"File"`, `"Device"`, and `"Sample"`.
#' @param ff_columns_relevant A character vector specifying the relevant markers for analysis.
#' @param transformlist A transformation function or a named list of functions for transforming marker intensities
#'        (default: `asinh(x / 1e3)`).
#' @param gatingsets A named list of gating sets for each dataset.
#' @param gatename_primary A character string specifying the primary gating population.
#' @param max_events_postgate An integer specifying the maximum number of events to keep post-gating.
#' @param marker_to_gate A named vector mapping marker names to their corresponding gates.
#' @param device_colors A named vector or function that provides colors for each device.
#'        If a function is provided, it should take the number of devices as input and return a vector of colors.
#' @param nClus An integer specifying the number of clusters for FlowSOM clustering (default: 5).
#' @param scale A logical indicating whether to scale the data in FlowSOM clustering (default: `FALSE`).
#' @param xdim An integer specifying the x-dimension of the FlowSOM grid (default: 3).
#' @param ydim An integer specifying the y-dimension of the FlowSOM grid (default: 3).
#' @param seed An integer specifying the random seed for FlowSOM clustering (default: `3711283`).
#' @param ... Additional parameters passed to `plot_flowsom()`.
#'
#' @return A named list of ggplot2 objects containing:
#'   \item{"Samples over time per device"}{A plot showing the number of samples collected over time per device.}
#'   \item{"Counts and percentages"}{Plots of gated cell counts and percentages per sample.}
#'   \item{"Positive population MFI"}{Plots showing the median fluorescence intensity (MFI) of positive gated populations.}
#'   \item{"Density plots"}{Density distributions of marker intensities across devices and samples.}
#'   \item{"Flowsom_PCA"}{Principal Component Analysis (PCA) plots of FlowSOM clustering results.}
#'   \item{"Flowsom_MA"}{MA plots comparing cluster proportions between devices.}
#'
#' @export

cycompare <- function(
    flowframes,
    df,
    ff_columns_relevant,
    transformlist = function(x) asinh(x / 1e3),
    gatingsets,
    gatename_primary,
    max_events_postgate = 10e3,
    marker_to_gate,
    device_colors = function(n) {
        RColorBrewer::brewer.pal(n, "Dark2")
    },
    # OTD parameters
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
    nClus = 5,
    scale = FALSE,
    xdim = 3,
    ydim = 3,
    seed = 3711283) {
    unique_devices <- unique(df[["Device"]])
    if (!all(unique_devices %in% names(device_colors))) {
        if (is.function(device_colors)) {
            # Automatically assign colors
            device_colors <- setNames(
                # RColorBrewer takes a minimum of 3 colors, thus the seq_along
                device_colors(length(unique_devices))[seq_along(unique_devices)],
                unique_devices
            )
            warning("The following device colors were automatically assigned:")
            warning(paste0(device_colors, collapse = ", "))
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

    if (quantile(counts_joint[pop == gatename_primary][["count"]], .9) < 100) {
        stop(
            "The primary gate has less than 100 cells in the 90th percentile of samples. Did you select the right gate for these samples?"
        )
    }

    # Plots of counts and percentages
    p1.2_3 <- plot_counts(counts_joint, gatename_primary, device_colors)

    #### Figure 2 plots
    # 1. Positive population MFI
    p2.1 <- plot_MFI_positivegates(
        dt_count_mfi = counts_joint,
        marker_to_gate = marker_to_gate,
        device_colors = device_colors,
        transformlist = transformlist
    )
    p2.1_ratio <- plot_MFI_positivegates(
        dt_count_mfi = counts_joint,
        marker_to_gate = marker_to_gate,
        device_colors = device_colors,
        transformlist = transformlist,
        meanratio = TRUE
    )
    gated_ff <- lapply(gated_ff, function(x) x[["flowset_gated"]][[1]])
    # Limit the number of events post-gating and select only the relevant columns
    gated_ff <- lapply(gated_ff, function(x) {
        nevents <- min(flowCore::nrow(x), max_events_postgate)
        x[sample(flowCore::nrow(x), nevents, replace = FALSE), ff_columns_relevant]
    })

    # 2. Density plots
    p2.2 <- plot_densities(
        ff_gated = gated_ff,
        df = df,
        device_colors = device_colors,
        transformlist = transformlist
    )

    # 3. OTD
    p3.1 <- plot_otd(
        ff_gated = gated_ff,
        df = df,
        device_colors = device_colors,
        transformlist = transformlist,
        n_mastersample = max_events_postgate,
        kwargs_loss = OTD_kwargs_loss
    )

    # 4 Clustering with FlowSOM
    p_flowsom <- plot_flowsom(
        ff_gated = gated_ff,
        df = df,
        device_colors = device_colors,
        transformlist = transformlist,
        nClus = nClus,
        scale = scale,
        xdim = xdim,
        ydim = ydim,
        seed = seed
    )

    # 5. Outcome prediction

    return(
        list(
            "Samples over time per device" = p1.1,
            "Counts and percentages" = p1.2_3,
            "Positive population MFI" = p2.1,
            "Positive population MFI ratio" = p2.1_ratio,
            "Density plots" = p2.2,
            "OTD to mastersample" = p3.1,
            "Flowsom_PCA" = p_flowsom[["plots_pca"]],
            "Flowsom_MA" = p_flowsom[["p_MA"]]
        )
    )
}
