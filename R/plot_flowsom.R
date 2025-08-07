#' Perform FlowSOM Clustering and Generate Device-Comparative Plots
#'
#' This function clusters gated flow cytometry data using FlowSOM and generates PCA and MA plots
#' to compare cluster distributions across devices. Optionally, only a training subset of the data
#' can be used for clustering.
#'
#' @param ff_gated A list of `flowFrame` objects containing gated flow cytometry data.
#' @param df
#' A `data.table` with metadata, containing at least the columns `"File"`,
#' `dfcol_grouping_samples`, and `"Sample"`.
#' @param device_colors A named vector assigning a color to each device. Used in plots.
#' @param transformlist
#' A list of transformation functions per marker. If a single function is given, it is applied to all markers.
#' @param nClus Integer specifying the number of clusters for FlowSOM.
#' @param scale Logical. If `TRUE`, scales the data before clustering.
#' @param xdim Integer. X-dimension of the FlowSOM grid.
#' @param ydim Integer. Y-dimension of the FlowSOM grid.
#' @param seed Integer. Random seed for reproducibility.
#' @param dfcol_grouping_samples
#' Character scalar specifying the column in `df` used to group samples in PCA plots (default: `"Device"`).
#' @param dfcol_train_validation_other
#' Optional character scalar. If set, this column is used to select the training subset
#' (`"train"`) for clustering and to shape/facet the plots.
#' @return A named list with two elements:
#'   \item{plots_pca}{A list of PCA plots showing clustering results colored by device.}
#'   \item{p_MA}{A list of MA plots comparing cluster proportions across devices.}
#'
#' @export
#'
#' @examples
#' # Simulate flow cytometry data
#' fs <- cytobench::simulate_fs(
#'     n_samples = 4,
#'     ncells = 250,
#'     columns = c("CD4", "CD8", "CD3")
#' )
#'
#' # Create metadata for each sample
#' df_meta <- data.table::data.table(
#'     File = flowCore::sampleNames(fs),
#'     Device = c("A", "A", "B", "B"),
#'     Sample = c("S1", "S2", "S1", "S2"),
#'     SuperSample = "Study_Z",
#'     TrainTest = c("train", "test", "train", "test")
#' )
#'
#' # Define color palette for devices
#' device_cols <- c("A" = "steelblue", "B" = "firebrick")
#'
#' # Run FlowSOM clustering and generate plots
#' result <- plot_flowsom(
#'     ff_gated = fs,
#'     df = df_meta,
#'     device_colors = device_cols,
#'     transformlist = function(x) asinh(x / 1000),
#'     nClus = 4,
#'     dfcol_grouping_samples = "Device",
#'     dfcol_train_validation_other = "TrainTest"
#' )
#'
#' # Show first PCA plot
#' print(result$plots_pca)
#'
#' # Show first MA plot (first device combination, first clustering)
#' print(result$p_MA)
plot_flowsom <- function(ff_gated,
                         df,
                         device_colors = NULL,
                         transformlist,
                         nClus = 5,
                         scale = FALSE,
                         xdim = 3,
                         ydim = 3,
                         seed = 3711283,
                         dfcol_grouping_samples = "Device",
                         dfcol_train_validation_other = NULL,
                         MA_horizontal_lines_FC = c(2, 10, 25),
                         MA_bins = 100) {
    File <- SuperSample <- Sample <- cluster_id <- proportion <- count <- NULL # R CMD check compatibility
    # Convert list of flowFrames into a flowSet for transformation
    if ("flowSet" %in% class(ff_gated)) {
        gated_fs <- ff_gated
    } else {
        gated_fs <- flowCore::flowSet(ff_gated)
    }

    # Apply the transformation: use same function for all markers if a single function is provided
    fc_transformlist <- flowCore::transformList(
        flowCore::colnames(gated_fs[[1]]),
        transformlist_named(transformlist, flowCore::colnames(gated_fs[[1]]))
    )

    # Transform the flowSet using the provided transformlist
    gated_fs_transformed <- flowCore::transform(gated_fs, fc_transformlist)

    # Subset to training samples if specified
    if (!all(is.null(dfcol_train_validation_other))) {
        gated_fs_transformed_train <- gated_fs_transformed[
            dplyr::filter(df, !!rlang::sym(dfcol_train_validation_other) == "train") |>
                dplyr::pull(File)
        ]
    } else {
        gated_fs_transformed_train <- gated_fs_transformed
    }

    # Run FlowSOM clustering on the transformed training data
    flowsom_all <- FlowSOM::FlowSOM(
        input = gated_fs_transformed_train,
        transform = FALSE, # Already transformed above
        transformList = NULL,
        nClus = nClus,
        scale = scale,
        xdim = xdim,
        ydim = ydim,
        seed = seed
    )

    # Predict cluster assignments on all transformed data
    fs_pred <- cytobench::flowSOM_predict(flowsom_all, gated_fs_transformed)

    # --- PCA Plots ---
    plots_pca <- lapply(names(fs_pred[["ncells_per_x"]]), function(x) {
        x_data <- fs_pred[["ncells_per_x"]][[x]]

        # Center & scale data (except for sample column)
        x_data_scaled <- scale(dplyr::select(x_data, -sample), center = TRUE, scale = TRUE)

        # Remove columns with constant values (zero variance)
        x_data_scaled_noconstant.cols <- x_data_scaled[, !is.na(apply(x_data_scaled, 2, var))]

        # PCA on normalized cluster proportions
        res_pca <- stats::prcomp(x_data_scaled_noconstant.cols, scale = FALSE)

        # Generate PCA plot using ggfortify
        p0 <- ggfortify:::autoplot.prcomp(
            res_pca,
            data = df,
            colour = dfcol_grouping_samples[[1]],
            shape = dfcol_train_validation_other[[1]]
        ) +
            ggpubr::theme_pubr() +
            ggplot2::theme(legend.position = "top") +
            ggplot2::labs(
                title = "PCA of FlowSOM",
                subtitle = x
            )

        # Apply manual color palette if given
        if (!all(is.null(device_colors))) {
            p0 <- p0 + ggplot2::scale_color_manual(values = device_colors)
        }

        return(p0)
    })

    # --- MA Plots ---
    p_MA <- lapply(names(fs_pred[["ncells_per_x"]]), function(x) {
        x_data <- fs_pred[["ncells_per_x"]][[x]]
        x_data_numeric <- dplyr::select(x_data, -sample)

        # Compute within-sample proportions
        proportions <- sweep(x_data_numeric, 1, rowSums(x_data_numeric), "/")
        x_data[, -1] <- proportions

        # Join with metadata
        x_data_df <- dplyr::left_join(x_data, df, by = c("sample" = "File"))

        # Pivot to long format
        x_data_df_long <- tidyr::pivot_longer(
            x_data_df,
            cols = tidyr::all_of(grep(colnames(x_data_df), pattern = "[cC]luster", value = TRUE)),
            names_to = "cluster_id",
            values_to = "proportion"
        )

        # Wide format for per-device comparison
        x_data_df_persample <- x_data_df_long |>
            dplyr::select(
                !!rlang::sym(dfcol_grouping_samples[[1]]),
                !!!rlang::syms(dfcol_train_validation_other),
                SuperSample,
                Sample,
                cluster_id,
                proportion
            ) |>
            tidyr::pivot_wider(
                names_from = dfcol_grouping_samples[[1]],
                values_from = "proportion"
            )

        # All pairwise device comparisons
        device_combinations <- combn(names(device_colors), 2)
        plotlist <- list()

        for (combination_i in seq_len(ncol(device_combinations))) {
            device_combination <- device_combinations[, combination_i]
            d1 <- device_combination[1]
            d2 <- device_combination[2]

            # Generate MA plot (log2 fold-change vs average)
            p0 <- ggplot2::ggplot(
                x_data_df_persample,
                ggplot2::aes(
                    y = log2(!!ggplot2::sym(d1) / !!ggplot2::sym(d2)),
                    x = ((!!ggplot2::sym(d1) + !!ggplot2::sym(d2)) / 2)
                )
            ) +
                ggplot2::ggtitle(paste0(d1, " vs ", d2), subtitle = x) +
                ggplot2::geom_abline(intercept = 0, slope = 0) +
                ggplot2::theme_bw() +
                ggplot2::scale_x_log10() +
                ggplot2::theme(
                    axis.text.x = ggplot2::element_text(size = 6),
                    legend.position = "right",
                    legend.key.size = ggplot2::unit(.5, "cm"),
                    legend.title = ggplot2::element_text(angle = -90)
                )
            minimum_x_value <- min(
                x_data_df_persample[[d1]],
                x_data_df_persample[[d2]],
                na.rm = TRUE
            )
            # Facet by train/test if applicable
            if (!all(is.null(dfcol_train_validation_other))) {
                p0 <- p0 + ggplot2::facet_wrap(dfcol_train_validation_other)
            }

            # Add horizontal reference lines for fold-changes (e.g., x2, x10, x25)
            for (specific_lines in MA_horizontal_lines_FC) {
                p0 <- p0 +
                    ggplot2::geom_hline(
                        yintercept = c(-log2(specific_lines), log2(specific_lines)),
                        linetype = "dashed"
                    ) +
                    ggplot2::annotate(
                        "text",
                        x = minimum_x_value, y = log2(specific_lines),
                        label = paste0("x", specific_lines),
                        hjust = 0, vjust = 0, size = 5 # 2.5
                    )
            }

            # Add density information via 2D binning
            plotlist[[combination_i]] <- p0 +
                ggplot2::stat_bin_2d(
                    ggplot2::aes(fill = log10(ggplot2::after_stat(count) + 1)),
                    bins = MA_bins,
                    geom = "tile"
                ) +
                ggplot2::labs(fill = "Log10(n cells +1)") +
                ggplot2::scale_fill_viridis_c(na.value = NA)
        }

        return(plotlist)
    })

    # Return plots
    return(list(
        plots_pca = plots_pca,
        p_MA = p_MA
    ))
}
