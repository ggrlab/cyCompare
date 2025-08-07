#' Plot PCA of Cluster Proportions (FlowSOM)
#'
#' @param fs_pred Output from `run_flowsom_clustering()`.
#' @param df Metadata table.
#' @param device_colors Named color vector.
#' @param dfcol_grouping_samples Column used to color samples.
#' @param dfcol_train_validation_other Optional column used for shape/facets.
#'
#' @return A list of ggplot2 PCA plots.
#' @export
plot_flowsom_pca <- function(fs_pred,
                             df,
                             device_colors = NULL,
                             dfcol_grouping_samples = "Device",
                             dfcol_train_validation_other = NULL) {
    lapply(names(fs_pred[["ncells_per_x"]]), function(x) {
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
}
