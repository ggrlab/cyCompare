#' Plot multiple ROC curves with ggplot2
#'
#' This function generates a combined ROC plot from a named list of `pROC::roc` objects.
#' It allows customization by mapping the curve identity to either color or linetype and supports additional arguments passed to `geom_line`.
#'
#' The function returns a list containing:
#' - A full ROC plot with legend
#' - A ROC plot without the legend
#' - A standalone ggplot object for the legend
#'
#' @param named_proc_list
#' Named list of `pROC::roc` objects. Names will be used in the legend and appended with AUC and CI.
#' @param col_or_linetype
#' Character vector of length 1, either `"linetype"` or `"color"`.
#' Determines whether ROC curves are distinguished by line type or color.
#' @param geom_line_args
#' List of additional arguments passed to `ggplot2::geom_line()` (e.g., `list(size = 1)`).
#'
#' @return A named list with three ggplot2 objects:
#' \describe{
#'   \item{plot}{Full ROC plot with legend}
#'   \item{onlyROC}{ROC plot without the legend}
#'   \item{onlyLEGEND}{Legend as a standalone plot}
#' }
#'
#' @export
#' @examples
#' data(aSAH, package = "pROC")
#' r1 <- pROC::roc(aSAH$outcome, aSAH$s100b, levels = c("Good", "Poor"), direction = "<")
#' r2 <- pROC::roc(aSAH$outcome, -aSAH$s100b, levels = c("Good", "Poor"), direction = "<")
#' r3 <- pROC::roc(aSAH$outcome, aSAH$s100b + rnorm(nrow(aSAH), 0, .1), levels = c("Good", "Poor"), direction = "<")
#' named_proc_list <- list(
#'     "ROC_1" = r1,
#'     "ROC_2" = r2,
#'     "ROC_3" = r3
#' )
#'
#' print(plot_ggroc_multiple(
#'     named_proc_list,
#'     col_or_linetype = "color", geom_line_args = list(size = 1)
#' ))
#' print(plot_ggroc_multiple(named_proc_list))
#' print(plot_ggroc_multiple(
#'     named_proc_list,
#'     col_or_linetype = "linetype", geom_line_args = list(size = 1)
#' ))
plot_ggroc_multiple <- function(named_proc_list, col_or_linetype = c("linetype", "color"), geom_line_args = NULL) {
    # Append AUC and confidence intervals to each curve's name
    names(named_proc_list) <- paste0(
        names(named_proc_list),
        ". AUC=",
        sapply(named_proc_list, function(x) {
            paste0(
                round(pROC::auc(x), 2),
                " (",
                round(pROC::ci.auc(x)[1], 2),
                ", ",
                round(pROC::ci.auc(x)[3], 2),
                ")"
            )
        })
    )

    # Convert ROC list to data frame for ggplot
    proc_data <- pROC::ggroc(named_proc_list)

    # Set up ggplot mapping aesthetics based on user preference
    if (col_or_linetype[1] == "color") {
        joint_proc <- ggplot2::ggplot(
            proc_data$data,
            ggplot2::aes(x = `1-specificity`, y = sensitivity, color = name)
        )
    } else {
        joint_proc <- ggplot2::ggplot(
            proc_data$data,
            ggplot2::aes(x = `1-specificity`, y = sensitivity, linetype = name)
        )
    }

    # Add ROC lines and styling
    if (!all(is.null(geom_line_args))) {
        joint_proc <- joint_proc +
            do.call(ggplot2::geom_line, geom_line_args) # Add geom_line with user-supplied arguments
    }
    joint_proc <- joint_proc +
        ggplot2::labs(
            x = "False Positive Rate",
            y = "True Positive Rate"
        ) +
        ggpubr::theme_pubr() + # Use publication-ready theme
        ggplot2::coord_fixed(ratio = 1) + # Fix aspect ratio
        ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dotted") + # Diagonal line
        ggplot2::theme(
            title = ggplot2::element_text(size = 6)
        )

    # Return full plot, plot without legend, and standalone legend
    return(
        list(
            "plot" = joint_proc,
            "onlyROC" = joint_proc + ggplot2::theme(legend.position = "none"),
            "onlyLEGEND" = ggpubr::as_ggplot(ggpubr::get_legend(joint_proc +
                ggplot2::guides(
                    color = ggplot2::guide_legend(ncol = 1)
                )))
        )
    )
}
