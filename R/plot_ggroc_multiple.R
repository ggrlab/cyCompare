plot_ggroc_multiple <- function(named_proc_list, col_or_linetype = c("linetype", "color"), geom_line_args = NULL) {
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
    
    proc_data <- pROC::ggroc(named_proc_list)
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
    joint_proc <- joint_proc +
        do.call(ggplot2::geom_line, geom_line_args) +
        ggplot2::labs(
            x = "False Positive Rate",
            y = "True Positive Rate"
        ) +
        ggpubr::theme_pubr() +
        ggplot2::coord_fixed(ratio = 1) +
        ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
        ggplot2::theme(
            title = ggplot2::element_text(size = 6)
        )

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
