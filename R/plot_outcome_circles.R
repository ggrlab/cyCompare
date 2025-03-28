plot_outcome_circles <- function(df,
                                 dfcol_grouping_supersamples,
                                 dfcol_outcomes) {
    df[["col_supersamples"]] <- interaction(
        # I tested the following, works as intended with multiple columns
        df[, dfcol_grouping_supersamples]
    )
    plotlist <- sapply(dfcol_outcomes, simplify = FALSE, function(outcome_x) {
        ggplot2::ggplot(df, ggplot2::aes(x = col_supersamples, fill = factor(!!ggplot2::sym(outcome_x)))) +
            ggplot2::geom_bar(width = 1) +
            ggplot2::coord_radial(theta = "y", expand = TRUE) +
            ggpubr::theme_pubr() +
            ggplot2::stat_count(
                geom = "text",
                ggplot2::aes(label = ggplot2::after_stat(count)),
                position = ggplot2::position_stack(vjust = 0.5),
                color = "white"
            ) +
            ggplot2::labs(
                # title = paste0("Outcome: ", outcome_x),
                x = "Supersample grouping",
                y = "Count",
                fill = outcome_x
            )
    })
    return(plotlist)
}
