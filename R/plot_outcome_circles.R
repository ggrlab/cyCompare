#' Plot Outcome Proportions by Supersample Groupings as Circular Bar Charts
#'
#' This function creates circular bar plots for one or more outcome variables,
#' grouped by combinations of columns (supersample groupings). Each plot shows
#' the distribution of outcome values across the grouping.
#'
#' @param df A `data.frame` containing both the grouping and outcome columns.
#' @param dfcol_grouping_supersamples A character vector specifying the column names
#'   to be used for creating supersample groups via interaction.
#' @param dfcol_outcomes A character vector specifying the column names of the
#'   outcome variables to be plotted.
#'
#' @return A named list of `ggplot` objects (one per outcome), each displaying
#' a radial bar chart of outcome frequencies grouped by supersample combination.
#'
#' @examples
#' \dontrun{
#' df <- data.frame(
#'     batch = rep(c("A", "B"), each = 50),
#'     device = rep(c("X", "Y"), times = 50),
#'     outcome1 = sample(c("pos", "neg"), 100, replace = TRUE),
#'     outcome2 = sample(c("high", "low"), 100, replace = TRUE)
#' )
#' plots <- plot_outcome_circles(df, c("batch", "device"), c("outcome1", "outcome2"))
#' plots[["outcome1"]]
#' }
#'
#' @export
plot_outcome_circles <- function(df,
                                 dfcol_grouping_supersamples,
                                 dfcol_outcomes) {
    col_supersamples <- count <- NULL
    # Create a new column encoding supersample groupings
    df[["col_supersamples"]] <- interaction(
        df[, dfcol_grouping_supersamples]
    )

    # Generate a list of circular bar plots (one per outcome)
    plotlist <- sapply(dfcol_outcomes, simplify = FALSE, function(outcome_x) {
        ggplot2::ggplot(
            df,
            ggplot2::aes(
                x = col_supersamples,
                fill = factor(!!ggplot2::sym(outcome_x))
            )
        ) +
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
                x = "Supersample grouping",
                y = "Count",
                fill = outcome_x
            )
    })

    return(plotlist)
}
