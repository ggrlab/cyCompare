plot_samples_by_time <- function(df, dfcol_grouping_supersamples = NULL) {
    if (!all(is.null(dfcol_grouping_supersamples))) {
        df <- df |>
            dplyr::group_by(!!!rlang::syms(dfcol_grouping_supersamples))
    }
    df_cumulative <- df |>
        dplyr::arrange(Time) |>
        dplyr::mutate(
            cumulative_sID = cumsum(!duplicated(Sample)),
            cumulative_ssID = cumsum(!duplicated(SuperSample)),
            date_minus_first = Time - min(Time)
        ) |>
        dplyr::ungroup() |>
        tidyr::pivot_longer(
            cols = c("cumulative_sID", "cumulative_ssID"),
            values_to = "n",
            names_to = "cumulative_x"
        ) |>
        dplyr::mutate(
            `Number of` = ifelse(grepl("ssID", cumulative_x), "SuperSamples", "Samples")
        )
    p0 <-
        ggplot2::ggplot(
            df_cumulative,
            ggplot2::aes(x = Time, y = n, col = `Number of`, group = !!rlang::sym(dfcol_grouping_supersamples[[1]]))
        ) +
        ggplot2::geom_line(linewidth = .6) +
        ggplot2::geom_point(shape = 3) +
        ggpubr::theme_pubclean() +
        ggplot2::labs(title = "Donors and samples over time", x = "Date", y = "Number of donors") +
        ggplot2::scale_color_brewer(palette = "Set2") +
        ggplot2::theme(
            legend.text = ggplot2::element_text(size = 5),
            legend.title = ggplot2::element_text(size = 7)
        ) +
        ggplot2::guides(
            color = ggplot2::guide_legend(nrow = 2, byrow = TRUE),
            linetype = ggplot2::guide_legend(nrow = 3, byrow = TRUE)
        )
    return(p0)
}
