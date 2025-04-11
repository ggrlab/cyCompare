plot_ROC_grouped <- function(calced_roc,
                             dfcol_grouping_supersamples = c("Study"),
                             dfcol_grouping_samples = "Device",
                             dfcol_outcomes = c("outcome_1", "outcome_2"),
                             dfcol_train_validation_other = "train_validation_test",
                             dv_class_positive = c("outcome_1" = "A", "outcome_2" = 5.1),
                             device_colors = c(
                                 "Fortessa" = "#FF0000",
                                 "Cytoflex" = "#00FF00",
                                 "Aurora" = "#0000FF"
                             )) {
    roc_plots_df_total <- tibble::tibble()
    for (group_trained_on_i in seq_len(nrow(calced_roc[["groups"]]))) {
        res_current <- calced_roc[["results"]][[group_trained_on_i]]
        group_current <- calced_roc[["groups"]][group_trained_on_i, ]
        group_current_training <- group_current
        colnames(group_current_training) <- paste0("train.", colnames(group_current_training))
        rocs_df <- res_current |>
            dplyr::group_by(
                !!!rlang::syms(
                    c(
                        "outcome_",
                        "model_",
                        "clustering_",
                        dfcol_grouping_supersamples,
                        dfcol_grouping_samples
                    )
                )
            ) |>
            dplyr::group_map(
                ~ {
                    tmp <- .x[["proc"]]
                    names(tmp) <- .x[[dfcol_train_validation_other]]
                    ordered_tmp <- list(
                        "test" = tmp[["test"]],
                        "validation" = tmp[["validation"]],
                        "train" = tmp[["train"]]
                    )
                    ordered_tmp <- c(
                        ordered_tmp,
                        tmp[setdiff(names(tmp), names(ordered_tmp))]
                    )
                    # Remove NULL values
                    ordered_tmp <- ordered_tmp[!sapply(ordered_tmp, is.null)]
                    ggrocs <- plot_ggroc_multiple(
                        ordered_tmp,
                        col_or_linetype = "linetype",
                        geom_line_args = list(col = device_colors[group_current[[dfcol_grouping_samples[[1]]]]])
                    )
                    for (name_x in names(ggrocs)) {
                        .y[[paste0("roc.", name_x)]] <- list(
                            ggrocs[[name_x]] +
                                ggplot2::ggtitle(
                                    "Trained on",
                                    subtitle = paste0(
                                        names(group_current), ": ",
                                        group_current,
                                        collapse = ", "
                                    )
                                )
                        )
                    }
                    return(.y)
                }
            )
        rocs_df <- do.call(rbind, rocs_df)
        roc_plots_df_total <- rbind(
            roc_plots_df_total,
            cbind(group_current_training, rocs_df)
        )
    }
    return(roc_plots_df_total)
}
