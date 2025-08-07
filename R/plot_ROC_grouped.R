#' Plot ROC Curves Grouped by Metadata for Multiple Models and Outcomes
#'
#' This function takes a structured ROC result object (as returned by your pipeline)
#' and generates grouped ROC plots for each model/outcome combination. ROC curves are grouped
#' by training configuration and colored using the specified device colors.
#'
#' @param calced_roc
#' A list with `results` and `groups` entries, where `results` contains ROC objects
#' and `groups` describes the training configuration for each result.
#' @param dfcol_grouping_supersamples
#' Character vector of column names used to group across supersamples (e.g., `"Study"`).
#' @param dfcol_grouping_samples
#' Character scalar naming the sample grouping column (e.g., `"Device"`).
#' @param dfcol_outcomes
#' Character vector of column names corresponding to different outcome types (e.g., `"outcome_1"`, `"outcome_2"`).
#' @param dfcol_train_validation_other
#' Character scalar for column indicating `"train"`, `"validation"`, or `"test"` set.
#' @param dv_class_positive
#' A named vector indicating the positive class for each outcome variable.
#' @param device_colors
#' A named character vector assigning colors to device names.
#'
#' @return A `data.frame` (tibble) where each row contains:
#'   - Grouping metadata (e.g., `"train.Device"`, `"train.Study"`, ...)
#'   - Outcome/model/clustering info
#'   - ROC plot objects (`roc.train`, `roc.validation`, `roc.test`, etc.) as `ggplot` objects.
#'
#' @export
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
    # Initialize result table
    roc_plots_df_total <- tibble::tibble()

    # Loop over all trained groups
    for (group_trained_on_i in seq_len(nrow(calced_roc[["groups"]]))) {
        res_current <- calced_roc[["results"]][[group_trained_on_i]]
        group_current <- calced_roc[["groups"]][group_trained_on_i, ]

        # Copy and rename group columns to indicate training context
        group_current_training <- group_current
        colnames(group_current_training) <- paste0("train.", colnames(group_current_training))

        # Group by outcome/model/clustering/Device/Study (or similar)
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
                    tmp <- .x[["proc"]] # Contains named list: train/test/validation ROC objects
                    names(tmp) <- .x[[dfcol_train_validation_other]]

                    # Reorder to always display in test/validation/train order
                    ordered_tmp <- list(
                        "test" = tmp[["test"]],
                        "validation" = tmp[["validation"]],
                        "train" = tmp[["train"]]
                    )
                    # Add any additional subsets not listed above
                    ordered_tmp <- c(
                        ordered_tmp,
                        tmp[setdiff(names(tmp), names(ordered_tmp))]
                    )
                    # Drop NULL entries
                    ordered_tmp <- ordered_tmp[!sapply(ordered_tmp, is.null)]

                    # Generate ggroc plots
                    ggrocs <- plot_ggroc_multiple(
                        ordered_tmp,
                        col_or_linetype = "linetype",
                        geom_line_args = list(
                            col = device_colors[group_current[[dfcol_grouping_samples[[1]]]]]
                        )
                    )

                    # Add each plot as new list-column entry
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

        # Combine group metadata and ROC plot objects
        rocs_df <- do.call(rbind, rocs_df)
        roc_plots_df_total <- rbind(
            roc_plots_df_total,
            cbind(group_current_training, rocs_df)
        )
    }

    return(roc_plots_df_total)
}
