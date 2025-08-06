#' Extract ROC Coordinates and Classification Metrics
#'
#' Computes detailed performance metrics from a `pROC::roc` object at a specified threshold, including AUC,
#' confidence intervals, and classification metrics such as accuracy, kappa, and log loss.
#'
#' @param proc_result A `pROC::roc` object resulting from a call to `pROC::roc()`.
#' @param threshold A numeric value specifying the classification threshold. If `NA` (default),
#'   the optimal threshold based on the "closest.topleft" method is selected.
#'
#' @return A tibble with one row containing:
#' \itemize{
#'   \item Threshold-based ROC metrics (sensitivity, specificity, etc.)
#'   \item AUC and its confidence interval
#'   \item The original `pROC::roc` object (in list-column `proc`)
#'   \item Classification metrics derived via `caret::confusionMatrix()`
#'   \item Log loss from `mlr3measures::logloss()`
#'   \item Class frequencies
#' }
#'
#' @export
#'
#' @examples
#' data(aSAH, package = "pROC")
#' r <- pROC::roc(aSAH$outcome, aSAH$s100b,
#'     levels = c("Good", "Poor")
#' )
#' coords_helper(proc_result = r)
#' coords_helper(proc_result = r, threshold = 0.5)
coords_helper <- function(proc_result, threshold = NA) {
    # Choose threshold: use optimal by default, or specified threshold
    if (is.na(threshold)) {
        coords <- pROC::coords(
            proc_result,
            "best",
            ret = "all",
            best.method = "closest.topleft"
        )
    } else {
        coords <- pROC::coords(
            proc_result,
            x = threshold,
            input = "threshold",
            ret = "all"
        )
    }

    coords <- tibble::as_tibble(coords)
    coords <- coords[1, ] # In case multiple rows returned (e.g., multiple "best" points)

    # Add ROC object and AUC values
    coords[["proc"]] <- list(proc_result)
    coords[["auc"]] <- pROC::auc(proc_result)[[1]]

    # Compute AUC confidence interval
    # Suppressed warning for perfect AUC=1 case
    auc_ci <- suppressWarnings(as.numeric(pROC::ci.auc(proc_result)))
    coords[["auc_ci.low"]] <- auc_ci[1]
    coords[["auc_ci.high"]] <- auc_ci[3]

    # Reorder for readability
    coords <- coords |> dplyr::relocate(proc, auc, auc_ci.low, auc_ci.high)

    # Now add caret results
    # confusionMatrix(pred, truth)
    #  “>” (default for multivariate curves): if the predictor values for the control group are higher
    # than the values of the case group (controls > t >= cases). “<”: if the predictor values for the
    # control group are lower or equal than the values of the case group (controls < t
    # <= cases).
    prediction_is_control <- eval(parse(text = paste0(
        "proc_result$predictor ", proc_result$direction, " coords$threshold"
    )))
    predicted_labels <- ifelse(prediction_is_control, proc_result$levels[1], proc_result$levels[2])

    # Compute confusion matrix using caret
    ca_cfm <- caret::confusionMatrix(
        factor(predicted_labels, levels = proc_result$levels),
        factor(proc_result$response, levels = proc_result$levels)
    )

    # Check agreement between caret and pROC's reported accuracy
    if (abs(ca_cfm$overall[["Accuracy"]] - coords$accuracy) > 1e-6) {
        stop("Mismatch between pROC and caret accuracy")
    }

    # Format predictions for log loss (mlr3measures expects both class probabilities)
    tmp_pred <- data.frame("positive" = proc_result$predictor, "negative" = 1 - proc_result$predictor)
    colnames(tmp_pred) <- rev(proc_result$levels) # Ensure correct column order

    logloss_mlr3 <- tryCatch(
        {
            logloss_mlr3 <- mlr3measures::logloss(
                truth = factor(proc_result$response, levels = proc_result$levels),
                prob = as.matrix(tmp_pred)
            )
            logloss_mlr3
        },
        error = function(e) {
            logloss_mlr3 <- NA_real_ # Handle cases where log loss cannot be computed
            logloss_mlr3
        }
    )

    # Append class counts and additional metrics
    return(
        tibble::add_column(
            coords,
            "n_classA" = sum(proc_result$response == proc_result$levels[1]),
            "n_classB" = sum(proc_result$response == proc_result$levels[2]),
            "n_total" = length(proc_result$response),
            "Kappa" = ca_cfm$overall[["Kappa"]],
            "AccuracyNull" = ca_cfm$overall[["AccuracyNull"]],
            "Accuracy.ci95.lower" = ca_cfm$overall[["AccuracyLower"]],
            "Accuracy.ci95.upper" = ca_cfm$overall[["AccuracyUpper"]],
            "AccuracyPValue" = ca_cfm$overall[["AccuracyPValue"]],
            "logloss" = logloss_mlr3
        )
    )
}
