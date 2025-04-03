coords_helper <- function(proc_result, threshold = NA) {
    if (is.na(threshold)) {
        coords <- pROC::coords(
            proc_result,
            "best",
            ret = "all",
            best.method = "closest.topleft",
        )
    } else {
        coords <- pROC::coords(
            proc_result,
            x = threshold,
            input = "threshold",
            ret = "all",
        )
    }
    coords <- tibble::as_tibble(coords)
    coords <- coords[1, ] # In case of multiple best points, take the first
    coords[["proc"]] <- list(proc_result)
    coords[["auc"]] <- pROC::auc(proc_result)[[1]]
    # auc_ci <- as.numeric(pROC::ci.auc(proc_res, method = "bootstrap", boot.n = 10000))
    # 8: In ci.auc.roc(proc_result) :
    # ci.auc() of a ROC curve with AUC == 1 is always 1-1 and can be misleading.
    auc_ci <- suppressWarnings(as.numeric(pROC::ci.auc(proc_result)))
    coords[["auc_ci.low"]] <- auc_ci[1]
    coords[["auc_ci.high"]] <- auc_ci[3]
    coords <- coords |> dplyr::relocate(proc, auc, auc_ci.low, auc_ci.high)

    # Now add caret results
    # confusionMatrix(pred, truth)
    #  “>” (default for multivariate curves): if the predictor values for the control group are higher
    # than the values of the case group (controls > t >= cases). “<”: if the predictor values for the
    # control group are lower or equal than the values of the case group (controls < t
    # <= cases).
    prediction_is_control <- eval(parse(text = paste0(
        "proc_result$predictor ", (proc_result$direction), " coords$threshold"
    )))
    predicted_labels <- ifelse(prediction_is_control, proc_result$levels[1], proc_result$levels[2])
    ca_cfm <- caret::confusionMatrix(
        factor(predicted_labels, levels = proc_result$levels),
        factor(proc_result$response, levels = proc_result$levels)
    )
    if (abs(ca_cfm$overall[["Accuracy"]] - coords$accuracy) > 1e-6) {
        stop("Mismatch between pROC and caret accuracy")
    }

    tmp_pred <- data.frame("positive" = proc_result$predictor, "negative" = 1 - proc_result$predictor)
    colnames(tmp_pred) <- rev(proc_result$levels)
    logloss_mlr3 <- mlr3measures::logloss(
        truth = factor(proc_result$response, levels = proc_result$levels),
        prob = as.matrix(tmp_pred)
    )

    return(
        tibble::add_column(
            coords,
            "n_classA" = sum(proc_result$response == proc_result$levels[1]),
            "n_classB" = sum(proc_result$response == proc_result$levels[2]),
            "n_total"= length(proc_result$response),
            "Kappa" = ca_cfm$overall[["Kappa"]],
            "AccuracyNull" = ca_cfm$overall[["AccuracyNull"]],
            "Accuracy.ci95.lower" = ca_cfm$overall[["AccuracyLower"]],
            "Accuracy.ci95.upper" = ca_cfm$overall[["AccuracyUpper"]],
            "AccuracyPValue" = ca_cfm$overall[["AccuracyPValue"]],
            "logloss" = logloss_mlr3
            # "Accuracy_caret" = ca_cfm$overall[["Accuracy"]],
        )
    )
}
