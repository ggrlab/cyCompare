transformlist_named <- function(transformlist, relevant_columns, flowcore = FALSE) {
    if (length(transformlist) == 1 && !is.null(transformlist)) {
        if (is.function(transformlist)) {
            transformlist <- list(transformlist)
        }
        transformlist <- setNames(rep(transformlist, length(relevant_columns)), relevant_columns)
    }
    if (flowcore) {
        return(
            flowCore::transformList(
                relevant_columns,
                transformlist
            )
        )
    } else {
        return(transformlist)
    }
}
