transformlist_named <- function(transformlist, relevant_columns) {
    if (length(transformlist) == 1 && !is.null(transformlist)) {
        if (is.function(transformlist)) {
            transformlist <- list(transformlist)
        }
        transformlist <- setNames(rep(transformlist, length(relevant_columns)), relevant_columns)
    }
    return(transformlist)
}
