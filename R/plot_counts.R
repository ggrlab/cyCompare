plot_counts <- function(dt_counts, populations = NULL) {
    if(all(is.null(populations))) {
        populations <- unique(dt_counts$name)
    }
    browser()
    library(ggplot2)
    # ggplot(
    #     dt_counts[Population %in% populations], 
    #     aes(x = Time, y = Count, col = Device)) + 
    #     geom_point() +
}