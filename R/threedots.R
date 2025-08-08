#' Access Internal Functions in a Package Namespace
#'
#' This infix operator provides access to non-exported (internal) functions
#' from a package's namespace, similar to using `pkg:::fun` in base R.
#' Gotten from:
#'  https://stat.ethz.ch/pipermail/r-devel/2013-August/067210.html
#'
#' @param pkg A character string. The name of the package.
#' @param fun A character string. The name of the internal function to access.
#' @return The function object retrieved from the specified package namespace.
#'
#' @name REEXPORT-ns-triplecolon
#'
#' @details
#' I need this in for the functions defined in the this very script "threedots.R".
#' This is a workaround to access internal functions that are not exported by the package.
#' It allows you to use the syntax `pkg %:::% fun` to retrieve the function `fun` from the package `pkg`.
#'
#' @examples
#' \dontrun{
#' # Access the internal .unitize function from the grid package
#' "grid" %:::% ".unitize"
#' }
#' @keywords re-export
`%:::%` <- function(pkg, fun) {
    # Retrieve an internal (non-exported) function from a package's namespace
    get(fun,
        envir = asNamespace(pkg),
        inherits = FALSE
    )
}

#' autoplot.prcomp
#'
#' This function is a wrapper for the `autoplot.prcomp` function from the `ggfortify` package.
#'
#' @param ... Arguments passed to the `autoplot.prcomp` function.
#' @return The result of the `autoplot.prcomp` function.
#' @keywords re-export
autoplot.prcomp <- function(...) {
    `%:::%`("ggfortify", "autoplot.prcomp")(...)
}
