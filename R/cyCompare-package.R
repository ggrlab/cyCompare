#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom data.table :=
#' @importFrom data.table .BY
#' @importFrom data.table .EACHI
#' @importFrom data.table .GRP
#' @importFrom data.table .I
#' @importFrom data.table .N
#' @importFrom data.table .NGRP
#' @importFrom data.table .SD
#' @importFrom data.table data.table
## usethis namespace: end
NULL


.onLoad <- function(libname, pkgname) {
    # Load namespace quietly; triggers mlr3learners::.onLoad() -> registers learners
    if (!requireNamespace("mlr3learners", quietly = TRUE)) {
        # nothing else needed; mlr3learners::.onLoad() handles the registration
    }
}
