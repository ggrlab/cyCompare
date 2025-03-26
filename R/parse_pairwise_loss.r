#' Parse, (read) and Format Pairwise Loss Data
#'
#' Reads or processes a data frame or CSV file containing pairwise distance (loss) values
#' between indexed samples and optionally reshapes the data into matrix form.
#'
#' @param bound_or_file
#' A character string giving the path to a CSV file, or a data.frame-like
#' object with columns `sample_A_i`, `sample_B_j`, `dist`, and optionally `time`. The indices
#' `sample_A_i` and `sample_B_j` are interpreted as explicit matrix positions.
#' Usually a result of `loss_pairwise`.
#' @param return_as_matrix
#' Logical. If `TRUE` (default), the function returns a list containing
#' numeric matrices of distances and times. If `FALSE`, the original data is returned as-is.
#' @param length_a
#' Integer. Required if `bound_or_file` is a data.frame. Otherwise, inferred
#' from the maximum `sample_A_i` value if not provided.
#' @param length_b
#' Integer. Required if `bound_or_file` is a data.frame. Otherwise, inferred
#' from the maximum `sample_B_j` value if not provided.
#'
#' @return If `return_as_matrix = TRUE`, returns a list with two components:
#' \describe{
#' \item{`dist`}{A numeric matrix of pairwise distances, indexed by `sample_A_i` and `sample_B_j`.}
#' \item{`time`}{A numeric matrix of associated times, if available.}
#' }
#' If `return_as_matrix = FALSE`, returns a data.frame in long format.
#'
#' @details
#' If the input data contains `sample_A` and `sample_B` columns (e.g., character sample labels),
#' these are used to assign row and column names to the resulting distance matrix.
#'
#' Matrix entries that are not represented in the input data are set to `NA`.
#'
#' @examples
#' \dontrun{
#' # From file
#' parse_pairwise_loss("distances_intermediate.csv")
#'
#' # From in-memory data.frame
#' df <- data.frame(
#'     sample_A_i = c(1, 1, 2), sample_B_j = c(2, 3, 3),
#'     dist = c(0.4, 0.3, 0.5), time = c(0.1, 0.2, 0.15)
#' )
#' parse_pairwise_loss(df, return_as_matrix = TRUE, length_a = 2, length_b = 3)
#' }
#'
#' @export
parse_pairwise_loss <- function(
    bound_or_file = "distances_intermediate.csv",
    return_as_matrix = TRUE,
    length_a = NULL,
    length_b = NULL) {
    if (is.character(bound_or_file)) {
        if (!file.exists(bound_or_file)) {
            stop("Intermediate file does not exist")
        }
        bound_matrix <- data.table::fread(bound_or_file)
        if (return_as_matrix) {
            if (is.null(length_a)) {
                warning("length_a not specified. Trying to infer it from the data., If the calculation was incomplete, this might be wrong.")
                length_a <- max(bound_matrix[["sample_A_i"]])
            }
            if (is.null(length_b)) {
                warning("length_b not specified. Trying to infer it from the data. If the calculation was incomplete, this might be wrong.")
                length_b <- max(bound_matrix[["sample_A_i"]])
            }
        }
    } else {
        bound_matrix <- bound_or_file

        if (is.null(length_a)) {
            stop("If bound_or_file is a matrix, length_a must be specified")
        }
        if (is.null(length_b)) {
            stop("If bound_or_file is a matrix, length_b must be specified")
        }
    }

    if (return_as_matrix) {
        retmat <- matrix(NA, nrow = length_a, ncol = length_b)
        timemat <- retmat
        for (row_i in seq_len(nrow(bound_matrix))) {
            retmat[bound_matrix[row_i, "sample_A_i"], bound_matrix[row_i, "sample_B_j"]] <- bound_matrix[row_i, "dist"]
            timemat[bound_matrix[row_i, "sample_A_i"], bound_matrix[row_i, "sample_B_j"]] <- bound_matrix[row_i, "time"]
        }
        if ("sample_A" %in% colnames(bound_matrix)) {
            rownames(retmat) <- bound_matrix[["sample_A"]]
        }
        if ("sample_B" %in% colnames(bound_matrix)) {
            colnames(retmat) <- bound_matrix[["sample_B"]]
        }
        return(list("dist" = retmat, "time" = timemat))
    } else {
        return(bound_matrix)
    }
}
