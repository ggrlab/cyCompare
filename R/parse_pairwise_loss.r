#' Parse and Format Pairwise Loss Data from File or Data Frame
#'
#' This function reads pairwise loss or distance data—typically computed via `loss_pairwise()`-
#' and optionally reshapes it into numeric matrix form. It supports both file-based input (CSV)
#' and in-memory `data.frame`/`data.table` input.
#'
#' @param bound_or_file
#' Either a path to a CSV file (as character) or a `data.frame`-like object.
#' Must contain columns `sample_A_i`, `sample_B_j`, and `dist`. The `time` column is optional.
#' Index columns `sample_A_i` and `sample_B_j` define matrix positions.
#' If `sample_A` and `sample_B` columns (e.g. names) are present, they are used as row and column names.
#'
#' Usually a result of `loss_pairwise()`.
#' @param return_as_matrix
#' Logical. If `TRUE` (default), returns a list of matrices (`dist`, `time`).
#'   If `FALSE`, returns the input as-is in long format.
#' @param length_a
#' Integer. Required when `bound_or_file` is a `data.frame` and `return_as_matrix = TRUE`.
#'   Defines number of rows in resulting matrix. Inferred from `sample_A_i` if `bound_or_file` is a file.
#' @param length_b
#' Integer. Required when `bound_or_file` is a `data.frame` and `return_as_matrix = TRUE`.
#'   Defines number of columns in resulting matrix. Inferred from `sample_B_j` if `bound_or_file` is a file.
#'
#' @return If `return_as_matrix = TRUE`, returns a named list:
#' \describe{
#'   \item{`dist`}{Numeric matrix of pairwise distances.}
#'   \item{`time`}{Numeric matrix of corresponding computation times (if available).}
#' }
#' If `return_as_matrix = FALSE`, returns the original `data.frame` in long format.
#'
#' @details
#' If the input data contains `sample_A` and `sample_B` columns (e.g., character sample labels),
#' these are used to assign row and column names to the resulting distance matrix.
#'
#' Matrix entries that are not represented in the input data are set to `NA`.
#'
#' @examples
#' df <- data.frame(
#'     sample_A_i = c(1, 1, 2), sample_B_j = c(2, 3, 3),
#'     dist = c(0.4, 0.3, 0.5), time = c(0.1, 0.2, 0.15)
#' )
#' parse_pairwise_loss(df, return_as_matrix = TRUE, length_a = 2, length_b = 3)
#' \dontrun{
#' # Example 1: Read and convert from CSV file
#' parse_pairwise_loss("distances.csv")
#' }
#'
#' @export
parse_pairwise_loss <- function(
    bound_or_file = "distances_intermediate.csv",
    return_as_matrix = TRUE,
    length_a = NULL,
    length_b = NULL) {
    # Step 1: Read or use in-memory data
    if (is.character(bound_or_file)) {
        if (!file.exists(bound_or_file)) {
            stop("Intermediate file does not exist")
        }
        bound_matrix <- data.table::fread(bound_or_file)

        # Infer matrix dimensions if not specified
        if (return_as_matrix) {
            wrongtext <- " not specified. Trying to infer it from the data., If the calculation was incomplete, this might be wrong."
            if (is.null(length_a)) {
                warning(
                    paste0("length_a", wrongtext)
                )
                length_a <- max(bound_matrix[["sample_A_i"]])
            }
            if (is.null(length_b)) {
                warning(paste0("length_b", wrongtext))
                length_b <- max(bound_matrix[["sample_B_j"]])
            }
        }
    } else {
        bound_matrix <- bound_or_file

        if (is.null(length_a)) {
            stop("If bound_or_file is a matrix or data.frame, length_a must be specified")
        }
        if (is.null(length_b)) {
            stop("If bound_or_file is a matrix or data.frame, length_b must be specified")
        }
    }

    # Step 2: Return long format directly
    if (!return_as_matrix) {
        return(bound_matrix)
    }

    # Step 3: Build output matrices
    retmat <- matrix(NA_real_, nrow = length_a, ncol = length_b)
    rownames(retmat) <- rep(NA_character_, length_a)
    colnames(retmat) <- rep(NA_character_, length_b)
    timemat <- retmat

    for (row_i in seq_len(nrow(bound_matrix))) {
        row_i_rowA <- bound_matrix[row_i, "sample_A_i"]
        row_i_colB <- bound_matrix[row_i, "sample_B_j"]

        retmat[row_i_rowA, row_i_colB] <- bound_matrix[row_i, "dist"]
        if ("time" %in% names(bound_matrix)) {
            timemat[row_i_rowA, row_i_colB] <- bound_matrix[row_i, "time"]
        }

        if ("sample_A" %in% colnames(bound_matrix)) {
            rownames(retmat)[row_i_rowA] <- bound_matrix[["sample_A"]][row_i_rowA]
        }
        if ("sample_B" %in% colnames(bound_matrix)) {
            colnames(retmat)[row_i_colB] <- bound_matrix[["sample_B"]][row_i_colB]
        }
    }

    dimnames(timemat) <- dimnames(retmat)
    return(list("dist" = retmat, "time" = timemat))
}
