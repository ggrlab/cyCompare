#' Compute Pairwise Loss Between Two Lists of Data
#'
#' Computes a pairwise loss matrix between two lists of data objects using a user-defined loss function.
#' Optionally supports parallelization, skipping combinations, verbose output, and intermediate result writing.
#'
#' @param datalist_A
#'  A list of data objects (e.g., matrices, data frames, vectors).
#'  Should be such that the loss can be computed between any two elements.
#'  Usually I expect a list of numeric matrices.
#'
#' This list corresponds to the rows of the output matrix.
#' @param datalist_B
#'  A list of data objects (e.g., matrices, data frames, vectors).
#'  Should be such that the loss can be computed between any two elements.
#'  Usually I expect a list of numeric matrices.
#'
#' This list corresponds to the columns of the output matrix.
#' @param loss
#'  A function taking two arguments (elements from `datalist_A` and `datalist_B`) and
#'  returning a (numeric) loss.
#' @param verbose Logical; if `TRUE`, progress messages during execution.
#' @param write_intermediate
#'  Logical; if `TRUE`, writes each computed result to a CSV file (`distances_intermediate.csv`).
#' The file looks like this:
#' \preformatted{
#' sample_A_i,sample_B_j,dist,time
#' 1,1,0.123,0.456
#' 1,2,0.789,0.012
#' ...
#' }
#' @param should_skip
#' A function with signature `function(i, j)` returning `TRUE` if the pair (i, j)
#' should be skipped and its loss NOT calculated.
#' By default, the function skips all pairs where `i >= j`. This is useful for symmetric loss functions
#' where loss(A,A) = 0 and loss(A,B) = loss(B,A).
#'
#' Useful should_skip functions are:
#' - `function(i, j) i >= j` to skip the upper triangle and diagonal of all possible pairs.
#' - `function(i, j) FALSE` to calculate all possible pairs.
#' - `function(i, j) i == j` to skip the diagonal
#' - `function(i, j) i == 1` to skip the first row
#' - `function(i, j) i == 1 || j == 1` to skip the first row and column
#'
#' @param take_time Logical; if `TRUE`, take timing information.
#' @param ... Additional arguments passed to the `loss` function.
#'
#' @return Invisibly returns `NULL`. Results are written to `distances_intermediate.csv` if `write_intermediate = TRUE`.
#' @details
#' For each combination of `datalist_A[[i]]` and `datalist_B[[j]]`, unless skipped by `should_skip(i, j)`,
#' the function calculates the loss and optionally stores the result with computation time in an intermediate CSV file.
#'
#' If `parallel = TRUE`, the function uses nested `future_lapply` calls from the `future.apply` package.
#'
#' @importFrom dplyr filter
#' @importFrom data.table data.table fwrite
#' @importFrom future.apply future_lapply
#' @examples
#' # Define dummy data and a simple loss function
#' A <- list(a = 1:5, b = 6:10)
#' B <- list(c = 1:5, d = 6:10)
#' loss_fn <- function(x, y) sum((x - y)^2)
#'
#' # Compute pairwise loss
#' loss_pairwise(A, B, loss = loss_fn, verbose = FALSE, write_intermediate = FALSE)
#'
#' @export
loss_pairwise <- function(datalist_A,
                          datalist_B,
                          loss = lossfun_hist,
                          verbose = FALSE,
                          write_intermediate = FALSE,
                          should_skip = function(i, j) i >= j,
                          take_time = FALSE,
                          return_as_matrix = TRUE,
                          ...) {
    ### 1. Checks
    # 1.1. Check if both datalists are actually lists
    if (!is.list(datalist_A) || !is.list(datalist_B)) {
        stop("Both datalists must be lists")
    }
    # 1.2. Check if the loss function is a function
    if (!is.function(loss)) {
        stop("The loss function must be a function")
    }
    # 1.3. Check if the skipping function is a function
    if (!is.function(should_skip)) {
        stop("The skipping function must be a function")
    }

    ### 2. Preparation
    # 2.1 Prepare the output matrices
    dist_mat <- matrix(NA, nrow = length(datalist_A), ncol = length(datalist_B))
    rownames(dist_mat) <- names(datalist_A)
    colnames(dist_mat) <- names(datalist_B)
    # time_mat is a matrix of the same size as dist_mat, but it contains the time it took to compute each distance
    time_mat <- dist_mat

    # 2.2 Generate the combinations which we will actually compute based on the skipping function
    filtered_combinations <- expand.grid(
        sample_A_i = seq_along(datalist_A),
        sample_B_j = seq_along(datalist_B)
    ) |>
        dplyr::filter(!should_skip(sample_A_i, sample_B_j)) |>
        dplyr::mutate(
            # We use indices anyways, but if the names are NULL, the columns are also NULL --> not present
            sample_A = names(datalist_A)[sample_A_i],
            sample_B = names(datalist_B)[sample_B_j]
        )


    ### 3. Compute the distances
    intermediate_file <- "distances_intermediate.csv"
    dt_res <- data.table::data.table(
        "sample_A_i" = NA_integer_,
        "sample_B_j" = NA_integer_,
        "dist" = NA_real_,
        "time" = NA_real_,
        "sample_A" = NA_character_,
        "sample_B" = NA_character_
    )
    if (write_intermediate) {
        if (file.exists(intermediate_file)) {
            warning("Intermediate file already exists. Overwriting it with an empty file.")
        }
        data.table::fwrite(
            dt_res[0, ],
            intermediate_file,
            # Append = False results in a new, empty intermediate file!
            append = FALSE
        )
    }

    distances <- future.apply::future_apply(
        filtered_combinations[, c("sample_A_i", "sample_B_j")],
        MARGIN = 1, # Apply over rows
        function(sample_Ai_Bj) {
            sample_A_i <- sample_Ai_Bj[["sample_A_i"]]
            sample_B_j <- sample_Ai_Bj[["sample_B_j"]]
            if (verbose) {
                message(paste0(Sys.time(), " ", sample_A_i, " vs ", sample_B_j), quote = FALSE, end = "   ")
            }
            calculation_time <- NA
            if (take_time) {
                start <- Sys.time()
            }
            current_dist <- loss(datalist_A[[sample_A_i]], datalist_B[[sample_B_j]], ...)
            if (take_time) {
                end <- Sys.time()
                calculation_time <- as.numeric(end - start, units = "secs")
            }
            # Create a copy of the data.table
            dt_current <- data.table::data.table(dt_res)
            dt_current[["sample_A_i"]] <- sample_A_i
            dt_current[["sample_B_j"]] <- sample_B_j
            dt_current[["dist"]] <- current_dist
            dt_current[["time"]] <- calculation_time
            try(dt_current[["sample_A"]] <- names(datalist_A)[sample_A_i], silent = TRUE)
            try(dt_current[["sample_B"]] <- names(datalist_B)[sample_B_j], silent = TRUE)
            if (verbose) {
                message(paste0(
                    "(t=", calculation_time, ")   ",
                    current_dist
                ))
            }

            if (write_intermediate) {
                data.table::fwrite(
                    dt_current,
                    intermediate_file,
                    append = TRUE,
                    col.names = FALSE
                )
                message("    Wrote distances_intermediate.csv")
            }
            return(dt_current)
        }
    )

    distances_bound <- data.table::rbindlist(distances)
    # If names were NULL, distances_bound and distances_bound_named are identical
    distances_bound_named <- dplyr::left_join(
        filtered_combinations,
        distances_bound,
        by = c("sample_A_i", "sample_B_j")
    )
    return(
        parse_pairwise_loss(
            bound_or_file = distances_bound_named,
            length_a = length(datalist_A),
            length_b = length(datalist_B),
            return_as_matrix = return_as_matrix
        )
    )
}


# metric_batches <- function(..., reference_mat_meta = NULL, single_metric = single_metric_emd, parallel = TRUE, verbose = TRUE, intermediate_dir = NA) {
#     if (parallel) {
#         metric_pairwise_fun <- metric_pairwise_foreach
#     } else {
#         metric_pairwise_fun <- metric_pairwise
#     }
#     vals <- list(...)
#     if (!is.null(reference_mat_meta)) {
#         if ("reference" %in% names(vals)) {
#             stop("Cannot specify both reference_mat_meta and reference")
#         } else {
#             vals[["reference"]] <- reference_mat_meta
#         }
#     }
#     if (any(names(vals) == "")) {
#         stop("All arguments must be named")
#     }
#     if (anyDuplicated(names(vals)) != 0) {
#         stop("All argument names must be unique")
#     }

#     # I assume that every element of vals is a list of length 2
#     # The first element is a list of multiple matrices \in R^{n_cells x p_parameters}
#     # The second element is a vector of meta values
#     vals_formatted <- lapply(vals, function(vals_X) {
#         if (length(vals_X) != 2) {
#             stop("Each given argument must be a list of length 2: A list of matrices (n_cells x p_parameters) and a vector of meta values (N matrices)")
#         }
#         if (!"list" %in% class(vals_X[[1]])) {
#             stop("The first element of each argument must be a list of matrices")
#         }
#         if (length(vals_X[[2]]) != length(vals_X[[1]])) {
#             stop("The second element of each argument must be a vector of length ", length(x[[1]]))
#         }
#         mats_meta <- list()
#         for (mat_i in seq_len(length(vals_X[[1]]))) {
#             mats_meta[[mat_i]] <- data.table::data.table(vals_X[[1]][[mat_i]])
#             mats_meta[[mat_i]][, "batch_ID" := vals_X[[2]][mat_i]]
#         }
#         return(data.table::rbindlist(mats_meta))
#     })

#     if (intermediate_dir == TRUE) {
#         intermediate_dir <- paste0("intermediate_pairwise_", format(Sys.time(), "%Y-%m-%d_%H.%M.%S"))
#     }
#     if (!is.na(intermediate_dir)) {
#         dir.create(intermediate_dir, recursive = TRUE)
#     }

#     res <- list()
#     do_break <- "reference" %in% names(vals)
#     for (arg_i in seq_len(length(vals_formatted))) {
#         res[[names(vals_formatted)[arg_i]]] <- list()
#         for (arg_j in seq_len(length(vals_formatted))) {
#             if (!is.na(intermediate_dir)) {
#                 current_intermediate_dir <- file.path(intermediate_dir, paste0(names(vals_formatted)[arg_i], "_vs_", names(vals_formatted)[arg_j]))
#                 dir.create(current_intermediate_dir, recursive = TRUE)
#             }

#             if ("reference" %in% names(vals)) {
#                 arg_j <- which(names(vals) == "reference")
#                 # Then compare all list elements to the reference
#                 # Otherwise compare all list elements to all other list elements
#             }

#             res[[names(vals_formatted)[arg_i]]][[names(vals_formatted)[arg_j]]] <- metric_pairwise_fun(
#                 mat_1 = vals_formatted[[arg_i]][, !"batch_ID", with = FALSE],
#                 meta_1 = vals_formatted[[arg_i]][, "batch_ID", with = FALSE],
#                 mat_2 = vals_formatted[[arg_j]][, !"batch_ID", with = FALSE],
#                 meta_2 = vals_formatted[[arg_j]][, "batch_ID", with = FALSE],
#                 single_metric = single_metric,
#                 verbose = verbose,
#                 intermediate_dir = current_intermediate_dir
#             )
#             if (do_break) {
#                 # This is the case where we only want to compare to the reference
#                 break
#             }
#         }
#     }
#     return(res)
# }
