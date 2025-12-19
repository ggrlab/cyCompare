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
#' @param return_as_matrix Logical; if `TRUE`, returns a matrix of distances.
#'   If `FALSE`, returns a `data.table` with columns `sample_A_i`, `sample_B_j`, `dist`, `time`, `sample_A`, `sample_B`.
#' @param intermediate_file
#'  A character string specifying the file to write intermediate results to.
#'  Default is "distances_intermediate.csv".
#'  If `write_intermediate = TRUE`, this file will be created or overwritten.
#'  If `write_intermediate = FALSE`, this parameter is ignored.
#'  If the file exists, it will be overwritten.
#'  If the file does not exist, it will be created.
#'  If `write_intermediate = TRUE`, the file will be written in append mode
#'  (i.e., new results will be added to the end of the file).
#'  If `write_intermediate = FALSE`, the file will not be written.
#' @param ... Additional arguments passed to the `loss` function.
#'
#' @return Invisibly returns `NULL`. Results are written to `distances_intermediate.csv` if `write_intermediate = TRUE`.
#' @details
#' For each combination of `datalist_A[[i]]` and `datalist_B[[j]]`, unless skipped by `should_skip(i, j)`,
#' the function calculates the loss and optionally stores the result with computation time in an intermediate CSV file.
#'
#' If `parallel = TRUE`, the function uses nested `future_lapply` calls from the `future.apply` package.
#'
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
                          intermediate_file = "distances_intermediate.csv",
                          should_skip = function(i, j) i >= j,
                          take_time = FALSE,
                          return_as_matrix = TRUE,
                          ...) {
    sample_A_i <- sample_B_j <- NULL # R CMD check compatibility
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
    # check if future is available
    distances <- NA
    distances <- tryCatch(
        {
            distances <- distances_wrapper(
                applyfun = future.apply::future_apply,
                filtered_combinations = filtered_combinations,
                datalist_A = datalist_A,
                datalist_B = datalist_B,
                loss = loss,
                dt_res = dt_res,
                verbose = verbose,
                write_intermediate = write_intermediate,
                intermediate_file = intermediate_file,
                take_time = take_time,
                ...
            )
        },
        error = function(e) {
            distances <- distances_wrapper(
                applyfun = apply,
                filtered_combinations = filtered_combinations,
                datalist_A = datalist_A,
                datalist_B = datalist_B,
                loss = loss,
                dt_res = dt_res,
                verbose = verbose,
                write_intermediate = write_intermediate,
                intermediate_file = intermediate_file,
                take_time = take_time,
                ...
            )
            warning("parallel distances_wrapper failed, falling back to sequential: \n", e)
            return(distances)
        }
    )

    distances_bound <- data.table::rbindlist(distances)
    # If names were NULL, distances_bound and distances_bound_named are identical
    distances_bound_named <- dplyr::left_join(
        filtered_combinations,
        distances_bound,
        by = c("sample_A_i", "sample_B_j"),
        suffix = c("", "_removeme")
    )
    distances_bound_named <- distances_bound_named[, !grepl("_removeme", names(distances_bound_named))]
    return(
        parse_pairwise_loss(
            bound_or_file = distances_bound_named,
            length_a = length(datalist_A),
            length_b = length(datalist_B),
            return_as_matrix = return_as_matrix
        )
    )
}

distances_wrapper <- function(
    applyfun,
    filtered_combinations,
    datalist_A,
    datalist_B,
    loss,
    sample_Ai_Bj,
    verbose,
    write_intermediate,
    intermediate_file,
    dt_res,
    take_time,
    ...) {
    return(applyfun(
        filtered_combinations[, c("sample_A_i", "sample_B_j")],
        MARGIN = 1, # Apply over rows
        simplify = FALSE,
        FUN = function(sample_Ai_Bj) {
            distfun_pairwise(
                datalist_A = datalist_A,
                datalist_B = datalist_B,
                sample_Ai_Bj = sample_Ai_Bj,
                loss = loss,
                verbose = verbose,
                write_intermediate = write_intermediate,
                intermediate_file = intermediate_file,
                dt_res = dt_res,
                take_time = take_time,
                ...
            )
        }
    ))
}

distfun_pairwise <- function(
    datalist_A,
    datalist_B,
    sample_Ai_Bj,
    loss,
    verbose,
    write_intermediate,
    intermediate_file,
    dt_res,
    take_time,
    ...) {
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
