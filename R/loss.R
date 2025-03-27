#' Compute Earth Mover's Distance (EMD) Between Matrices via Histograms
#'
#' Computes the Earth Mover's Distance (EMD) between corresponding columns of two matrices
#' by binning values into histograms and comparing them using a specified EMD function.
#'
#' @param mat1 A numeric matrix, where each column represents a variable or measurement.
#' @param mat2 A second numeric matrix of the same dimensions as `mat1`.
#' @param n_breaks Integer. Number of histogram bins to use (default is 100).
#' @param emd_fun A function to compute the Earth Mover's Distance between two binned distributions.
#'   Must accept two matrices of equal shape and return a numeric scalar. Default is `emdist::emd2d`.
#' @param ... Additional arguments passed to `emd_fun`.
#'
#' @return A numeric value: the EMD between `mat1` and `mat2` over all columns.
#'   Currently returns only the EMD between the two matrices as a scalar (element [2,1] of the matrix).
#'
#' @details
#' Each matrix is converted into a set of histograms (one per column) using globally shared bin breaks.
#' EMD is then computed between the two matrices by applying `emd_fun()` to the binned data.
#'
#' @examples
#' \dontrun{
#' mat1 <- matrix(rnorm(1000), ncol = 10)
#' mat2 <- matrix(rnorm(1000, mean = 1), ncol = 10)
#' lossfun_emd(mat1, mat2)
#' }
#'
#' @export
lossfun_hist <- function(
    mat1,
    mat2,
    n_breaks = 100,
    lossfun = lossfun_hist_weighted, ...) {
    if (n_breaks < 3) {
        stop("n_breaks must be at least 3. With 2 breaks (-> min and max), all values are in the same bin.")
    }

    # Handle potential override of n_breaks via ...
    extra_args <- list(...)
    if ("n_breaks" %in% names(extra_args)) {
        n_breaks <- extra_args$n_breaks
    }


    # Combine matrices and compute global histogram breaks
    matlist <- list(mat1, mat2)
    global_breaks <- seq(min(sapply(matlist, min)), max(sapply(matlist, max)), length.out = n_breaks)
    dist_between_breaks <- global_breaks[2] - global_breaks[1]

    # Bin each column of each matrix using shared breaks
    broken_histogrammed <- lapply(matlist, function(mat) {
        apply(mat, 2, function(col) {
            # $counts: n integers; for each cell, the number of x[] inside.
            graphics::hist(col, breaks = global_breaks, plot = FALSE)$counts
        })
    })
    # broken_histogrammed is now a list of two matrices which both are
    #   n_breaks X ncol(mat)
    # matrices. Each column represents a histogram of the corresponding column in the input matrix.

    # Compute EMD between the two sets of histograms
    # This is also why the breaks must be shared for the two matrices,
    # otherwise the calculation would not be meaningful.
    return(lossfun(
        broken_histogrammed[[1]],
        broken_histogrammed[[2]],
        ydist = dist_between_breaks,
        ...
    ))
}

#' Wrapper Functions for Histogram-Based Loss Calculation
#'
#' These functions provide specific wrappers around \code{lossfun_hist()} for evaluating
#' differences between two matrices using Earth Mover's Distance (EMD). Both use
#' \code{emdist::emd2d()} with different argument configurations.
#'
#' \code{\link[emdist]{emd2d}} interprets the two matrices A and B as a distibution
#' over a two-dimensional grid. The distance between the grid points in each direction
#' is defined by xdist and ydist.
#' @inheritParams lossfun_hist
#'
#' @details
#' \describe{
#'   \item{\code{lossfun_hist_cytonorm}}{
#'     Mimics the behavior of the CytoNorm implementation by discarding the \code{ydist}
#'     argument before calling \code{emd2d()}.
#'   }
#'   \item{\code{lossfun_hist_weighted}}{
#'     A generic weighted EMD computation that retains all arguments, including \code{ydist},
#'     allowing full control over the distance calculation.
#'   }
#' }
#' @return A numeric value representing the loss (EMD) between the input matrices.
#' @seealso \code{\link[emdist]{emd2d}}, \code{\link{lossfun_hist}}
#'
#' @export
lossfun_hist_cytonorm <- function(x, y, ...) {
    paramlist <- list(...)
    # ydist holds the width of a single bin
    paramlist[["ydist"]] <- NULL
    do.call(emdist::emd2d, c(list(x, y, dist = "euclidean"), paramlist))
}
#' @rdname lossfun_hist_cytonorm
#' @export
lossfun_hist_weighted <- function(x, y, ...) {
    paramlist <- list(...)
    do.call(emdist::emd2d, c(list(x, y, dist = "euclidean"), paramlist))
}



single_metric_kBET <- function(cellmat, cell_ids, k0 = NULL, testSize = NULL, ...) {
    # part_testing_ids <- c(1:100, nrow(cellmat):(nrow(cellmat) - 100))
    res <- kBET::kBET(
        df = as.matrix(cellmat[part_testing_ids, ]),
        batch = as.numeric(factor(cell_ids))[part_testing_ids],
        #  number of nearest neighbours to test on (neighbourhood size)
        k0 = k0,
        # number of data points to test, (10 percent sample size default, but at least 25)
        testSize = testSize,
        do.pca = FALSE
    )
    return(res$average.pval)
}
single_metric_kBET_fast <- function(cellmat, cell_ids, n_subsample_per_id = 2500, k0 = NULL, testSize = NULL, ...) {
    # get a random subset of row numbers per unique cell id
    if (!is.na(n_subsample_per_id)) {
        unique_cells_per_cellID <- unlist(tapply(cell_ids, cell_ids, function(x) {
            return(sample(length(x), size = n_subsample_per_id, replace = FALSE))
        }), use.names = FALSE)
    } else {
        unique_cells_per_cellID <- 1:nrow(cellmat)
    }

    res <- kBET::kBET(
        # Only for speeding up during debugging
        # df = as.matrix(cellmat[part_testing_ids, ]),
        # batch = as.numeric(factor(cell_ids))[part_testing_ids],
        df = as.matrix(cellmat[unique_cells_per_cellID, ]),
        batch = as.numeric(factor(cell_ids))[unique_cells_per_cellID],
        #  number of nearest neighbours to test on (neighbourhood size)
        k0 = k0,
        # number of data points to test, (10 percent sample size default, but at least 25)
        testSize = testSize,
        do.pca = FALSE,
        adapt = FALSE,
        heuristic = FALSE,
        plot = FALSE
    )
    return(as.matrix(res$average.pval))
}
