#' Compute Histogram-Based Distance Between Two Matrices
#'
#' This function converts each column of two numeric matrices into histograms using a shared global binning scheme,
#' then computes a distance (e.g., Earth Mover's Distance) between the histogrammed representations using a supplied
#' `lossfun`. It is useful for comparing sample distributions.
#'
#' See also \code{\link{lossfun_hist_columnwise}} for a column-wise variant. Here, the loss
#' uses emdist::emd2d() to compute the distance between the two binned matrices jointly,
#' understanding each matrix as a 2D distribution over the grid defined by the columns and bins.
#'
#' This makes only sense if the columns are somewhat related, if you cannot directly compare
#' the columns of the two matrices, I _think_ that \code{lossfun_hist_columnwise()} is more useful.
#'
#' @param mat1
#' A numeric matrix. Each column represents a variable (e.g., marker) and rows are observations (e.g., cells).
#' @param mat2
#' A second numeric matrix, same shape and interpretation as `mat1`.
#' @param n_breaks Integer (default: 100). Number of histogram bins to use. Must be at least 3.
#' @param lossfun
#' A function that takes two matrices of binned counts and computes a loss. Default is `lossfun_hist_weighted`.
#' @param ...
#' Additional arguments passed to `lossfun`.

#' @return A numeric value: the lossfun between `mat1` and `mat2` over all columns.
#'   Currently returns only the lossfun between the two matrices as a scalar (element `[2,1]` of the matrix).
#'
#' @details
#' Each matrix is converted into a set of histograms (one per column) using globally shared bin breaks.
#' lossfun is then computed between the two matrices by applying `lossfun()` to the binned data.
#'
#' @examples
#' mat1 <- matrix(rnorm(50), ncol = 10)
#' mat2 <- matrix(rnorm(50, mean = 1), ncol = 10)
#' lossfun_hist(mat1, mat2)
#' lossfun_hist(mat1, mat2, n_breaks = 100, lossfun = lossfun_hist_cytonorm)
#'
#' @export
lossfun_hist <- function(
    mat1,
    mat2,
    n_breaks = 100,
    lossfun = lossfun_hist_weighted,
    quantilebreaks = FALSE,
    ...) {
    if (n_breaks < 3) {
        stop("n_breaks must be at least 3. With 2 breaks (-> min and max), all values are in the same bin.")
    }

    # Combine matrices and compute global histogram breaks
    matlist <- list(mat1, mat2)
    if (quantilebreaks) {
        # Note: This breaks at quantiles ACROSS COLUMNS and MATRICES!
        global_breaks <- quantile(
            unlist(matlist),
            probs = seq(0, 1, length.out = n_breaks)
        )
        dist_between_breaks <- diff(global_breaks)
    } else {
        global_breaks <- seq(min(sapply(matlist, min)), max(sapply(matlist, max)), length.out = n_breaks)
        dist_between_breaks <- global_breaks[2] - global_breaks[1]
    }

    # Bin each column of each matrix using shared breaks
    broken_histogrammed <- lapply(matlist, function(mat) {
        apply(mat, 2, function(col) {
            # $counts: n integers; for each cell, the number of x[] inside.
            # graphics::hist(col, breaks = global_breaks, plot = FALSE)$counts
            graphics::hist(col, breaks = global_breaks, plot = FALSE)$density
        })
    })
    # graphcs::hist says:
    #   values \hat f(x_i), as estimated density values. If all(diff(breaks) == 1),
    # they are the relative frequencies counts/n and in general satisfy \sum_i \hat
    # f(x_i) (b_{i+1}-b_i) = 1, where b_i = breaks[i].

    # test if the densities actually sum up to 1 using dist_between_breaks
    if (!all(
        sapply(broken_histogrammed, function(x) {
            # Multiply each row (bin) with the bin width and sum over all bins
            if (length(dist_between_breaks) > 1) {
                # variable bin widths
                apply(diag(dist_between_breaks) %*% x, 2, sum)
            } else {
                # constant bin width
                apply(x * dist_between_breaks, 2, sum)
            }
        }) - 1 < .Machine$double.eps * 10
    )) {
        warning("Histogram densities do not sum to 1. Check if the breaks are set correctly.")
    }
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
#' In contrast to \code{lossfun_hist()}, which computes the histogram-based loss
#' across all columns jointly, this function computes the loss column-wise and sums
#' the results.
#'
#' NOTE: I am NOT sure if this is better or worse than \code{lossfun_hist()}.
#'
#' @param result_percolumn
#' Logical (default: FALSE). If TRUE, returns a numeric vector of losses per column
#'
#' @inheritParams lossfun_hist
#' @export
lossfun_hist_columnwise <- function(
    mat1,
    mat2,
    n_breaks = 100,
    lossfun = lossfun_hist_weighted,
    quantilebreaks = FALSE,
    result_percolumn = FALSE,
    ...) {
    if (n_breaks < 3) {
        stop("n_breaks must be at least 3. With 2 breaks (-> min and max), all values are in the same bin.")
    }
    dist_per_column <- sapply(seq_len(ncol(mat1)), function(col_i) {
        values <- c(mat1[, col_i], mat2[, col_i])
        if (quantilebreaks) {
            # Note: This breaks at quantiles ACROSS COLUMNS and MATRICES!
            breaks <- quantile(
                values,
                probs = seq(0, 1, length.out = n_breaks)
            )
            dist_between_breaks <- diff(breaks)
        } else {
            breaks <- seq(min(values), max(values), length.out = n_breaks)
            dist_between_breaks <- breaks[2] - breaks[1]
        }

        # Bin each column of each matrix using shared breaks
        histogrammed_1 <- matrix(graphics::hist(mat1[, col_i], breaks = breaks, plot = FALSE)$density, ncol = 1)
        histogrammed_2 <- matrix(graphics::hist(mat2[, col_i], breaks = breaks, plot = FALSE)$density, ncol = 1)

        lossfun(
            histogrammed_1,
            histogrammed_2,
            ydist = dist_between_breaks,
            ...
        )
    })
    if (result_percolumn) {
        return(dist_per_column)
    } else {
        return(sum(dist_per_column))
    }
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
#' @examples
#' mat1 <- matrix(rnorm(50), ncol = 10)
#' mat2 <- matrix(rnorm(50, mean = 1), ncol = 10)
#' lossfun_hist_cytonorm(mat1, mat2)
#' # Usually you would not use this function directly, but rather
#' # use lossfun_hist() which calls this internally.
#' lossfun_hist(mat1, mat2, n_breaks = 100, lossfun = lossfun_hist_cytonorm)
lossfun_hist_cytonorm <- function(mat1, mat2, ...) {
    paramlist <- list(...)
    # ydist holds the width of a single bin
    paramlist[["ydist"]] <- NULL
    do.call(emdist::emd2d, c(list(mat1, mat2, dist = "euclidean"), paramlist))
}
#' @rdname lossfun_hist_cytonorm
#' @inheritParams lossfun_hist
#' @export
#' @examples
#' mat1 <- matrix(rnorm(50), ncol = 10)
#' mat2 <- matrix(rnorm(50, mean = 1), ncol = 10)
#' lossfun_hist_weighted(mat1, mat2)
#' # Usually you would not use this function directly, but rather
#' # use lossfun_hist() which calls this internally.
#' lossfun_hist(mat1, mat2, n_breaks = 100, lossfun = lossfun_hist_weighted)
lossfun_hist_weighted <- function(mat1, mat2, ...) {
    paramlist <- list(...)
    do.call(emdist::emd2d, c(list(mat1, mat2, dist = "euclidean"), paramlist))
}
