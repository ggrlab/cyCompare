lossfun_emd <- function(
    mat1,
    mat2,
    n_breaks = 100,
    emd_fun = function(x, y, ...) {
        emdist::emd2d(x, y, dist = "euclidean", ...)
    }, ...) {
    if ("n_breaks" %in% names(list(...))) {
        n_breaks <- list(...)$n_breaks
    }
    matlist <- list(mat1, mat2)
    global_breaks <- seq(min(sapply(matlist, min)), max(sapply(matlist, max)), length.out = n_breaks)
    dist_between_breaks <- global_breaks[2] - global_breaks[1]
    broken_histogrammed <- list()
    for (mat_i in seq_along(matlist)) {
        broken_histogrammed[[mat_i]] <- apply(
            matlist[[mat_i]], 2, function(x) {
                hist(x, breaks = global_breaks, plot = FALSE)$counts
            }
        )
    }
    emd_dists <- matrix(NA, nrow = length(broken_histogrammed), ncol = length(broken_histogrammed))
    for (cell_id_i in seq_along(broken_histogrammed)) {
        for (cell_id_j in seq_along(broken_histogrammed)) {
            if (cell_id_i > cell_id_j) {
                emd_dists[cell_id_i, cell_id_j] <- emd_fun(
                    broken_histogrammed[[cell_id_i]],
                    broken_histogrammed[[cell_id_j]],
                    ydist = dist_between_breaks,
                    ...
                )
            }
        }
    }
    # diag(emd_dists) <- 0
    # emd_dists[upper.tri(emd_dists)] <- t(emd_dists)[upper.tri(emd_dists)]
    # return(emd_dists)
    return(emd_dists[2, 1, drop = TRUE])
}
lossfun_emd_cytonorm <- function(mat1, mat2, n_breaks = 100, ...) {
    lossfun_emd(mat1, mat2, n_breaks, emd_fun = function(x, y, ...) {
        paramlist <- list(...)
        paramlist[["ydist"]] <- NULL
        do.call(emdist::emd2d, c(list(x, y, dist = "euclidean"), paramlist))
    }, ...)
}
lossfun_emd_weighted <- function(mat1, mat2, n_breaks = 100, ...) {
    lossfun_emd(mat1, mat2, n_breaks, emd_fun = function(x, y, ...) {
        paramlist <- list(...)
        do.call(emdist::emd2d, c(list(x, y, dist = "euclidean"), paramlist))
    }, ...)
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
