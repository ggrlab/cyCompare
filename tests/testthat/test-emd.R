devtools::load_all()
test_that("EMD basic test", {
    n_dims <- 3
    n_points <- 1e2 * n_dims
    datalist_1 <- list(
        matrix(rnorm(n_points), ncol = n_dims),
        matrix(rnorm(n_points), ncol = n_dims)
    )
    datalist_2 <- list(
        matrix(rnorm(n_points), ncol = n_dims),
        matrix(rnorm(n_points), ncol = n_dims),
        matrix(rnorm(n_points), ncol = n_dims)
    )
    l1 <- loss_pairwise(
        datalist_1,
        datalist_2,
        loss = lossfun_hist,
        verbose = FALSE
    )
    # both time and distance matrices should be of the same size (here 2x3)
    testthat::expect_equal(
        dim(l1[[1]]),
        c(length(datalist_1), length(datalist_2))
    )
    testthat::expect_equal(
        dim(l1[[2]]),
        c(2, 3)
    )
    testthat::expect_true(is.numeric(l1[[1]]))
    testthat::expect_true(is.matrix(l1[[1]]))
})

test_that("EMD potential parameter test", {
    n_dims <- 3
    n_points <- 1e2 * n_dims
    datalist_1 <- list(
        matrix(rnorm(n_points), ncol = n_dims),
        matrix(rnorm(n_points), ncol = n_dims)
    )
    datalist_2 <- list(
        matrix(rnorm(n_points), ncol = n_dims),
        matrix(rnorm(n_points), ncol = n_dims),
        matrix(rnorm(n_points), ncol = n_dims)
    )
    param_options <- list(
        "datalist_A" = list(datalist_1),
        "datalist_B" = list(datalist_1, datalist_2),
        "lossfun" = list(lossfun_hist_cytonorm, lossfun_hist_weighted),
        "verbose" = c(TRUE, FALSE),
        "write_intermediate" = c(TRUE, FALSE),
        "should_skip" = list(\(i, j) i >= j, \(i, j) FALSE),
        "take_time" = c(TRUE, FALSE),
        "return_as_matrix" = c(TRUE, FALSE),
        "n_breaks" = c(3, 100)
    )
    all_parameter_combinations <- expand.grid(sapply(param_options, function(x) {
        seq_along(x)
    }))
    for (
        combination_i in seq_len(nrow(all_parameter_combinations))
    ) {
        cat(sprintf("Combination %d/%d\n", combination_i, nrow(all_parameter_combinations)))
        current_selection <- all_parameter_combinations[combination_i, ]

        current_parameters <- sapply(names(param_options), simplify = FALSE, function(param_name) {
            param_options[[param_name]][current_selection[[param_name]]][[1]]
        })
        # not necessary, but otherwise a warning is thrown
        unlink("distances_intermediate.csv", force = TRUE)
        tryCatch(
            {
                do.call(loss_pairwise, current_parameters)
            },
            error = function(e) {
                message(current_parameters)
                stop(e)
            }
        )
    }
    testthat::expect_true(TRUE)
})

test_that("EMD n_breaks minimum", {
    n_dims <- 3
    n_points <- 1e2 * n_dims
    datalist_1 <- list(
        matrix(rnorm(n_points), ncol = n_dims),
        matrix(rnorm(n_points), ncol = n_dims)
    )
    datalist_2 <- list(
        matrix(rnorm(n_points), ncol = n_dims),
        matrix(rnorm(n_points), ncol = n_dims),
        matrix(rnorm(n_points), ncol = n_dims)
    )
    testthat::expect_error(
        {
            loss_pairwise(
                datalist_1,
                datalist_2,
                loss = lossfun_hist,
                n_breaks = 2,
                verbose = FALSE
            )
        },
        regexp = "n_breaks must be at least 3."
    )
    loss_pairwise(
        datalist_1,
        datalist_2,
        loss = lossfun_hist,
        n_breaks = 3,
        verbose = FALSE
    )
    testthat::expect_true(TRUE)
})
