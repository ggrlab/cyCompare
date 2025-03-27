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



# single_metric_kBET <- function(cellmat, cell_ids, k0 = NULL, testSize = NULL, ...) {
#     # part_testing_ids <- c(1:100, nrow(cellmat):(nrow(cellmat) - 100))
#     res <- kBET::kBET(
#         df = as.matrix(cellmat[part_testing_ids, ]),
#         batch = as.numeric(factor(cell_ids))[part_testing_ids],
#         #  number of nearest neighbours to test on (neighbourhood size)
#         k0 = k0,
#         # number of data points to test, (10 percent sample size default, but at least 25)
#         testSize = testSize,
#         do.pca = FALSE
#     )
#     return(res$average.pval)
# }
# single_metric_kBET_fast <- function(cellmat, cell_ids, n_subsample_per_id = 2500, k0 = NULL, testSize = NULL, ...) {
#     # get a random subset of row numbers per unique cell id
#     if (!is.na(n_subsample_per_id)) {
#         unique_cells_per_cellID <- unlist(tapply(cell_ids, cell_ids, function(x) {
#             return(sample(length(x), size = n_subsample_per_id, replace = FALSE))
#         }), use.names = FALSE)
#     } else {
#         unique_cells_per_cellID <- 1:nrow(cellmat)
#     }

#     res <- kBET::kBET(
#         # Only for speeding up during debugging
#         # df = as.matrix(cellmat[part_testing_ids, ]),
#         # batch = as.numeric(factor(cell_ids))[part_testing_ids],
#         df = as.matrix(cellmat[unique_cells_per_cellID, ]),
#         batch = as.numeric(factor(cell_ids))[unique_cells_per_cellID],
#         #  number of nearest neighbours to test on (neighbourhood size)
#         k0 = k0,
#         # number of data points to test, (10 percent sample size default, but at least 25)
#         testSize = testSize,
#         do.pca = FALSE,
#         adapt = FALSE,
#         heuristic = FALSE,
#         plot = FALSE
#     )
#     return(as.matrix(res$average.pval))
# }
