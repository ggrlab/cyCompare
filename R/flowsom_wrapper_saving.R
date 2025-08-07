#' Save FlowSOM Prediction Results to Disk
#'
#' Saves the output of a FlowSOM prediction run, including cluster assignments,
#' cell counts, cluster proportions, and phenotype metadata to the specified output directory.
#'
#' @param fs_applied
#' A list of FlowSOM prediction results, typically returned by `flowsom_repeatsubsampling_apply()`.
#'    Expected to contain the elements:
#'    - `"cells_clusters_from_train"`: List of predicted clusters per cell.
#'    - `"ncells_per_x"`: Named list of data.tables with cluster cell counts per sample.
#'    - `"proportions_per_x"`: Named list of data.tables with normalized cluster proportions.
#' @param dir_out
#' A character string specifying the directory to which files should be saved.
#' Will be created if it doesn't exist.
#' @param modelname
#' A string used as the base name for the saved `.qs` file of cell-level cluster
#' assignments (default: `"cells_asinh.random10k.cd3"`).
#'
#' @return No return value. Side effect: writes files to disk in the specified output folder.
#'
#' @export
flowsom_wrapper_saving <- function(fs_applied, dir_out, modelname = "cells_asinh.random10k.cd3") {
    # Create output directory if it doesn't exist
    dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)

    # Save cell-level predicted cluster assignments as .qs file
    qs::qsave(
        fs_applied[["cells_clusters_from_train"]],
        file.path(dir_out, paste0(modelname, ".qs"))
    )

    # Save cell count tables (cluster cell counts per sample) as CSVs
    for (name_x in names(fs_applied[["ncells_per_x"]])) {
        data.table::fwrite(
            fs_applied[["ncells_per_x"]][[name_x]],
            file.path(dir_out, paste0("ncells_per_", name_x, ".csv"))
        )
    }

    # Save normalized cluster proportions per sample as CSVs
    for (name_x in names(fs_applied[["proportions_per_x"]])) {
        data.table::fwrite(
            fs_applied[["proportions_per_x"]][[name_x]],
            file.path(dir_out, paste0("proportions_per_", name_x, ".csv"))
        )
    }

}
