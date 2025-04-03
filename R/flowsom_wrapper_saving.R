flowsom_wrapper_saving <- function(fs_applied, dir_out, modelname = "cells_asinh.random10k.cd3") {
    dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)
    qs::qsave(fs_applied[["cells_clusters_from_train"]], file.path(dir_out, paste0(modelname, ".qs")))
    for (name_x in names(fs_applied[["ncells_per_x"]])) {
        data.table::fwrite(fs_applied[["ncells_per_x"]][[name_x]], file.path(dir_out, paste0("ncells_per_", name_x, ".csv")))
    }
    for (name_x in names(fs_applied[["proportions_per_x"]])) {
        data.table::fwrite(fs_applied[["proportions_per_x"]][[name_x]], file.path(dir_out, paste0("proportions_per_", name_x, ".csv")))
    }
    data.table::fwrite(fs_applied[["pheno"]], file.path(dir_out, "pheno.csv"))
}
