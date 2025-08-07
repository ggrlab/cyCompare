prepare_testdata <- function() {
    File <- Cytometer <- Donor <- Sample <- NULL

    cn <- c("FITC-A", "PE-A", "ECD-A", "PB-A", "KrO-A")
    set.seed(42)
    cytobench::simulate_ff
    fs_1 <- cytobench::simulate_fs(n_samples = 3, ncells = 500, flowcore = TRUE, columns = cn)
    flowCore::sampleNames(fs_1) <- sprintf(
        "StudyA_d%02d_s%02d_CyA", 1:3, 1:3
    )
    set.seed(43)
    fs_2 <- cytobench::simulate_fs(n_samples = 3, ncells = 500, flowcore = TRUE, columns = cn)
    fs_2 <- flowCore::fsApply(fs_2, function(x) {
        flowCore::exprs(x) <- flowCore::exprs(x) + 5
        x
    })
    flowCore::sampleNames(fs_2) <- sprintf(
        "StudyA_d%02d_s%02d_CyB", 1:3, 1:3
    )
    df <- tibble::tibble(
        File = c(flowCore::sampleNames(fs_1), flowCore::sampleNames(fs_2))
    ) |>
        tidyr::separate_wider_delim(
            File,
            delim = "_",
            names = c("Study", "Donor", "Sample", "Cytometer"),
            cols_remove = FALSE
        ) |>
        dplyr::rename(
            Device = Cytometer,
            SuperSample = Donor,
            Sample = Sample
        )
    time <- as.POSIXct(
        c(1742229352, 1742230803, 1742232183, 1742225703, 1742227011, 1742228356)
    )

    df[["Time"]] <- as.POSIXct(time)
    df[["outcome_1"]] <- rep(rep(c("A", "B"), length.out = nrow(df) / 2), 2)
    df[["outcome_2"]] <- rep(rep(c(2.3, 5.1), length.out = nrow(df) / 2), 2)
    df[["train_validation_test"]] <- rep(c("train", "train", "test"), 2)


    ### Actual testing of cycompare
    ff_list <- c(flowCore::flowSet_to_list(fs_1), flowCore::flowSet_to_list(fs_2))
    fs <- flowCore::flowSet(ff_list)
    flowCore::markernames(fs) <- stats::setNames(
        paste0(
            "M", seq_along(flowCore::markernames(fs)),
            "-", flowCore::markernames(fs)
        ),
        flowCore::markernames(fs)
    )
    gs <- flowWorkspace::GatingSet(fs)
    flowWorkspace::gs_pop_add(
        gs,
        flowCore::rectangleGate(
            "FITC-A" = c(0, 1000),
            "PE-A" = c(0, 1000),
            filterId = "FITC_PE_gate"
        )
    )
    flowWorkspace::gs_pop_add(
        gs,
        flowCore::rectangleGate(
            "FITC-A" = c(0, 1000),
            "KrO-A" = c(-1000, 1000),
            filterId = "FITC_KrO_gate"
        ),
        parent = "FITC_PE_gate"
    )
    flowWorkspace::gs_pop_add(
        gs,
        flowCore::rectangleGate(
            "PE-A" = c(0, 1000),
            "PB-A" = c(-1000, 1000),
            filterId = "PE_PB_gate"
        ),
        parent = "FITC_PE_gate"
    )
    flowWorkspace::gs_get_pop_paths(gs)

    gs_per_sample <- sapply(names(ff_list), simplify = FALSE, function(x) {
        gs[1]
    })

    relevant_mn <- flowCore::markernames(fs)
    cn <- c("FITC-A", "PE-A", "ECD-A", "PB-A", "KrO-A")
    marker_to_gate <- list(
        "FITC-A" = c("/FITC_PE_gate"),
        # Twice the same gate possible?
        "PE-A" = c("/FITC_PE_gate/FITC_KrO_gate", "/FITC_PE_gate/FITC_KrO_gate"),
        # The following two gates dont make a lot of sense, but nice for testing
        "ECD-A" = c(
            "/FITC_PE_gate/FITC_KrO_gate",
            "/FITC_PE_gate"
        ),
        "PB-A" = c(
            "/FITC_PE_gate/PE_PB_gate",
            "/FITC_PE_gate/FITC_KrO_gate",
            "/FITC_PE_gate"
        ),
        "KrO-A" = c("/FITC_PE_gate")
    )
    return(
        list(
            "fs" = fs,
            "df" = df,
            "relevant_mn" = relevant_mn,
            "marker_to_gate" = marker_to_gate,
            "gatingsets" = gs_per_sample
        )
    )
}


prepare_testdata_fcs <- function(n_repeats = 6, tempdir) {
    tmp <- prepare_testdata()
    tmp$df[["old_FILE"]] <- tmp$df$File

    df_artificial <- tmp$df[1:3, ]
    df_artificial$Device <- "NEW_aurora"
    df_artificial$SuperSample <- "NEW_Donor_1"
    df_artificial$File <- gsub(pattern = "CyA", replacement = "CyA.NEW", df_artificial$File)
    df_artificial$train_validation_test <- "prospective"
    df_complete <- rbind(tmp$df, df_artificial)
    df_complete$Study <- c(rep("StudyA", nrow(tmp$df)), rep("Study_B", nrow(df_artificial)))

    for (set_i in seq_len(n_repeats)) {
        for (i in seq_len(nrow(df_complete))) {
            newpath <- file.path(tempdir, paste0(df_complete$File[i], "_", set_i, ".fcs"))
            flowCore::write.FCS(tmp$fs[[df_complete[["old_FILE"]][i]]], newpath)
        }
    }


    ff_files <- list.files(tempdir, recursive = TRUE, full.names = TRUE)
    ff_list <- sapply(ff_files, flowCore::read.FCS, simplify = FALSE)

    df_complete_v2 <- do.call(rbind, rep(list(df_complete), n_repeats))
    df_complete_v2$File <- ff_files

    df_complete_v2[["Device"]] <- rep(
        c("CyA", "CyB", "CyC"),
        each = 3 * n_repeats
    )
    df_complete_v2[["SuperSample"]] <- c(
        paste0("Donor_", c(1:(3 * n_repeats))),
        paste0("Donor_", c(1:(3 * n_repeats))),
        paste0("NEW_Donor_", c(1:(3 * n_repeats)))
    )
    df_complete_v2[["Sample"]] <- c(
        paste0("Sample_", 1:(3 * n_repeats * 2)),
        paste0("Sample_", 1:(3 * n_repeats))
    )

    df_complete_v2[["train_validation_test"]] <- c(
        rep(
            c(
                rep("train", n_repeats),
                rep("validation", n_repeats),
                rep("test", n_repeats)
            ),
            # repeat twice; once for each AURORA, once for FORTESSA
            2
        ),
        rep("prospective", 3 * n_repeats)
    )
    df_complete_v2[["Study"]] <- c(
        rep("BAD", 3 * n_repeats * 2),
        rep("NEW", 3 * n_repeats)
    )
    gslist <- sapply(names(ff_list), simplify = FALSE, function(x) {
        tmp$gatingsets[[1]]
    })
    return(
        list(
            "df" = df_complete_v2,
            "fs" = flowCore::flowSet(ff_list),
            "relevant_mn" = tmp$relevant_mn,
            "marker_to_gate" = tmp$marker_to_gate,
            "gatingsets" = gslist
        )
    )
}
