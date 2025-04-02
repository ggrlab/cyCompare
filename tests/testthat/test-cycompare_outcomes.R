devtools::load_all()

test_that("CyCompare Badalona samples, outcomes", {
    options(warn = 1)
    ## Preparation of a proper dataset
    ff_files <- list.files("example_data/2025-03-17_BAD", recursive = TRUE, full.names = TRUE)
    ff_files_panel <- grep("panel", ff_files, value = TRUE)
    df <- tibble::tibble(
        File = ff_files_panel
    ) |>
        dplyr::mutate(
            Device = basename(dirname(dirname(File))),
            # You could have multiple SAMPLES for each SUPERSAMPLE
            SuperSample = rep(paste0("Donor_", c(1:3)), 2),
            Sample = rep(paste0("Sample_", c(1:3)), 2)
        )
    time <- lapply(ff_files_panel, function(x) {
        fcs_head <- flowCore::read.FCSheader(x)[[1]]
        begin_time <- fcs_head[["$BTIM"]]
        begin_date <- fcs_head[["$DATE"]]
        as.POSIXct(paste0(begin_date, "_", begin_time), format = "%d-%b-%Y_%H:%M:%OS")
    }) |> unlist()

    df[["Time"]] <- as.POSIXct(time)
    df[["outcome_1"]] <- rep(rep(c("A", "B"), length.out = nrow(df) / 2), 2)
    df[["outcome_2"]] <- rep(rep(c(2.3, 5.1), length.out = nrow(df) / 2), 2)
    df[["train_validation_test"]] <- rep(c("train", "validation", "test"), 2)
    df_artificial <- df[1:3, ]
    df_artificial$Device <- "NEW_aurora"
    df_artificial$SuperSample <- "NEW_Donor_1"
    df_artificial$File <- gsub(pattern = "/2025-03-17_BAD/aurora/", replacement = "/NEW_aurora/", df_artificial$File)
    lapply(dirname(df_artificial$File), dir.create, recursive = TRUE, showWarnings = FALSE)
    file.copy(df[1:3, ]$File, df_artificial$File)
    df_artificial$train_validation_test <- "prospective"
    df_complete <- rbind(df, df_artificial)
    df_complete$Study <- c(rep("BAD", nrow(df)), rep("NEW", nrow(df_artificial)))
    for (x in df_complete$File) {
        newpath <- sub("example_data", "example_data/small", x)
        dir.create(dirname(newpath), recursive = TRUE, showWarnings = FALSE)
        file.copy(x, file.path(newpath))
    }
    n_times <- 6
    for (set_i in 2:n_times) {
        for (x in df_complete$File) {
            newpath <- sub("example_data", "example_data/small", x)
            newpath <- gsub("(BAD[0-9]+)", paste0("\\1_v", set_i), newpath)
            dir.create(dirname(newpath), recursive = TRUE, showWarnings = FALSE)
            file.copy(x, file.path(newpath))
        }
    }

    ff_files <- list.files("example_data/small", recursive = TRUE, full.names = TRUE)
    lapply(ff_files, function(x) {
        a <- flowCore::read.FCS(x)
        flowCore::write.FCS(a[sample(nrow(a), min(1e3, nrow(a))), ], x)
    })

    ### Actual testing of cycompare
    ff_files <- list.files("example_data/small", recursive = TRUE, full.names = TRUE)
    ff_files_panel <- grep("panel", ff_files, value = TRUE)
    ff_list <- ff_files_panel |> sapply(flowCore::read.FCS, simplify = FALSE)
    ff_list_downsampled <- lapply(ff_list, function(x) {
        # x[sample(nrow(x), 1e3), ]
        x
    })
    flowCore::flowSet(ff_list_downsampled)
    df_complete_v2 <- do.call(rbind, rep(list(df_complete), n_times))
    df_complete_v2$File <- ff_files
    df_complete_v2 <- dplyr::mutate(
        df_complete_v2,
        Device = basename(dirname(dirname(File))),
    )
    df_complete_v2[["SuperSample"]] <- c(
        paste0("Donor_", c(1:(3 * n_times))),
        paste0("Donor_", c(1:(3 * n_times))),
        paste0("NEW_Donor_", c(1:(3 * n_times)))
    )
    df_complete_v2[["Sample"]] <- c(
        paste0("Sample_", 1:(3 * n_times * 2)),
        paste0("Sample_", 1:(3 * n_times))
    )

    df_complete_v2[["train_validation_test"]] <- c(
        rep(
            c(
                rep("train", n_times),
                rep("validation", n_times),
                rep("test", n_times)
            ),
            # repeat twice; once for each AURORA, once for FORTESSA
            2
        ),
        rep("prospective", 3 * n_times)
    )
    df_complete_v2[["Study"]] <- c(
        rep("BAD", 3 * n_times * 2),
        rep("NEW", 3 * n_times)
    )
    # with(df_complete_v2, table(outcome_1, train_validation_test, Study))

    prot_aurora <- cytoKal::read_protocol("example_data/2025-03-17_BAD/2025-03-18_aurora_renamed_TcellFull.protocol")
    prot_fortessa <- cytoKal::read_protocol("example_data/2025-03-17_BAD/2025-03-19_fortessa_BADALONA_renamed_TcellFull.protocol")

    gslist <- sapply(names(ff_list_downsampled), simplify = FALSE, function(x) {
        if (grepl("aurora", x)) {
            prot_aurora[["gatingsets"]][[1]]
        } else {
            prot_fortessa[["gatingsets"]][[1]]
        }
    })
    # flowWorkspace::gh_get_pop_paths(gslist[[4]])

    relevant_mn <- flowCore::flowSet(ff_list_downsampled) |> flowCore::markernames()
    marker_to_gate <- list(
        "CD45RA-FITC-A" = c("/Singlets/CD45+/CD3+/CD45RA+CCR7-"),
        "CCR7-PE-A" = c("/Singlets/CD45+/CD3+/CD45RA+CCR7+"),
        "CD28-ECD-A" = c(
            "/Singlets/CD45+/CD3+/CD27+CD28+",
            "/Singlets/CD45+/CD3+/CD27+CD28+"
        ),
        "PD1-PC5.5-A" = c("/Singlets/CD45+/CD3+/CD57-PD1+", "/Singlets/CD45+/CD3+/CD57-PD1+"),
        "CD27-PC7-A" = c(
            "/Singlets/CD45+/CD3+/CD27+CD28+",
            "/Singlets/CD45+/CD3+/CD27+CD28+"
        ),
        "CD4-APC-A" = "/Singlets/CD45+/CD3+/CD4+",
        "CD8-AF700-A" = "/Singlets/CD45+/CD3+/CD8+",
        "CD3-AA750-A" = "/Singlets/CD45+/CD3+",
        "CD57-PB-A" = c(
            "/Singlets/CD45+/CD3+/CD57+PD1+",
            "/Singlets/CD45+/CD3+/CD57+PD1-",
            "/Singlets/CD45+/CD3+/CD57+PD1+",
            "/Singlets/CD45+/CD3+/CD57+PD1-"
        ),
        "CD45-KrO-A" = c("/Singlets/CD45+")
    )


    # options(
    #     future.globals.maxSize = 10 * 1024^3
    # )
    pacman::p_load("mlr3learners")

    result_analysis <- cycompare_outcomes_analyse(
        flowframes = ff_list_downsampled,
        ff_columns_relevant = names(relevant_mn),
        df = df_complete_v2,
        dfcol_outcomes = c("outcome_1", "outcome_2"),
        outcome_models = list("glmnet" = glmnet::cv.glmnet),
        n_events_postgate = 1e4,
        gatingsets = gslist,
        gatename_primary = "/Singlets/CD45+/CD3+",
        marker_to_gate = marker_to_gate,
        dfcol_train_validation_other = "train_validation_test",
        clustering_outdir = "intermediate/clustering"
    )
    browser()
    pdf("removeme.pdf", width = 20)
    print(tmp)
    dev.off()
})
