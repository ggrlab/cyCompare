devtools::load_all()

test_that("CyCompare Badalona samples, outcomes", {
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
    for (x in df_complete$File) {
        newpath <- sub("example_data", "example_data/small", x)
        newpath <- gsub("(BAD[0-9]+)", "\\1_v2", newpath)
        dir.create(dirname(newpath), recursive = TRUE, showWarnings = FALSE)
        file.copy(x, file.path(newpath))
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
    df_complete_v2 <- rbind(df_complete, df_complete)
    df_complete_v2$File <- ff_files
    df_complete_v2 <- dplyr::mutate(
        df_complete_v2,
        Device = basename(dirname(dirname(File))),
    )
    df_complete_v2[["SuperSample"]] <- c(
        paste0("Donor_", c(1:6)),
        paste0("Donor_", c(1:6)),
        paste0("NEW_Donor_", c(1:3)),
        paste0("NEW_Donor_", c(1:3))
    )
    df_complete_v2[["Sample"]] <- c(
        paste0("Sample_", 1:12), 
        paste0("Sample_", 1:6)
    )
    
    df_complete_v2[["train_validation_test"]] <- c(
        c("train", "train", "validation", "validation", "test", "test"),
        c("train", "train", "validation", "validation", "test", "test"),
        rep("prospective", 6)
    )
    df_complete_v2[["Study"]] <- c(
        rep("BAD", 12),
        rep("NEW", 6)
    )

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


    tmp <- cycompare_outcomes(
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
    pdf("removeme.pdf", width = 20)
    print(tmp)
    dev.off()
})
