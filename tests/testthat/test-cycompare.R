devtools::load_all()

test_that("CyCompare Badalona samples", {
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
    df[["train_validation_test"]] <- rep(c("train", "train", "test"), 2)


    ### Actual testing of cycompare
    ff_list <- ff_files_panel |> sapply(flowCore::read.FCS, simplify = FALSE)
    ff_list_downsampled <- lapply(ff_list, function(x) {
        # x[sample(nrow(x), 1e3), ]
        x
    })
    flowCore::flowSet(ff_list_downsampled)


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

    tmp <- cycompare(
        flowframes = ff_list_downsampled,
        ff_columns_relevant = names(relevant_mn),
        df = df,
        max_events_postgate = 1e4,
        outcome_columns_df = c("outcome_1", "outcome_2"),
        outcome_models = list("glmnet" = glmnet::cv.glmnet),
        gatingsets = gslist,
        gatename_primary = "/Singlets/CD45+/CD3+",
        marker_to_gate = marker_to_gate
    )
    pdf("removeme.pdf", width = 20)
    print(tmp)
    dev.off()
})
