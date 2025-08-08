#' Prepare a small, reproducible flow cytometry test dataset
#'
#' This helper simulates two 3-sample \code{flowCore::flowSet}s measured on two
#' pseudo-devices ("CyA" and "CyB"), constructs simple rectangular gates for a
#' minimal gating hierarchy, and returns both the raw data and convenient
#' metadata for unit tests and examples.
#'
#' Imagine the scenario where you have three samples measured on both devices,
#' and the second device (CyB) has a constant shift of +5 in all channels to simulate a
#' device effect.
#'
#' @details
#' The second device (CyB) is constructed by adding a constant shift (+5) to all
#' expression channels, providing a simple between-device effect useful for
#' normalization or alignment tests.
#'
#' The returned list contains:
#' \itemize{
#'   \item \code{fs}: A \code{flowCore::flowSet} with six samples (three per device).
#'   \item \code{df}: A \code{tibble} of sample‑level metadata (Study, Device, SuperSample, Sample, Time, outcomes, split).
#'   \item \code{relevant_mn}: Character vector of marker names (prefixed as \code{"M#-<channel>"} for stability in tests).
#'   \item \code{marker_to_gate}: Named \code{list} mapping channels to example gating paths.
#'   \item \code{gatingsets}: \code{list} of per‑sample \code{flowWorkspace::GatingSet} objects
#'         (here, simple references to a shared template for speed).
#' }
#'
#' @return A named \code{list} with elements \code{fs}, \code{df}, \code{relevant_mn},
#' \code{marker_to_gate}, and \code{gatingsets}.
#'
#' @keywords testing datasets
#' @export
prepare_testdata <- function() {
    # Bindings created by tidyr/dplyr NSE; predeclare to avoid R CMD check NOTES.
    File <- Cytometer <- Donor <- Sample <- NULL

    # Define the channels to simulate. Keep short for fast examples/tests.
    cn <- c("FITC-A", "PE-A", "ECD-A", "PB-A", "KrO-A")

    # Simulate first device (CyA): 3 samples, 500 cells each.
    set.seed(42)
    # Note: A stray call to cytobench::simulate_ff in the original code had no effect and was removed.
    fs_1 <- cytobench::simulate_fs(n_samples = 3, ncells = 500, flowcore = TRUE, columns = cn)
    flowCore::sampleNames(fs_1) <- sprintf("StudyA_d%02d_s%02d_CyA", 1:3, 1:3)

    # Simulate second device (CyB) with a constant +5 shift to all channels to emulate a device effect.
    set.seed(43)
    fs_2 <- cytobench::simulate_fs(n_samples = 3, ncells = 500, flowcore = TRUE, columns = cn)
    fs_2 <- flowCore::fsApply(fs_2, function(x) {
        flowCore::exprs(x) <- flowCore::exprs(x) + 5
        x
    })
    flowCore::sampleNames(fs_2) <- sprintf("StudyA_d%02d_s%02d_CyB", 1:3, 1:3)

    # Build sample-level metadata from filenames and add outcomes + a simple train/test split.
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

    # Fixed timestamps (POSIXct) to keep tests deterministic across platforms.
    time <- as.POSIXct(
        c(1742229352, 1742230803, 1742232183, 1742225703, 1742227011, 1742228356),
        origin = "1970-01-01", tz = "UTC"
    )
    df[["Time"]] <- time

    # Add simple outcomes and a naive train/test split; adjust as your tests require.
    df[["outcome_1"]] <- rep(rep(c("A", "B"), length.out = nrow(df) / 2), 2)
    df[["outcome_2"]] <- rep(rep(c(2.3, 5.1), length.out = nrow(df) / 2), 2)
    df[["train_validation_test"]] <- rep(c("train", "train", "test"), 2)

    # Combine into a single flowSet and stabilize marker names with a synthetic prefix.
    ff_list <- c(flowCore::flowSet_to_list(fs_1), flowCore::flowSet_to_list(fs_2))
    fs <- flowCore::flowSet(ff_list)
    flowCore::markernames(fs) <- stats::setNames(
        paste0("M", seq_along(flowCore::markernames(fs)), "-", flowCore::markernames(fs)),
        flowCore::markernames(fs)
    )

    # Build a tiny gating hierarchy with rectangular gates; enough structure for downstream tests.
    gs <- flowWorkspace::GatingSet(fs)
    flowWorkspace::gs_pop_add(
        gs,
        flowCore::rectangleGate(
            "FITC-A" = c(0, 1000),
            "PE-A"   = c(0, 1000),
            filterId = "FITC_PE_gate"
        )
    )
    flowWorkspace::gs_pop_add(
        gs,
        flowCore::rectangleGate(
            "FITC-A" = c(0, 1000),
            "KrO-A"  = c(-1000, 1000),
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
    # Touch population paths to realize structure (useful in some downstream workflows).
    flowWorkspace::gs_get_pop_paths(gs)

    # For convenience, provide a per-sample list of GatingSets; reuse the template to keep runtime tiny.
    gs_per_sample <- sapply(names(ff_list), simplify = FALSE, function(x) gs[1])

    # Prepare outputs used across tests.
    relevant_mn <- flowCore::markernames(fs)
    marker_to_gate <- list(
        "FITC-A" = c("/FITC_PE_gate"),
        # Duplicate on purpose to test duplicate handling in downstream code.
        "PE-A"   = c("/FITC_PE_gate/FITC_KrO_gate", "/FITC_PE_gate/FITC_KrO_gate"),
        # Mixed paths expand coverage for path handling utilities.
        "ECD-A"  = c("/FITC_PE_gate/FITC_KrO_gate", "/FITC_PE_gate"),
        "PB-A"   = c("/FITC_PE_gate/PE_PB_gate", "/FITC_PE_gate/FITC_KrO_gate", "/FITC_PE_gate"),
        "KrO-A"  = c("/FITC_PE_gate")
    )

    # Return structure
    list(
        "fs"             = fs,
        "df"             = df,
        "relevant_mn"    = relevant_mn,
        "marker_to_gate" = marker_to_gate,
        "gatingsets"     = gs_per_sample
    )
}

#' Prepare and write artificial FCS files for testing
#'
#' This function extends the synthetic dataset from `prepare_testdata()` by
#' adding a third artificial device/sample set, writing the resulting samples to
#' `.fcs` files, and returning the updated data structures.
#'
#' It is primarily intended for unit tests or vignettes that require actual FCS
#' files on disk. The generated data include:
#' \itemize{
#'   \item multiple devices (\code{"CyA"}, \code{"CyB"}, and a \code{"CyC"}),
#'   \item multiple repetitions of each sample,
#'   \item modified sample IDs and metadata for broader test coverage.
#' }
#'
#'
#'
#'
#' @param n_repeats Integer scalar. Number of repeated copies to write for each
#'   sample (\emph{per device}). Defaults to 6.
#' @param tempdir Character scalar. Path to a writable directory where the `.fcs`
#'   files will be created.
#'
#' @return A named \code{list} with elements:
#' \describe{
#'   \item{\code{df}}{A \code{data.frame} of sample metadata, with file paths in \code{File}.}
#'   \item{\code{fs}}{A \code{flowCore::flowSet} of all written FCS files.}
#'   \item{\code{relevant_mn}}{Character vector of marker names (from [prepare_testdata()]).}
#'   \item{\code{marker_to_gate}}{Named list mapping channel names to example gating paths.}
#'   \item{\code{gatingsets}}{List of \code{flowWorkspace::GatingSet} objects (references to the template).}
#' }
#'
#' @details
#' In this example dataset, we simulate measurements from three cytometry devices:
#' **CyA**, **CyB**, and **CyC**.
#'
#' - Devices **CyA** and **CyB** belong to the same example study (*EX*), where
#'   `n_repeats * 3` donors are measured on **both** devices.
#'   For each sample, two outcomes are recorded:
#'   \itemize{
#'     \item `outcome_1`: categorical with groups `"A"` and `"B"`.
#'     \item `outcome_2`: numeric with values `2.3` or `5.1`.
#'   }
#'   The *EX* study samples are assigned to one of three splits:
#'   `"train"`, `"validation"`, or `"test"`.
#'
#' - Device **CyC** represents a *prospective* scenario: a new instrument with
#'   three unique donors, each measured `n_repeats` times.
#'   These samples are labeled `"prospective"` in the split column.
#'
#' This setup allows you to explore common analysis tasks such as:
#' \itemize{
#'   \item Comparing measurements across devices,
#'   \item Handling repeated measurements,
#'   \item Working with both categorical and numeric outcomes,
#'   \item Training and evaluating models using predefined dataset splits.
#' }
#'
#' All samples are written to `.fcs` files \emph{n\_repeats} times with a suffix
#' \code{"_<repeat_index>.fcs"}.
#'
#' @export
#' @keywords testing datasets
prepare_testdata_fcs <- function(n_repeats = 6, tempdir) {
    # Generate base synthetic dataset from prepare_testdata()
    tmp <- prepare_testdata()
    # Preserve the original file IDs for lookup in tmp$fs
    tmp$df[["old_FILE"]] <- tmp$df$File

    # Create an artificial third device using the first three samples
    df_artificial <- tmp$df[1:3, ]
    df_artificial$Device <- "CyA.NEW"
    df_artificial$SuperSample <- "NEW_Donor_1"
    df_artificial$File <- gsub(pattern = "CyA", replacement = "CyA.NEW", df_artificial$File)
    df_artificial$train_validation_test <- "prospective"

    # Combine original + artificial metadata
    df_complete <- rbind(tmp$df, df_artificial)
    df_complete$Study <- c(rep("StudyA", nrow(tmp$df)), rep("Study_B", nrow(df_artificial)))

    # Write each sample to disk n_repeats times
    for (set_i in seq_len(n_repeats)) {
        for (i in seq_len(nrow(df_complete))) {
            newpath <- file.path(tempdir, paste0(df_complete$File[i], "_", set_i, ".fcs"))
            flowCore::write.FCS(tmp$fs[[df_complete[["old_FILE"]][i]]], newpath)
        }
    }

    # Read the generated FCS files back in
    ff_files <- list.files(tempdir, recursive = TRUE, full.names = TRUE)
    ff_list <- sapply(ff_files, flowCore::read.FCS, simplify = FALSE)

    # Replicate metadata for each repeat and update File paths
    df_complete_v2 <- do.call(rbind, rep(list(df_complete), n_repeats))
    df_complete_v2$File <- ff_files

    # Overwrite device and donor IDs to create variety for testing
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

    # Assign train/validation/test/prospective splits
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

    # Assign study labels for each block
    df_complete_v2[["Study"]] <- c(
        rep("EX", 3 * n_repeats * 2),
        rep("NEW", 3 * n_repeats)
    )

    # Build per-sample GatingSet list (reuse the first template for simplicity)
    gslist <- sapply(names(ff_list), simplify = FALSE, function(x) {
        tmp$gatingsets[[1]]
    })

    # Return results
    list(
        "df" = df_complete_v2,
        "fs" = flowCore::flowSet(ff_list),
        "relevant_mn" = tmp$relevant_mn,
        "marker_to_gate" = tmp$marker_to_gate,
        "gatingsets" = gslist
    )
}
