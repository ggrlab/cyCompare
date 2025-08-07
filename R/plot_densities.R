#' Plot Density Distributions for Gated Flow Cytometry Data
#'
#' This function generates density plots for multiple flow cytometry markers across different devices and samples.
#' It applies optional transformations, merges metadata, and creates a faceted plot of density curves.
#'
#' Every (unique!) flowframe must have a unique name present in `df$File`, which will be
#' matched to the `"Sample"` column. One sample would often (not necessarily!)
#' correspond to multiple flowframes, e.g., if multiple files were acquired
#' for the same sample on multiple devices.
#'
#' @param ff_gated
#' Named list of `flowFrame` objects, each representing gated flow cytometry data for one sample.
#' Names must correspond to entries in `df$File`. Can be a flowset or a named list of flowFrames.
#' @param df
#' A `data.table` containing metadata for samples. Must include `"File"`, `dfcol_grouping_samples`, and `"Sample"` columns.
#' @param device_colors
#' Optional named vector of colors for each device (e.g., `c("Fortessa" = "blue", "Aurora" = "red")`).
#' @param transformlist Either a single transformation function or a named list of functions, one per marker.
#'                      If NULL, no transformation is applied.
#' @param density_n Integer; number of points for kernel density estimation (default: 500).
#' @param dfcol_grouping_samples
#' Column name (character) in `df` used for device/sample grouping (default: `"Device"`).
#' @param relevant_columns
#' Character vector of markers to include in the plot. If NULL, all markers are included.
#' @param limit_density_quantile
#' Numeric between 0 and 1; limits plotted densities to a quantile threshold (e.g., 0.95).
#' Useful to suppress extreme density peaks which make all other densities barely visible.
#' Set to `NA` to disable.
#'
#' @return A `ggplot2` object showing density curves across markers, devices, and samples.
#' @export
#' @examples
#' fs <- cytobench::simulate_fs(n_samples = 3, ncells = 250, columns = c("CD4", "CD8"))
#'
#' # Metadata
#' df_meta <- data.table::data.table(
#'     File = flowCore::sampleNames(fs),
#'     Device = c("A", "B", "A"),
#'     Sample = c("S1", "S2", "S2")
#' )
#'
#' # Colors for devices
#' device_cols <- c("A" = "steelblue", "B" = "firebrick")
#'
#' # Apply density plot function
#' plot <- plot_densities(
#'     ff_gated = fs,
#'     df = df_meta,
#'     device_colors = device_cols,
#'     transformlist = function(x) asinh(x / 1000),
#'     relevant_columns = c("CD4", "CD8")
#' )
#'
#' # Show plot
#' print(plot)
plot_densities <- function(
    ff_gated,
    df,
    device_colors,
    transformlist = NULL,
    density_n = 500,
    dfcol_grouping_samples = "Device",
    relevant_columns = NULL,
    limit_density_quantile = NA) {
    . <- variable <- value <- File <- y <- V1 <- NULL # R CMD check compatibility
    if ("flowSet" %in% class(ff_gated)) {
        ff_gated <- flowCore::flowSet_to_list(ff_gated)
    }
    # Ensure ff_gated is a named list
    if (!all(!is.null(names(ff_gated)))) {
        stop("ff_gated must be a named list of flowFrame objects")
    }
    if (all(is.null(relevant_columns))) {
        relevant_columns <- flowCore::colnames(ff_gated[[1]])
    }

    # Use all columns (channels/markers) if no subset is given
    if (all(is.null(relevant_columns))) {
        relevant_columns <- flowCore::colnames(ff_gated[[1]])
    }

    # Convert single transform function to named list if needed
    transformlist <- transformlist_named(transformlist, relevant_columns)
    # Apply identity function to unused columns to keep them unchanged
    for (col_x in flowCore::colnames(ff_gated[[1]])[!flowCore::colnames(ff_gated[[1]]) %in% relevant_columns]) {
        transformlist[[col_x]] <- identity
    }


    # Extract expression data and melt into long format
    gated_dt <- lapply(ff_gated, function(x) {
        flowCore::exprs(x) |> data.table::as.data.table()
    }) |>
        data.table::rbindlist(idcol = "File", fill = TRUE) |>
        data.table::melt(id.vars = "File")

    # Keep only relevant markers
    gated_dt <- gated_dt[variable %in% relevant_columns]

    # Density estimation helper
    compute_density <- function(x, transform_x) {
        x_noNA <- x[!is.na(x)]
        if (length(x_noNA) == 0) {
            res <- data.table::data.table(x = NA_real_, y = NA_real_)
        } else {
            d <- density(transform_x(x_noNA), n = density_n)
            res <- data.table::data.table(x = d$x, y = d$y)
        }
        return(res)
    }

    # Compute densities for each File × Marker
    densities <- gated_dt[,
        {
            compute_density(value, transformlist[[as.character(variable[[1]])]])
        },
        by = .(File, variable)
    ]

    # Optionally clamp densities to avoid extreme peaks
    if (!is.na(limit_density_quantile)) {
        # Limit the density to a certain quantile within each variable
        # calculate the top quantile for each variable
        top_quantile <- densities[, quantile(y, limit_density_quantile), by = variable]
        densities <- densities[top_quantile, on = "variable"][, y := pmin(y, V1)][, V1 := NULL]
    }

    # Merge metadata: File, Device, Sample
    df_part <- data.table::data.table(df)
    df_part <- df_part[, c("File", dfcol_grouping_samples[[1]], "Sample"), with = FALSE]
    densities <- densities[df_part, on = "File"]

    # Plot densities by Sample and Marker (x = intensity, y = density)
    p0 <- ggplot2::ggplot(
        densities,
        ggplot2::aes(x = x, y = y, color = !!rlang::sym(dfcol_grouping_samples[[1]]))
    ) +
        ggplot2::geom_line() +
        ggplot2::geom_ribbon(
            ggplot2::aes(ymin = 0, ymax = y, fill = !!rlang::sym(dfcol_grouping_samples[[1]])),
            alpha = 0.2
        ) +
        ggh4x::facet_grid2(Sample ~ variable, scales = "free", independent = "y") +
        ggpubr::theme_pubr() +
        ggplot2::theme(
            axis.line.y = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank()
        ) +
        ggplot2::ylab("Density") +
        ggplot2::xlab("Transformed MFI")

    # Add manual colors if provided
    if (!all(is.null(device_colors))) {
        p0 <- p0 +
            ggplot2::scale_color_manual(values = device_colors) +
            ggplot2::scale_fill_manual(values = device_colors)
    }

    return(p0)
}
