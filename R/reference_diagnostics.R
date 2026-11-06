#' Plot stable-reference profiles before and after normalization
#'
#' Gray violins show every protein observed in all samples. Blue lines, points,
#' and boxplots show the selected reference proteins. Each protein is centered
#' on its own mean across samples, separately before and after normalization.
#' Both panels use the same y-axis limits. Selected references are part of the
#' background distribution, so this is a descriptive check rather than an
#' independent validation of biological invariance.
#'
#' @param se A SummarizedExperiment returned by a stable-reference
#'   normalization function.
#' @param method `"genorm"` or `"pca"`.
#' @param style `"BPG"` (default) or `"ggplot2"`.
#' @param ylimits Shared y-axis limits for both panels; default `c(-2, 2)`.
#' @return A BPG-style trellis plot or ggplot object.
#' @export
plot_stable_reference_diagnostics <- function(se,
    method = c("genorm", "pca"), style = c("BPG", "ggplot2"),
    ylimits = c(-2, 2)) {
    if (!methods::is(se, "SummarizedExperiment"))
        stop("se must be a SummarizedExperiment", call. = FALSE)
    method <- match.arg(method)
    style <- match.arg(style)
    if (!is.numeric(ylimits) || length(ylimits) != 2L ||
        any(!is.finite(ylimits)) || ylimits[[1L]] >= ylimits[[2L]])
        stop("ylimits must contain increasing finite bounds", call. = FALSE)
    info <- S4Vectors::metadata(se)[[paste0("reference_normalization_", method)]]
    if (is.null(info) || !all(c("selected", "input_assay", "output_assay") %in%
        names(info)))
        stop("Run the selected stable-reference normalization method first",
            call. = FALSE)
    if (!all(c(info$input_assay, info$output_assay) %in%
        SummarizedExperiment::assayNames(se)) ||
        !all(info$selected %in% rownames(se)))
        stop("Reference normalization provenance does not match se", call. = FALSE)
    raw <- SummarizedExperiment::assay(se, info$input_assay)
    adjusted <- SummarizedExperiment::assay(se, info$output_assay)
    complete <- rowSums(is.finite(raw) & raw > 0) == ncol(raw)
    ids <- rownames(se)[complete]
    if (!length(ids) || !all(info$selected %in% ids))
        stop("Selected references must be observed in every sample", call. = FALSE)
    labels <- colnames(se)
    if ("sample_id" %in% names(SummarizedExperiment::colData(se))) {
        short <- as.character(SummarizedExperiment::colData(se)$sample_id)
        if (all(!is.na(short) & nzchar(short)) && !anyDuplicated(short))
            labels <- short
    }
    profiles <- do.call(rbind, lapply(c("Before", "After"), function(stage) {
        values <- if (stage == "Before") raw else adjusted
        logged <- log2(values[ids, , drop = FALSE])
        if (any(!is.finite(logged)))
            stop("Complete proteins have nonfinite normalized values", call. = FALSE)
        deviation <- sweep(logged, 1L, rowMeans(logged), "-")
        data.frame(stage = stage, protein = rep(ids, times = ncol(se)),
            sample_index = rep(seq_len(ncol(se)), each = length(ids)),
            sample_name = rep(labels, each = length(ids)),
            deviation = as.vector(deviation),
            selected = rep(ids %in% info$selected, times = ncol(se)))
    }))
    profiles$stage <- factor(profiles$stage, levels = c("Before", "After"))
    profiles$sample_name <- factor(profiles$sample_name, levels = labels)
    references <- profiles[profiles$selected, , drop = FALSE]
    title <- if (method == "genorm") "geNorm" else "PCA"
    y_label <- "Centered log2 abundance"
    if (style == "ggplot2") {
        return(ggplot2::ggplot() +
            ggplot2::geom_violin(data = profiles,
                ggplot2::aes(x = .data$sample_name, y = .data$deviation,
                    group = .data$sample_name), fill = "grey90",
                colour = "grey70", width = 0.85, scale = "width") +
            ggplot2::geom_line(data = references,
                ggplot2::aes(x = .data$sample_name, y = .data$deviation,
                    group = .data$protein), colour = "#287D9C",
                alpha = 0.12, linewidth = 0.25) +
            ggplot2::geom_boxplot(data = references,
                ggplot2::aes(x = .data$sample_name, y = .data$deviation,
                    group = .data$sample_name), fill = "white",
                colour = "#176B88", width = 0.32,
                outlier.shape = NA, alpha = 0.8) +
            ggplot2::geom_point(data = references,
                ggplot2::aes(x = .data$sample_name, y = .data$deviation),
                colour = "#176B88", alpha = 0.18, size = 0.45) +
            ggplot2::facet_wrap(ggplot2::vars(.data$stage), nrow = 1) +
            ggplot2::coord_cartesian(ylim = ylimits) +
            ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(angle = 45)) +
            ggplot2::labs(title = title, x = NULL, y = y_label) +
            ggplot2::theme_classic(base_size = 12) +
            ggplot2::theme(
                panel.border = ggplot2::element_rect(colour = "black",
                    fill = NA, linewidth = 0.7),
                axis.line = ggplot2::element_blank(),
                axis.title = ggplot2::element_text(face = "bold"),
                axis.text = ggplot2::element_text(colour = "black", face = "bold"),
                plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
                strip.background = ggplot2::element_rect(fill = "white",
                    colour = "black"),
                strip.text = ggplot2::element_text(face = "bold"),
                panel.grid = ggplot2::element_blank()))
    }
    n_samples <- ncol(se)
    lattice::xyplot(deviation ~ sample_index | stage, data = profiles,
        layout = c(2, 1), xlim = c(0.5, n_samples + 0.5), ylim = ylimits,
        xlab = NULL, ylab = y_label, main = title,
        main.just = "center",
        scales = list(x = list(at = seq_len(n_samples), labels = labels,
                rot = 45, cex = 0.8, tck = c(0.5, 0), alternating = FALSE,
                relation = "free"),
            y = list(cex = 0.8, tck = c(0.5, 0),
                alternating = FALSE)),
        par.settings = list(
            axis.line = list(col = "black", lwd = 1),
            strip.background = list(col = "white"),
            par.xlab.text = list(font = 2),
            par.ylab.text = list(font = 2)),
        panel = function(x, y, subscripts, ...) {
            panel_data <- profiles[subscripts, , drop = FALSE]
            for (sample in seq_len(n_samples)) {
                values <- panel_data$deviation[panel_data$sample_index == sample]
                density <- stats::density(values, na.rm = TRUE)
                width <- density$y / max(density$y) * 0.42
                lattice::panel.polygon(
                    c(sample - width, rev(sample + width)),
                    c(density$x, rev(density$x)),
                    col = "grey90", border = "grey70")
            }
            selected <- panel_data[panel_data$selected, , drop = FALSE]
            for (protein in unique(selected$protein)) {
                rows <- selected[selected$protein == protein, , drop = FALSE]
                rows <- rows[order(rows$sample_index), , drop = FALSE]
                lattice::panel.lines(rows$sample_index, rows$deviation,
                    col = grDevices::adjustcolor("#287D9C", alpha.f = 0.12),
                    lwd = 0.5)
            }
            for (sample in seq_len(n_samples)) {
                values <- selected$deviation[selected$sample_index == sample]
                q <- stats::quantile(values, c(0.25, 0.5, 0.75))
                lattice::panel.rect(sample - 0.16, q[[1L]], sample + 0.16,
                    q[[3L]], col = "white", border = "#176B88", lwd = 1.5)
                lattice::panel.segments(sample - 0.16, q[[2L]],
                    sample + 0.16, q[[2L]], col = "#176B88", lwd = 1.5)
                lattice::panel.points(rep(sample, length(values)), values,
                    pch = 19, cex = 0.25,
                    col = grDevices::adjustcolor("#176B88", alpha.f = 0.18))
            }
        })
}
