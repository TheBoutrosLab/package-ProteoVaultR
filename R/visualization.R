#' Plot the abundance of one protein group across samples
#'
#' Selects an exact `source_group_id` from the protein quantity table. The
#' default `matrix` method uses DIA-NN normalized protein-group quantities.
#' Attached sample metadata supplies shorter x axis labels when `sample_id`
#' is present. Missing quantities have no bar; zero quantities remain zero.
#'
#' @inheritParams tbl_protein_quant
#' @param protein Exact `source_group_id` of the protein group to plot, such as
#'   `"sp|O75083|WDR1_HUMAN"`.
#' @param method Protein-group quantification method; defaults to `"matrix"`.
#' @param samples Optional original sample names to select and order.
#' @param color_by Optional column name in the attached sample metadata used
#'   to color bars. Every selected sample must have a nonmissing value.
#' @param style Plotting style: `"BPG"` (default) or `"ggplot2"`.
#' @return A BPG trellis plot or ggplot object, according to `style`.
#' @export
plot_protein_abundance <- function(
    vault,
    protein,
    source_id = NULL,
    method = "matrix",
    samples = NULL,
    color_by = NULL,
    style = c("BPG", "ggplot2")
) {
    style <- match.arg(style)
    if (!is.character(protein) || length(protein) != 1L ||
        is.na(protein) || !nzchar(protein)) {
        stop("protein must be one nonempty source_group_id", call. = FALSE)
    }
    query <- tbl_protein_quant(vault = vault, source_id = source_id, method = method)
    available <- query %>%
        dplyr::distinct(.data$sample_name) %>%
        dplyr::collect()
    if (is.null(samples)) samples <- sort(available$sample_name)
    if (!is.character(samples) || !length(samples) || anyNA(samples) ||
        anyDuplicated(samples) || !all(samples %in% available$sample_name)) {
        stop("samples must be distinct names present in the quantification source", call. = FALSE)
    }
    quantities <- query %>%
        dplyr::filter(
            .data$source_group_id == .env$protein,
            .data$sample_name %in% .env$samples
        ) %>%
        dplyr::select(dplyr::all_of(c("sample_name", "quantity"))) %>%
        dplyr::collect()
    if (!nrow(quantities) || all(is.na(quantities$quantity))) {
        stop("No quantities found for this protein and sample selection", call. = FALSE)
    }
    quantities <- quantities %>%
        dplyr::group_by(.data$sample_name) %>%
        dplyr::summarise(
            quantity = aggregate_duplicate_quantities(.data$quantity, "best"),
            .groups = "drop"
        )
    quantities$sample_name <- factor(quantities$sample_name, levels = samples)
    labels <- stats::setNames(samples, samples)
    short_labels <- vault_sample_labels(vault)
    if (!is.null(short_labels)) labels[names(short_labels)] <- short_labels
    quantities$sample <- factor(
        unname(labels[as.character(quantities$sample_name)]),
        levels = rev(unname(labels[samples]))
    )
    colours <- "#4E9DA6"
    if (!is.null(color_by)) {
        metadata <- sample_metadata(vault)
        if (!is.character(color_by) || length(color_by) != 1L ||
            is.na(color_by) || !nzchar(color_by) || is.null(metadata) ||
            !color_by %in% names(metadata) || color_by == "sample_name") {
            stop("color_by must name a column in attached sample metadata", call. = FALSE)
        }
        groups <- metadata[[color_by]][match(samples, metadata$sample_name)]
        if (!(is.character(groups) || is.factor(groups)) ||
            anyNA(groups) || any(!nzchar(as.character(groups)))) {
            stop("color_by must have a nonmissing categorical value for every selected sample", call. = FALSE)
        }
        group_names <- unique(as.character(groups))
        colours <- stats::setNames(
            c("#4E9DA6", "#D98C70", "#6B9F59", "#596B82")[
                rep(seq_len(4L), length.out = length(group_names))
            ],
            group_names
        )
        if (length(group_names) > 4L) {
            colours <- stats::setNames(
                grDevices::hcl.colors(length(group_names), palette = "Dark 3"),
                group_names
            )
        }
        quantities$group <- factor(
            as.character(groups[match(as.character(quantities$sample_name), samples)]),
            levels = group_names
        )
    }
    value_range <- range(c(0, quantities$quantity), na.rm = TRUE)
    padding <- max(diff(value_range) * 0.05, 0.01)
    xlimits <- value_range + c(-padding, padding)
    if (style == "BPG") {
        plot_data <- quantities[!is.na(quantities$quantity), , drop = FALSE]
        return(BoutrosLab.plotting.general::create.barplot(
            formula = sample ~ quantity,
            data = plot_data,
            groups = if (is.null(color_by)) NULL else group,
            plot.horizontal = TRUE,
            col = unname(colours),
            border.col = "black",
            border.lwd = 1,
            key = if (is.null(color_by)) list(text = list(lab = "")) else list(
                rectangles = list(col = rev(unname(colours))),
                text = list(lab = rev(names(colours))),
                columns = 1,
                space = "right"
            ),
            xlimits = xlimits,
            xlab.label = "Protein-group abundance",
            ylab.label = NULL,
            xlab.cex = 1,
            xaxis.cex = 0.8,
            yaxis.cex = 0.8,
            left.padding = 3,
            right.padding = 0.5,
            xaxis.tck = c(0.5, 0),
            yaxis.tck = c(0.5, 0)
        ))
    }
    plot <- ggplot2::ggplot(quantities,
        ggplot2::aes(x = .data$sample, y = .data$quantity))
    if (is.null(color_by)) {
        plot <- plot + ggplot2::geom_col(width = 0.7, fill = colours,
            colour = "black", linewidth = 0.7, na.rm = TRUE)
    } else {
        plot <- plot + ggplot2::geom_col(ggplot2::aes(fill = group), width = 0.7,
            colour = "black", linewidth = 0.7, na.rm = TRUE) +
            ggplot2::scale_fill_manual(values = colours, name = color_by)
    }
    plot +
        ggplot2::scale_x_discrete(
            limits = rev(unname(labels[samples])),
            expand = ggplot2::expansion(add = c(0.5, 0.5))
        ) +
        ggplot2::coord_flip(ylim = xlimits) +
        ggplot2::labs(x = NULL, y = "Protein-group abundance") +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::theme(
            panel.grid = ggplot2::element_blank(),
            panel.border = ggplot2::element_rect(
                colour = "black", fill = NA, linewidth = 0.5
            ),
            axis.title = ggplot2::element_text(face = "bold"),
            axis.text = ggplot2::element_text(colour = "black", face = "bold")
        )
}

#' Plot differential abundance as a volcano plot
#'
#' Plots log2 fold change against `-log10` of a limma q-value (BH-adjusted
#' p-value) or p-value. The y-axis displays the transformed values, such as
#' 0, 1, and 2. When showing p-values, `color_by` can use q-values for
#' significance coloring. Points at or below `cutoff` for the color metric are
#' colored by negative versus nonnegative log2 fold change; other points are gray. Missing
#' statistics are omitted. An exact zero probability is plotted one log10 unit
#' above the smallest positive probability in the contrast, or at machine
#' precision when all plotted probabilities are zero. Annotated points are
#' highlighted in orange, drawn over other points, and connected to their labels.
#'
#' Protein labels can match a feature ID, complete source protein-group name,
#' group member accession or name, or a mapped identifier. A term matching
#' several protein groups labels each matched point. Unmatched terms error.
#' `label_by` selects a character column from the fitted result's rowData
#' annotations to display in place of the matched accession.
#'
#' @param x An object returned by [fit_limma()].
#' @param contrast Name of the contrast to plot. May be omitted only if the fit
#'   contains one contrast.
#' @param significance `"q_value"` (default) for `adj.P.Val`, or `"p_value"`
#'   for `P.Value`.
#' @param color_by Probability used for point colors and the legend:
#'   `"q_value"` or `"p_value"`. Defaults to `significance`. When
#'   `significance = "q_value"`, only `"q_value"` is allowed.
#' @param proteins Optional character vector of protein identifiers to label.
#' @param label_by Optional rowData column name to use as the displayed label,
#'   for example `"leading_gene_name"`. The selected values must be present.
#' @param cutoff Inclusive threshold for `color_by`, default 0.05.
#' @param style `"BPG"` (default) or `"ggplot2"`.
#' @return A BPG trellis plot or a ggplot object, according to `style`.
#' @examples
#' \dontrun{
#' plot_limma_volcano(fit, contrast = "SBI_vs_DMSO",
#'   proteins = "O75083", label_by = "leading_gene_name")
#' plot_limma_volcano(fit, contrast = "SBI_vs_DMSO",
#'   significance = "p_value", color_by = "q_value", cutoff = 0.05)
#' }
#' @export
plot_limma_volcano <- function(
    x,
    contrast = NULL,
    significance = c("q_value", "p_value"),
    color_by = NULL,
    proteins = NULL,
    label_by = NULL,
    cutoff = 0.05,
    style = c("BPG", "ggplot2")
) {
    if (!inherits(x, "proteovault_limma") || !is.data.frame(x$results)) {
        stop("x must be a proteovault_limma result", call. = FALSE)
    }
    significance <- match.arg(significance)
    if (is.null(color_by)) {
        color_by <- significance
    }
    color_by <- match.arg(color_by, c("q_value", "p_value"))
    if (significance == "q_value" && color_by != "q_value") {
        stop("color_by must be q_value when significance is q_value",
            call. = FALSE)
    }
    style <- match.arg(style)
    if (!is.numeric(cutoff) || length(cutoff) != 1L ||
        !is.finite(cutoff) || cutoff <= 0 || cutoff >= 1) {
        stop("cutoff must be one number between 0 and 1", call. = FALSE)
    }
    if (!is.null(proteins) &&
        (!is.character(proteins) || anyNA(proteins) ||
            any(!nzchar(proteins)) || anyDuplicated(proteins))) {
        stop("proteins must contain distinct, nonempty identifiers", call. = FALSE)
    }
    if (!is.null(label_by) &&
        (!is.character(label_by) || length(label_by) != 1L ||
            is.na(label_by) || !nzchar(label_by))) {
        stop("label_by must name one rowData column", call. = FALSE)
    }
    available <- unique(as.character(x$results$contrast))
    if (is.null(contrast)) {
        if (length(available) != 1L) {
            stop("Select one contrast from: ", paste(available, collapse = ", "),
                call. = FALSE)
        }
        contrast <- available
    }
    if (!is.character(contrast) || length(contrast) != 1L ||
        is.na(contrast) || !contrast %in% available) {
        stop("contrast must name one fitted comparison: ",
            paste(available, collapse = ", "), call. = FALSE)
    }
    data <- as.data.frame(x$results[x$results$contrast == contrast, , drop = FALSE])
    probability_column <- if (significance == "q_value") "adj.P.Val" else "P.Value"
    probability <- data[[probability_column]]
    if (!is.numeric(probability) || any(probability < 0 | probability > 1,
        na.rm = TRUE)) {
        stop("Selected probabilities must be numeric values between 0 and 1",
            call. = FALSE)
    }
    data$probability <- probability
    color_column <- if (color_by == "q_value") "adj.P.Val" else "P.Value"
    color_probability <- data[[color_column]]
    if (!is.numeric(color_probability) ||
        any(color_probability < 0 | color_probability > 1, na.rm = TRUE)) {
        stop("Color probabilities must be numeric values between 0 and 1",
            call. = FALSE)
    }
    data$color_probability <- color_probability
    data <- data[is.finite(data$logFC) & !is.na(data$probability) &
        !is.na(data$color_probability), , drop = FALSE]
    if (!nrow(data)) {
        stop("No finite fold changes and probabilities for this contrast",
            call. = FALSE)
    }
    positive <- data$probability[data$probability > 0]
    floor_value <- if (length(positive)) min(positive) / 10 else .Machine$double.xmin
    floor_value <- max(floor_value, .Machine$double.xmin)
    data$negative_log10 <- -log10(pmax(data$probability, floor_value))
    cutoff_text <- format(cutoff, trim = TRUE, scientific = FALSE)
    legend_labels <- c(paste0("> ", cutoff_text),
        paste0("\u2264 ", cutoff_text, "; FC < 0"),
        paste0("\u2264 ", cutoff_text, "; FC \u2265 0"))
    legend_title <- if (color_by == "q_value") "q-value" else "p-value"
    data$significant <- factor(
        ifelse(data$color_probability > cutoff, legend_labels[[1L]],
            ifelse(data$logFC < 0, legend_labels[[2L]], legend_labels[[3L]])),
        levels = legend_labels
    )
    data$label <- volcano_protein_labels(data, proteins, label_by)
    annotated <- data[nzchar(data$label), , drop = FALSE]
    maximum <- max(data$negative_log10)
    top_tick <- max(1L, ceiling(maximum))
    ticks <- seq.int(0L, top_tick)
    ylimits <- c(0, max(maximum * 1.08 + 0.1, top_tick))
    if (nrow(annotated)) {
        positions <- volcano_label_positions(annotated, range(data$logFC), ylimits[[2L]])
        annotated$label_x <- positions$x
        annotated$label_y <- positions$y
        ylimits[[2L]] <- max(ylimits[[2L]], max(annotated$label_y) + 0.12)
    }
    y_label <- if (significance == "q_value") {
        "-log10(q-value)"
    } else {
        "-log10(p-value)"
    }
    colours <- stats::setNames(c("#777777", "#4C78A8", "#B24A59"),
        legend_labels)
    annotated_colour <- "#D98C70"
    if (style == "BPG") {
        significant <- data$significant
        plot <- BoutrosLab.plotting.general::create.scatterplot(
            formula = negative_log10 ~ logFC,
            data = data,
            groups = significant,
            col = unname(colours),
            pch = 19,
            cex = 0.8,
            alpha = 0.7,
            xlab.label = "Log2 fold change",
            ylab.label = y_label,
            xlab.cex = 1,
            ylab.cex = 1,
            xaxis.cex = 0.8,
            yaxis.cex = 0.8,
            yat = ticks,
            yaxis.lab = as.character(ticks),
            ylimits = ylimits,
            add.grid = FALSE,
            xaxis.tck = c(0.5, 0),
            yaxis.tck = c(0.5, 0),
            key = list(
                title = as.expression(bquote(bold(.(legend_title)))),
                cex.title = 1,
                points = list(col = rev(unname(colours)), pch = 19),
                text = list(lab = as.expression(list(
                    bquote(phantom(0) <= .(cutoff) ~ "; FC" >= 0),
                    bquote(phantom(0) <= .(cutoff) ~ "; FC" < 0),
                    bquote(phantom(0) > .(cutoff))
                )), cex = 0.8, font = 2),
                columns = 1,
                space = "right"
            ),
            left.padding = 3
        )
        if (nrow(annotated)) {
            original_panel <- plot$panel
            plot$panel <- function(x, y, ...) {
                original_panel(x, y, ...)
                lattice::panel.segments(annotated$logFC, annotated$negative_log10,
                    annotated$label_x, annotated$label_y, col = "#555555", lwd = 0.7)
                lattice::panel.points(annotated$logFC, annotated$negative_log10,
                    pch = 19, col = annotated_colour, cex = 0.8)
                lattice::panel.text(annotated$label_x, annotated$label_y,
                    labels = annotated$label, cex = 0.8, font = 2)
            }
        }
        return(plot)
    }
    plot <- ggplot2::ggplot(data,
        ggplot2::aes(x = .data$logFC, y = .data$negative_log10,
            colour = .data$significant)) +
        ggplot2::geom_point(size = 1.8, alpha = 0.7, show.legend = TRUE) +
        ggplot2::scale_colour_manual(values = colours, drop = FALSE,
            name = legend_title) +
        ggplot2::scale_y_continuous(breaks = ticks, labels = ticks,
            limits = ylimits) +
        ggplot2::labs(x = "Log2 fold change", y = y_label) +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::theme(
            panel.grid = ggplot2::element_blank(),
            panel.border = ggplot2::element_rect(
                colour = "black", fill = NA, linewidth = 0.5
            ),
            axis.title = ggplot2::element_text(face = "bold", size = 12),
            axis.text = ggplot2::element_text(colour = "black", face = "bold",
                size = 11),
            legend.title = ggplot2::element_text(face = "bold", size = 12),
            legend.text = ggplot2::element_text(face = "bold", size = 11),
            legend.position = "right"
        )
    if (nrow(annotated)) {
        plot <- plot + ggplot2::geom_segment(
            data = annotated,
            mapping = ggplot2::aes(xend = .data$label_x,
                yend = .data$label_y),
            colour = "#555555",
            linewidth = 0.35,
            show.legend = FALSE
        ) + ggplot2::geom_point(
            data = annotated,
            colour = annotated_colour,
            size = 1.8,
            alpha = 1,
            show.legend = FALSE
        ) + ggplot2::geom_text(
            data = annotated,
            mapping = ggplot2::aes(x = .data$label_x, y = .data$label_y,
                label = .data$label),
            size = 3.5,
            fontface = "bold",
            position = ggplot2::position_identity(),
            show.legend = FALSE
        )
    }
    plot
}

volcano_protein_labels <- function(data, proteins, label_by = NULL) {
    labels <- rep("", nrow(data))
    if (is.null(proteins) || !length(proteins)) return(labels)
    if (!is.null(label_by) && (!label_by %in% names(data) ||
        !is.character(data[[label_by]]))) {
        stop("label_by must name a character rowData column in the fitted result",
            call. = FALSE)
    }
    matches <- lapply(seq_len(nrow(data)), function(index) {
        terms <- as.character(data$feature_id[[index]])
        if ("source_group_id" %in% names(data)) {
            terms <- c(terms, as.character(data$source_group_id[[index]]))
        }
        for (column in c("leading_protein_accession", "member_accessions",
            "member_names")) {
            if (column %in% names(data)) terms <- c(terms, data[[column]][[index]])
        }
        for (column in grep("^(leading|member)_", names(data), value = TRUE)) {
            value <- data[[column]][[index]]
            if (is.character(value)) terms <- c(terms,
                unlist(strsplit(value[!is.na(value)], ";", fixed = TRUE),
                    use.names = FALSE))
        }
        proteins[proteins %in% terms]
    })
    found <- unique(unlist(matches, use.names = FALSE))
    missing <- setdiff(proteins, found)
    if (length(missing)) {
        stop("No plotted feature matches proteins: ", paste(missing, collapse = ", "),
            call. = FALSE)
    }
    if (is.null(label_by)) {
        return(vapply(matches, paste, collapse = ", ", FUN.VALUE = character(1)))
    }
    labels <- as.character(data[[label_by]])
    selected <- lengths(matches) > 0L
    if (any(selected & (is.na(labels) | !nzchar(labels)))) {
        stop("Selected features have missing values in label_by", call. = FALSE)
    }
    labels[!selected] <- ""
    labels
}

volcano_label_positions <- function(annotated, x_range, y_max) {
    x <- annotated$logFC
    label_x <- x
    y <- annotated$negative_log10 + 0.06 * y_max
    x_span <- max(diff(x_range), 1)
    y_step <- max(y_max * 0.065, 0.11)
    order <- order(annotated$negative_log10, annotated$logFC)
    placed <- integer()
    for (index in order) {
        label_x[[index]] <- x[[index]] + if (length(placed) %% 2L) {
            -0.07 * x_span
        } else {
            0.07 * x_span
        }
        nearby <- placed[abs(x[placed] - x[[index]]) < 0.28 * x_span]
        if (length(nearby)) {
            while (any(abs(y[[index]] - y[nearby]) < y_step)) {
                y[[index]] <- y[[index]] + y_step
            }
        }
        placed <- c(placed, index)
    }
    list(x = label_x, y = y)
}

#' Plot a limma p-value histogram
#'
#' Displays the raw `P.Value` values for one fitted contrast. Missing values
#' are omitted; p-values must lie between zero and one. Bars have black borders
#' in both plotting styles.
#'
#' @param x An object returned by [fit_limma()].
#' @param contrast Name of the contrast to plot. May be omitted only if the fit
#'   contains one contrast.
#' @param bins Number of equal-width bins between zero and one; default 40.
#' @param style `"BPG"` (default) or `"ggplot2"`.
#' @return A BPG trellis plot or a ggplot object, according to `style`.
#' @export
plot_limma_pvalue_histogram <- function(
    x,
    contrast = NULL,
    bins = 40L,
    style = c("BPG", "ggplot2")
) {
    if (!inherits(x, "proteovault_limma") || !is.data.frame(x$results)) {
        stop("x must be a proteovault_limma result", call. = FALSE)
    }
    style <- match.arg(style)
    if (!is.numeric(bins) || length(bins) != 1L || is.na(bins) ||
        !is.finite(bins) || bins < 1 || bins != as.integer(bins)) {
        stop("bins must be one positive integer", call. = FALSE)
    }
    available <- unique(as.character(x$results$contrast))
    if (is.null(contrast)) {
        if (length(available) != 1L) {
            stop("Select one contrast from: ", paste(available, collapse = ", "),
                call. = FALSE)
        }
        contrast <- available
    }
    if (!is.character(contrast) || length(contrast) != 1L ||
        is.na(contrast) || !contrast %in% available) {
        stop("contrast must name one fitted comparison: ",
            paste(available, collapse = ", "), call. = FALSE)
    }
    values <- x$results$P.Value[x$results$contrast == contrast]
    if (!is.numeric(values) || any(values < 0 | values > 1, na.rm = TRUE)) {
        stop("P.Value must contain numeric values between 0 and 1", call. = FALSE)
    }
    values <- values[is.finite(values)]
    if (!length(values)) {
        stop("No finite p-values for this contrast", call. = FALSE)
    }
    breaks <- seq(0, 1, length.out = as.integer(bins) + 1L)
    if (style == "BPG") {
        return(BoutrosLab.plotting.general::create.histogram(
            x = values,
            breaks = breaks,
            type = "count",
            col = "#4E9DA6",
            border.col = "black",
            xlimits = c(0, 1),
            xlab.label = "p-value",
            ylab.label = "Count",
            xlab.cex = 1,
            ylab.cex = 1,
            xaxis.cex = 0.8,
            yaxis.cex = 0.8,
            xaxis.tck = c(0.5, 0),
            yaxis.tck = c(0.5, 0),
            key = list(text = list(lab = ""))
        ))
    }
    data <- data.frame(p_value = values)
    ggplot2::ggplot(data, ggplot2::aes(x = .data$p_value)) +
        ggplot2::geom_histogram(breaks = breaks, fill = "#4E9DA6",
            colour = "black", linewidth = 0.5, closed = "right") +
        ggplot2::scale_x_continuous(limits = c(0, 1),
            breaks = seq(0, 1, by = 0.2), expand = ggplot2::expansion(mult = 0)) +
        ggplot2::labs(x = "p-value", y = "Count") +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::theme(
            panel.grid = ggplot2::element_blank(),
            panel.border = ggplot2::element_rect(
                colour = "black", fill = NA, linewidth = 0.5
            ),
            axis.title = ggplot2::element_text(face = "bold"),
            axis.text = ggplot2::element_text(colour = "black", face = "bold")
        )
}
