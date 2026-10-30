# QC plots use the cohort metrics stored in the vault.

#' Plot precursor quantity distributions by sample
#'
#' Normalized quantities come from the core DIA-NN precursor matrix unless
#' `max_q_value` is set. A q-value cutoff selects the optional DIA-NN report
#' (`quant-details` layer), as do raw quantities.
#' Nonpositive and missing quantities are excluded before taking log2.
#'
#' @param vault An open ProteoVault connection returned by [open_vault()].
#' @param quantity_type `"normalized"`, `"raw"`, or `"both"`.
#' @param plot_type `"boxplot"` or `"violin"`.
#' @param source_id Optional source identifier. For `"both"`, provide a named
#'   vector with `normalized` and `raw` source IDs if selection is ambiguous.
#' @param samples Optional sample names to select and order.
#' @param max_q_value Optional inclusive run-specific precursor q-value cutoff
#'   between 0 and 1. Requires the `quant-details` layer. Cannot recover rows
#'   omitted by DIA-NN when it generated the report.
#' @return A ggplot object. Quantities are displayed on a log2 scale.
#' @export
plot_precursor_quantity_distribution <- function(
    vault,
    quantity_type = c("normalized", "raw", "both"),
    plot_type = c("boxplot", "violin"),
    source_id = NULL,
    samples = NULL,
    max_q_value = NULL
) {
    quantity_type <- match.arg(arg = quantity_type)
    plot_type <- match.arg(arg = plot_type)
    if (!is.null(max_q_value) &&
        (!is.numeric(max_q_value) || length(max_q_value) != 1L ||
            !is.finite(max_q_value) || max_q_value < 0 || max_q_value > 1)) {
        stop("max_q_value must be one number between 0 and 1", call. = FALSE)
    }

    methods <- c(
        normalized = if (is.null(max_q_value)) {
            "matrix"
        } else {
            "Precursor.Normalised"
        },
        raw = "Precursor.Quantity"
    )
    selected_types <- if (quantity_type == "both") {
        c("raw", "normalized")
    } else {
        quantity_type
    }
    if (quantity_type == "both" && !is.null(source_id) &&
        (!all(c("raw", "normalized") %in% names(source_id)) ||
            length(source_id) != 2L)) {
        stop(
            "For both quantity types, source_id must name raw and normalized sources",
            call. = FALSE
        )
    }
    quantities <- lapply(selected_types, function(type) {
        method <- methods[[type]]
        selected_source <- if (quantity_type == "both") {
            if (is.null(source_id)) NULL else source_id[[type]]
        } else {
            source_id
        }
        query <- tbl_precursor_quant(
            vault = vault,
            source_id = selected_source,
            method = method
        )
        if (!is.null(max_q_value)) {
            query <- query %>%
                dplyr::filter(
                    !is.na(.data$q_value),
                    .data$q_value <= .env$max_q_value
                )
        }
        query %>%
            dplyr::select(dplyr::all_of(c(
                "sample_name", "precursor_id", "quantity"
            ))) %>%
            dplyr::collect() %>%
            dplyr::mutate(quantity_type = type)
    })
    quantities <- dplyr::bind_rows(quantities)
    plot_quantity_distribution(
        quantities = quantities,
        samples = samples,
        plot_type = plot_type,
        quantity_types = selected_types,
        label = "precursor"
    )
}

#' Plot protein-group quantity distributions by sample
#'
#' Uses the normalized DIA-NN protein-group matrix. Nonpositive and missing
#' quantities are excluded before taking log2.
#'
#' @param vault An open ProteoVault connection returned by [open_vault()].
#' @param plot_type `"boxplot"` or `"violin"`.
#' @param source_id Optional DIA-NN protein-group matrix source identifier.
#' @param samples Optional sample names to select and order.
#' @return A ggplot object. Quantities are displayed on a log2 scale.
#' @export
plot_protein_quantity_distribution <- function(
    vault,
    plot_type = c("boxplot", "violin"),
    source_id = NULL,
    samples = NULL
) {
    plot_type <- match.arg(arg = plot_type)
    quantities <- tbl_protein_quant(
        vault = vault,
        source_id = source_id,
        method = "matrix"
    ) %>%
        dplyr::select(dplyr::all_of(c(
            "sample_name", "group_id", "quantity"
        ))) %>%
        dplyr::collect() %>%
        dplyr::mutate(quantity_type = "normalized")
    plot_quantity_distribution(
        quantities = quantities,
        samples = samples,
        plot_type = plot_type,
        quantity_types = "normalized",
        label = "protein-group"
    )
}

plot_quantity_distribution <- function(
    quantities,
    samples,
    plot_type,
    quantity_types,
    label
) {
    available <- sort(unique(quantities$sample_name))
    if (is.null(samples)) {
        samples <- available
    }
    if (!is.character(samples) || !length(samples) || anyNA(samples) ||
        anyDuplicated(samples) || !all(samples %in% available)) {
        stop("samples must be distinct names present in the quantities", call. = FALSE)
    }

    quantities <- quantities %>%
        dplyr::filter(
            .data$sample_name %in% .env$samples,
            !is.na(.data$quantity),
            is.finite(.data$quantity),
            .data$quantity > 0
        )
    if (!nrow(quantities)) {
        stop("No positive ", label, " quantities match the selection", call. = FALSE)
    }
    quantities$sample_name <- factor(
        quantities$sample_name,
        levels = samples
    )
    quantities$quantity_type <- factor(
        quantities$quantity_type,
        levels = quantity_types
    )
    quantities$log2_quantity <- log2(quantities$quantity)

    plot <- ggplot2::ggplot(
        data = quantities,
        mapping = ggplot2::aes(
            x = sample_name,
            y = log2_quantity
        )
    )
    if (plot_type == "boxplot") {
        plot <- plot + ggplot2::geom_boxplot(
            fill = "#4E9DA6",
            colour = "black",
            linewidth = 0.3,
            outlier.alpha = 0.2
        )
    } else {
        plot <- plot + ggplot2::geom_violin(
            fill = "#4E9DA6",
            colour = "black",
            linewidth = 0.3,
            scale = "width"
        )
    }
    if (length(quantity_types) > 1L) {
        plot <- plot + ggplot2::facet_wrap(
            facets = ggplot2::vars(quantity_type),
            ncol = 1
        )
    }
    plot +
        ggplot2::labs(
            x = "Sample",
            y = paste("Log2", label, "quantity")
        ) +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::theme(
            panel.grid = ggplot2::element_blank(),
            panel.border = ggplot2::element_rect(
                colour = "black",
                fill = NA,
                linewidth = 0.5
            ),
            axis.title = ggplot2::element_text(face = "bold"),
            axis.text = ggplot2::element_text(
                colour = "black",
                face = "bold"
            ),
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
        )
}

filter_qc_rows <- function(
    vault,
    metric_names,
    source_file_pattern,
    layer,
    source_process,
    search_engine,
    fraction_id
) {
    qc_results <- tbl_qc_metrics(vault = vault) %>%
        dplyr::filter(.data$metric_name %in% .env$metric_names) %>%
        dplyr::collect()
    required_columns <- c(
        "dataset_id", "sample_name", "metric_name", "group",
        "value", "unit", "source_file"
    )
    if (!is.data.frame(qc_results) ||
        !all(required_columns %in% names(qc_results))) {
        stop("QC table is missing required columns", call. = FALSE)
    }
    selectors <- list(
        layer = layer,
        source_process = source_process,
        search_engine = search_engine,
        fraction_id = fraction_id
    )
    for (name in names(selectors)) {
        selected <- selectors[[name]]
        if (is.null(selected)) {
            next
        }
        if (!is.character(selected) || length(selected) != 1L ||
            is.na(selected) || !name %in% names(qc_results)) {
            stop(name, " must be one string and a QC table column", call. = FALSE)
        }
    }
    if (!is.null(source_file_pattern) &&
        (!is.character(source_file_pattern) || length(source_file_pattern) != 1L ||
            is.na(source_file_pattern) || !nzchar(source_file_pattern))) {
        stop("source_file_pattern must be one nonempty string", call. = FALSE)
    }

    rows <- qc_results[
        !is.na(qc_results$metric_name) &
            qc_results$metric_name %in% metric_names,
        ,
        drop = FALSE
    ]
    for (name in names(selectors)) {
        selected <- selectors[[name]]
        if (is.null(selected)) {
            next
        }
        values <- as.character(rows[[name]])
        values[is.na(values)] <- ""
        rows <- rows[values == selected, , drop = FALSE]
    }
    if (!is.null(source_file_pattern)) {
        pattern <- if (startsWith(source_file_pattern, "*")) {
            utils::glob2rx(source_file_pattern)
        } else {
            source_file_pattern
        }
        matches <- !is.na(rows$source_file) & grepl(
            pattern = pattern,
            x = rows$source_file,
            perl = TRUE
        )
        rows <- rows[matches, , drop = FALSE]
    }
    if ("status" %in% names(rows)) {
        rows <- rows[!is.na(rows$status) & rows$status == "available", , drop = FALSE]
    }
    if (nrow(rows) == 0L) {
        stop("No available QC rows found for: ", paste(metric_names, collapse = ", "), call. = FALSE)
    }
    rows
}

#' Plot a cohort QC metric across samples.
#'
#' @param vault An open ProteoVault connection returned by [open_vault()].
#' @param metric_name One `metric_name` value from the vault QC table.
#' @param source_file_pattern Optional regex for source_file, or a glob beginning
#'   with `*` (for example, `*_MSFragger.idXML`).
#' @param layer Optional QC layer to select.
#' @param source_process Optional process name, such as `MSBooster`.
#' @param search_engine Optional search-engine name.
#' @param fraction_id Optional fraction ID; pass it as a string.
#' @param style Plotting style: `BPG` or `ggplot2`.
#' @return A BPG trellis plot or a ggplot object, according to `style`.
#' @export
plot_qc_metric <- function(
    vault,
    metric_name,
    source_file_pattern = NULL,
    layer = NULL,
    source_process = NULL,
    search_engine = NULL,
    fraction_id = NULL,
    style = c("BPG", "ggplot2")
) {
    style <- match.arg(arg = style)
    if (!is.character(metric_name) || length(metric_name) != 1L ||
        is.na(metric_name) || !nzchar(metric_name)) {
        stop("metric_name must be one nonempty string", call. = FALSE)
    }
    metric <- filter_qc_rows(
        vault = vault,
        metric_names = metric_name,
        source_file_pattern = source_file_pattern,
        layer = layer,
        source_process = source_process,
        search_engine = search_engine,
        fraction_id = fraction_id
    )

    metric$sample <- ifelse(
        is.na(metric$sample_name) | !nzchar(metric$sample_name),
        metric$dataset_id,
        metric$sample_name
    )
    if (anyNA(metric$sample) || any(!nzchar(metric$sample)) ||
        anyNA(metric$group) || any(!nzchar(metric$group))) {
        stop("QC rows need a dataset or sample name and a group", call. = FALSE)
    }
    if (anyDuplicated(metric[c("sample", "group")])) {
        stop(
            "Metric has multiple rows for the same sample and group; ",
            "filter by layer, source_process, search_engine, fraction_id, ",
            "or source_file_pattern",
            call. = FALSE
        )
    }
    metric$value <- suppressWarnings(as.numeric(metric$value))
    if (any(!is.finite(metric$value))) {
        stop("Metric contains nonnumeric or nonfinite values", call. = FALSE)
    }

    groups <- unique(metric$group)
    if (length(groups) > 1L && !"all" %in% groups) {
        stop("Grouped QC metrics must contain an 'all' group", call. = FALSE)
    }
    group_order <- c(
        intersect(c("all", "target", "decoy", "unknown"), groups),
        sort(setdiff(groups, c("all", "target", "decoy", "unknown")))
    )
    sample_order <- sort(unique(metric$sample))
    metric$sample <- factor(metric$sample, levels = rev(sample_order))
    metric$group <- factor(metric$group, levels = group_order)

    units <- unique(stats::na.omit(metric$unit))
    axis_label <- if (length(units) == 1L) {
        paste0(metric_name, " (", units, ")")
    } else {
        metric_name
    }
    value_range <- range(c(0, metric$value))
    padding <- max(diff(value_range) * 0.05, 0.01)
    xlimits <- value_range + c(-padding, padding)
    colours <- c(
        all = "#4E9DA6",
        target = "#6B9F59",
        decoy = "#D98C70",
        unknown = "#596B82"
    )
    missing_colours <- setdiff(group_order, names(colours))
    if (length(missing_colours)) {
        colours <- c(
            colours,
            stats::setNames(
                grDevices::hcl.colors(length(missing_colours), palette = "Dark 3"),
                missing_colours
            )
        )
    }
    colours <- colours[group_order]

    if (style == "ggplot2") {
        return(
            ggplot2::ggplot(
                data = metric,
                mapping = ggplot2::aes(x = sample, y = value, fill = group)
            ) +
                ggplot2::geom_col(
                    position = ggplot2::position_dodge(width = 0.8),
                    width = 0.7,
                    colour = "black",
                    linewidth = 0.3
                ) +
                ggplot2::scale_fill_manual(values = colours, drop = FALSE) +
                ggplot2::scale_x_discrete(
                    expand = ggplot2::expansion(add = c(0.5, 0.5))
                ) +
                ggplot2::coord_flip(ylim = xlimits) +
                ggplot2::labs(x = NULL, y = axis_label, fill = NULL) +
                ggplot2::theme_bw(base_size = 12) +
                ggplot2::theme(
                    panel.grid = ggplot2::element_blank(),
                    panel.border = ggplot2::element_rect(
                        colour = "black", fill = NA, linewidth = 0.5
                    ),
                    axis.title = ggplot2::element_text(face = "bold"),
                    axis.text = ggplot2::element_text(
                        colour = "black", face = "bold"
                    ),
                    legend.position = if (length(group_order) > 1L) "top" else "none"
                )
        )
    }

    if (length(group_order) == 1L) {
        plot <- BoutrosLab.plotting.general::create.barplot(
            formula = sample ~ value,
            data = metric,
            plot.horizontal = TRUE,
            col = unname(colours),
            border.col = "black",
            border.lwd = 0.3,
            xlimits = xlimits,
            xlab.label = axis_label,
            ylab.label = NULL,
            xlab.cex = 1,
            xaxis.cex = 0.8,
            yaxis.cex = 0.8,
            left.padding = 3,
            right.padding = 0.5,
            xaxis.tck = c(0.5, 0),
            yaxis.tck = c(0.5, 0)
        )
    } else {
        plot <- BoutrosLab.plotting.general::create.barplot(
            formula = sample ~ value,
            data = metric,
            groups = group,
            plot.horizontal = TRUE,
            col = unname(colours),
            border.col = "black",
            border.lwd = 0.3,
            key = list(
                rectangles = list(col = rev(unname(colours))),
                text = list(lab = rev(group_order)),
                columns = 1,
                space = "right"
            ),
            xlimits = xlimits,
            xlab.label = axis_label,
            ylab.label = NULL,
            xlab.cex = 1,
            xaxis.cex = 0.8,
            yaxis.cex = 0.8,
            left.padding = 3,
            right.padding = 0.5,
            top.padding = 1,
            xaxis.tck = c(0.5, 0),
            yaxis.tck = c(0.5, 0)
        )
    }

    plot
}

#' Plot QC percentile summaries as boxes with 5th and 95th percentile whiskers.
#'
#' @param vault An open ProteoVault connection returned by [open_vault()].
#' @param metric_name Metric stem, such as `dia_cycle_time`.
#' @param source_file_pattern Optional regex for source_file, or a glob beginning
#'   with `*`.
#' @param layer Optional QC layer to select.
#' @param source_process Optional process name, such as `MSBooster`.
#' @param search_engine Optional search-engine name.
#' @param fraction_id Optional fraction ID; pass it as a string.
#' @param style Plotting style: `ggplot2` or `BPG`. The BPG style uses
#'   `lattice::bwplot()` with a panel drawn from the supplied quantiles.
#' @return A ggplot object or lattice trellis plot, according to `style`.
#' @export
plot_qc_boxplot <- function(
    vault,
    metric_name,
    source_file_pattern = NULL,
    layer = NULL,
    source_process = NULL,
    search_engine = NULL,
    fraction_id = NULL,
    style = c("ggplot2", "BPG")
) {
    style <- match.arg(arg = style)
    if (!is.character(metric_name) || length(metric_name) != 1L ||
        is.na(metric_name) || !nzchar(metric_name)) {
        stop("metric_name must be one nonempty string", call. = FALSE)
    }
    percentiles <- c("p05", "p25", "p50", "p75", "p95")
    metric_name <- sub("_p(05|25|50|75|95)$", "", metric_name)
    metric <- filter_qc_rows(
        vault = vault,
        metric_names = paste0(metric_name, "_", percentiles),
        source_file_pattern = source_file_pattern,
        layer = layer,
        source_process = source_process,
        search_engine = search_engine,
        fraction_id = fraction_id
    )

    metric$sample <- ifelse(
        is.na(metric$sample_name) | !nzchar(metric$sample_name),
        metric$dataset_id,
        metric$sample_name
    )
    if (anyNA(metric$sample) || any(!nzchar(metric$sample)) ||
        anyNA(metric$group) || any(!nzchar(metric$group))) {
        stop("QC rows need a dataset or sample name and a group", call. = FALSE)
    }
    metric$percentile <- sub("^.*_", "", metric$metric_name)
    if (anyDuplicated(metric[c("sample", "group", "percentile")])) {
        stop(
            "Metric has repeated percentiles for the same sample and group; ",
            "filter by layer, source_process, search_engine, fraction_id, ",
            "or source_file_pattern",
            call. = FALSE
        )
    }
    metric$value <- suppressWarnings(as.numeric(metric$value))
    if (any(!is.finite(metric$value))) {
        stop("Metric contains nonnumeric or nonfinite percentiles", call. = FALSE)
    }

    summary_groups <- split(
        metric,
        interaction(metric$sample, metric$group, drop = TRUE)
    )
    summaries <- lapply(summary_groups, function(rows) {
        values <- stats::setNames(rows$value, rows$percentile)
        if (!setequal(names(values), percentiles)) {
            stop(
                "Each sample and group needs p05, p25, p50, p75, and p95",
                call. = FALSE
            )
        }
        values <- values[percentiles]
        if (is.unsorted(values, strictly = FALSE)) {
            stop("Percentiles must be nondecreasing", call. = FALSE)
        }
        data.frame(
            sample = rows$sample[[1]],
            group = rows$group[[1]],
            p05 = values[["p05"]],
            p25 = values[["p25"]],
            p50 = values[["p50"]],
            p75 = values[["p75"]],
            p95 = values[["p95"]]
        )
    })
    summary <- do.call(rbind, summaries)
    rownames(summary) <- NULL

    group_order <- c(
        intersect(c("all", "target", "decoy", "unknown"), unique(summary$group)),
        sort(setdiff(
            unique(summary$group),
            c("all", "target", "decoy", "unknown")
        ))
    )
    sample_order <- rev(sort(unique(summary$sample)))
    summary$sample <- factor(summary$sample, levels = sample_order)
    summary$group <- factor(summary$group, levels = group_order)
    colours <- c(
        all = "#4E9DA6",
        target = "#6B9F59",
        decoy = "#D98C70",
        unknown = "#596B82"
    )
    missing_colours <- setdiff(group_order, names(colours))
    if (length(missing_colours)) {
        colours <- c(
            colours,
            stats::setNames(
                grDevices::hcl.colors(length(missing_colours), palette = "Dark 3"),
                missing_colours
            )
        )
    }
    colours <- colours[group_order]

    units <- unique(metric$unit[!is.na(metric$unit) & nzchar(metric$unit)])
    axis_label <- if (length(units) == 1L) {
        paste0(metric_name, " (", units, ")")
    } else {
        metric_name
    }
    value_range <- range(c(summary$p05, summary$p95))
    padding <- if (diff(value_range) > 0) {
        diff(value_range) * 0.05
    } else {
        max(abs(value_range) * 0.05, 0.01)
    }

    if (style == "BPG") {
        group_count <- length(group_order)
        return(lattice::bwplot(
            sample ~ p50,
            data = summary,
            horizontal = TRUE,
            xlim = value_range + c(-padding, padding),
            xlab = axis_label,
            ylab = NULL,
            box.ratio = 0.65,
            key = if (group_count > 1L) {
                list(
                    rectangles = list(col = unname(colours)),
                    text = list(lab = group_order),
                    columns = group_count,
                    space = "top"
                )
            } else {
                NULL
            },
            par.settings = list(
                axis.text = list(font = 2),
                par.xlab.text = list(font = 2)
            ),
            panel = function(x, y, subscripts, ...) {
                rows <- summary[subscripts, , drop = FALSE]
                positions <- as.numeric(y)
                groups <- match(as.character(rows$group), group_order)
                width <- 0.65 / group_count
                positions <- positions +
                    (groups - (group_count + 1) / 2) * width
                half_width <- width * 0.42
                fill <- unname(colours[as.character(rows$group)])
                lattice::panel.segments(
                    x0 = rows$p05, y0 = positions,
                    x1 = rows$p95, y1 = positions,
                    col = "black", lwd = 0.5
                )
                lattice::panel.segments(
                    x0 = rows$p05, y0 = positions - half_width * 0.5,
                    x1 = rows$p05, y1 = positions + half_width * 0.5,
                    col = "black", lwd = 0.5
                )
                lattice::panel.segments(
                    x0 = rows$p95, y0 = positions - half_width * 0.5,
                    x1 = rows$p95, y1 = positions + half_width * 0.5,
                    col = "black", lwd = 0.5
                )
                lattice::panel.rect(
                    xleft = rows$p25,
                    ybottom = positions - half_width,
                    xright = rows$p75,
                    ytop = positions + half_width,
                    col = fill,
                    border = "black",
                    lwd = 0.5
                )
                lattice::panel.segments(
                    x0 = rows$p50, y0 = positions - half_width,
                    x1 = rows$p50, y1 = positions + half_width,
                    col = "black", lwd = 1.2
                )
            }
        ))
    }

    ggplot2::ggplot(
        data = summary,
        mapping = ggplot2::aes(
            x = sample,
            ymin = p05,
            lower = p25,
            middle = p50,
            upper = p75,
            ymax = p95,
            fill = group
        )
    ) +
        ggplot2::geom_boxplot(
            stat = "identity",
            position = ggplot2::position_dodge2(
                width = 0.75,
                preserve = "single"
            ),
            width = 0.65,
            staplewidth = 0.5,
            colour = "black",
            linewidth = 0.3
        ) +
        ggplot2::scale_fill_manual(
            values = colours,
            breaks = group_order,
            drop = FALSE
        ) +
        ggplot2::scale_x_discrete(
            expand = ggplot2::expansion(add = c(0.7, 0.7))
        ) +
        ggplot2::scale_y_continuous(
            expand = ggplot2::expansion(mult = 0)
        ) +
        ggplot2::coord_flip(
            ylim = value_range + c(-padding, padding)
        ) +
        ggplot2::labs(
            x = NULL,
            y = axis_label,
            fill = NULL
        ) +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::theme(
            panel.grid = ggplot2::element_blank(),
            panel.border = ggplot2::element_rect(
                colour = "black",
                fill = NA,
                linewidth = 1
            ),
            axis.title = ggplot2::element_text(face = "bold"),
            axis.text = ggplot2::element_text(colour = "black", face = "bold"),
            axis.ticks = ggplot2::element_line(colour = "black"),
            legend.position = if (length(group_order) > 1L) "top" else "none",
            legend.direction = "horizontal",
            legend.key = ggplot2::element_blank()
        )
}

#' Plot PSM counts by precursor charge across samples.
#'
#' @param vault An open ProteoVault connection returned by [open_vault()].
#' @param charges Charge states to display, in order.
#' @param source_file_pattern Optional regex for source_file, or a glob beginning
#'   with `*`.
#' @param layer Optional QC layer to select.
#' @param source_process QC source process; defaults to MSBooster.
#' @param search_engine Optional search-engine name.
#' @param fraction_id Optional fraction ID; pass it as a string.
#' @param style Plotting style: `BPG` or `ggplot2`.
#' @return A BPG trellis plot or a ggplot object, according to `style`.
#' @export
plot_qc_charge_psm_counts <- function(
    vault,
    charges = 1:4,
    source_file_pattern = NULL,
    layer = NULL,
    source_process = "MSBooster",
    search_engine = NULL,
    fraction_id = NULL,
    style = c("BPG", "ggplot2")
) {
    style <- match.arg(arg = style)
    if (!is.numeric(charges) || !length(charges) || anyNA(charges) ||
        any(charges < 1 | charges != as.integer(charges)) ||
        anyDuplicated(charges)) {
        stop("charges must contain distinct positive integers", call. = FALSE)
    }
    charge_names <- as.character(charges)
    metric <- filter_qc_rows(
        vault = vault,
        metric_names = paste0("feature_charge_", charge_names, "_psm_count"),
        source_file_pattern = source_file_pattern,
        layer = layer,
        source_process = source_process,
        search_engine = search_engine,
        fraction_id = fraction_id
    )
    metric <- metric[metric$group %in% c("all", "target", "decoy"), , drop = FALSE]
    if (!nrow(metric)) {
        stop("No all, target, or decoy charge counts found", call. = FALSE)
    }
    metric$sample <- ifelse(
        is.na(metric$sample_name) | !nzchar(metric$sample_name),
        metric$dataset_id,
        metric$sample_name
    )
    if (anyNA(metric$sample) || any(!nzchar(metric$sample))) {
        stop("Charge counts need a dataset or sample name", call. = FALSE)
    }
    metric$charge <- sub(
        pattern = "^feature_charge_([0-9]+)_psm_count$",
        replacement = "\\1",
        x = metric$metric_name
    )
    if (anyDuplicated(metric[c("sample", "group", "charge")])) {
        stop(
            "Multiple charge counts found for a sample and group; filter by ",
            "layer, source_process, search_engine, fraction_id, or source_file_pattern",
            call. = FALSE
        )
    }
    metric$value <- suppressWarnings(as.numeric(metric$value))
    if (any(!is.finite(metric$value)) || any(metric$value < 0)) {
        stop("Charge counts must be finite, nonnegative numbers", call. = FALSE)
    }
    metric$sample <- factor(metric$sample, levels = rev(sort(unique(metric$sample))))
    metric$group <- factor(metric$group, levels = c("all", "decoy", "target"))
    metric$charge <- factor(metric$charge, levels = charge_names)
    metric <- as.data.frame(metric)
    colours <- grDevices::hcl.colors(length(charges), palette = "Dark 3")

    if (style == "ggplot2") {
        return(
            ggplot2::ggplot(
                data = metric,
                mapping = ggplot2::aes(x = sample, y = value, fill = charge)
            ) +
                ggplot2::geom_col(
                    position = ggplot2::position_dodge(width = 0.8),
                    width = 0.7,
                    colour = "black",
                    linewidth = 0.3
                ) +
                ggplot2::facet_grid(cols = ggplot2::vars(group), drop = FALSE) +
                ggplot2::scale_fill_manual(
                    values = stats::setNames(colours, charge_names),
                    labels = paste0("Charge ", charge_names),
                    drop = FALSE
                ) +
                ggplot2::scale_x_discrete(
                    expand = ggplot2::expansion(add = c(0.5, 0.5))
                ) +
                ggplot2::coord_flip(ylim = c(0, max(metric$value) * 1.05 + 1)) +
                ggplot2::labs(x = NULL, y = "PSM count", fill = NULL) +
                ggplot2::theme_bw(base_size = 12) +
                ggplot2::theme(
                    panel.grid = ggplot2::element_blank(),
                    panel.border = ggplot2::element_rect(
                        colour = "black", fill = NA, linewidth = 0.5
                    ),
                    axis.title = ggplot2::element_text(face = "bold"),
                    axis.text = ggplot2::element_text(
                        colour = "black", face = "bold"
                    ),
                    strip.text = ggplot2::element_text(face = "bold"),
                    legend.position = "top"
                )
        )
    }

    BoutrosLab.plotting.general::create.barplot(
        formula = sample ~ value | group,
        data = metric,
        groups = charge,
        layout = c(3, 1),
        plot.horizontal = TRUE,
        col = colours,
        border.col = "black",
        border.lwd = 0.3,
        key = list(
            rectangles = list(col = colours),
            text = list(lab = paste0("Charge ", charge_names)),
            columns = 1,
            space = "right"
        ),
        xlimits = c(0, max(metric$value) * 1.05 + 1),
        xlab.label = "PSM count",
        ylab.label = NULL,
        xlab.cex = 1,
        xaxis.cex = 0.8,
        yaxis.cex = 0.8,
        left.padding = 3,
        right.padding = 0.5,
        top.padding = 1,
        xaxis.tck = c(0.5, 0),
        yaxis.tck = c(0.5, 0)
    )
}
