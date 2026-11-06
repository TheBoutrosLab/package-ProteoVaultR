#' Plot precursor quantity line profiles for selected proteins
#'
#' Each line is a modified peptide and charge state. Protein values may be
#' source-group names or member accessions; names on `proteins` provide short
#' facet labels. Group-to-precursor and peptide-to-protein links identify
#' associated precursors. Only positive finite quantities are
#' plotted. Lines may span missing samples, as in the source notebook.
#'
#' @param vault An open ProteoVault.
#' @param proteins Distinct protein identifiers, optionally named with plot labels.
#' @param method Precursor quantification method, default `"Precursor.Normalised"`.
#' @param source_id Optional precursor quantification source identifier.
#' @param samples Optional original sample names, in plotting order.
#' @param colours Optional colors for precursors. Supply one color per plotted
#'   precursor in legend order, or a named vector keyed by precursor label.
#' @param style `"BPG"` (default) or `"ggplot2"`.
#' @return A BPG trellis plot or ggplot object.
#' @export
plot_protein_precursor_quantity_lineplot <- function(vault, proteins,
    method = "Precursor.Normalised", source_id = NULL, samples = NULL,
    colours = NULL, style = c("BPG", "ggplot2")) {
    style <- match.arg(style)
    if (!is.character(proteins) || !length(proteins) || anyNA(proteins) ||
        any(!nzchar(proteins)) || anyDuplicated(proteins))
        stop("proteins must contain distinct nonempty identifiers", call. = FALSE)
    labels <- names(proteins)
    if (is.null(labels)) labels <- proteins
    if (anyNA(labels) || any(!nzchar(labels)) || anyDuplicated(labels))
        stop("Protein labels must be distinct and nonempty", call. = FALSE)
    query <- tbl_precursor_quant(vault, source_id = source_id, method = method)
    available <- dplyr::collect(dplyr::distinct(query, .data$sample_name))$sample_name
    if (is.null(samples)) {
        metadata <- sample_metadata(vault)
        samples <- if (is.null(metadata)) sort(available) else
            metadata$sample_name[metadata$sample_name %in% available]
    }
    if (!is.character(samples) || !length(samples) || anyNA(samples) ||
        anyDuplicated(samples) || !all(samples %in% available))
        stop("samples must select distinct names in this quantification source",
            call. = FALSE)
    groups <- dplyr::collect(tbl_protein_groups(vault,
        group_type = "quantification"))
    members <- dplyr::collect(tbl_group_members(vault,
        group_type = "quantification"))
    links <- dplyr::collect(tbl_vault(vault, "precursor_groups"))
    peptide_links <- dplyr::collect(tbl_vault(vault, "peptide_proteins"))
    precursors <- dplyr::collect(tbl_vault(vault, "precursors"))
    peptidoforms <- dplyr::collect(tbl_vault(vault, "peptidoforms"))
    proteins_table <- dplyr::collect(tbl_vault(vault, "proteins"))
    precursor_ids <- lapply(proteins, function(protein) {
        group_ids <- groups$group_id[groups$source_group_id == protein]
        member_ids <- members$protein_id[members$group_id %in% group_ids]
        member_ids <- unique(c(member_ids,
            proteins_table$protein_id[proteins_table$accession == protein]))
        group_ids <- unique(c(group_ids,
            members$group_id[members$protein_id %in% member_ids]))
        from_groups <- links$precursor_id[links$group_id %in% group_ids]
        peptide_ids <- peptide_links$peptide_id[
            peptide_links$protein_id %in% member_ids]
        from_peptides <- precursors$precursor_id[
            precursors$peptidoform_id %in%
                peptidoforms$peptidoform_id[peptidoforms$peptide_id %in% peptide_ids]]
        unique(as.character(c(from_groups, from_peptides)))
    })
    if (any(lengths(precursor_ids) == 0L))
        stop("No linked precursors for proteins: ",
            paste(proteins[lengths(precursor_ids) == 0L], collapse = ", "),
            call. = FALSE)
    ids <- unique(unlist(precursor_ids, use.names = FALSE))
    quantities <- query %>%
        dplyr::filter(.data$precursor_id %in% .env$ids,
            .data$sample_name %in% .env$samples) %>%
        dplyr::select(dplyr::all_of(c("precursor_id", "sample_name",
            "modified_sequence", "charge", "quantity"))) %>%
        dplyr::collect()
    quantities$precursor_id <- as.character(quantities$precursor_id)
    quantities <- quantities[is.finite(quantities$quantity) &
        quantities$quantity > 0, , drop = FALSE]
    if (!nrow(quantities))
        stop("No positive precursor quantities for the selected proteins",
            call. = FALSE)
    plot_data <- do.call(rbind, lapply(seq_along(proteins), function(index) {
        rows <- quantities[quantities$precursor_id %in% precursor_ids[[index]],
            , drop = FALSE]
        rows$protein <- rep(labels[[index]], nrow(rows))
        rows
    }))
    plot_data <- plot_data[!is.na(plot_data$protein), , drop = FALSE]
    if (!nrow(plot_data)) stop("No plotted precursors", call. = FALSE)
    absent <- labels[!labels %in% plot_data$protein]
    if (length(absent))
        stop("No positive precursor quantities for proteins: ",
            paste(absent, collapse = ", "), call. = FALSE)
    plot_data$protein <- factor(plot_data$protein, levels = labels)
    plot_data$precursor <- paste(plot_data$modified_sequence,
        paste0("z", plot_data$charge))
    plot_data$log2_quantity <- log2(plot_data$quantity)
    plot_data$sample_index <- match(plot_data$sample_name, samples)
    short <- vault_sample_labels(vault)
    sample_labels <- samples
    if (!is.null(short)) {
        index <- match(samples, names(short))
        found <- !is.na(index)
        sample_labels[found] <- unname(short[index[found]])
    }
    if (anyNA(sample_labels) || any(!nzchar(sample_labels)) ||
        anyDuplicated(sample_labels)) sample_labels <- samples
    plot_data$sample <- factor(sample_labels[plot_data$sample_index],
        levels = sample_labels)
    precursor_names <- sort(unique(plot_data$precursor))
    if (is.null(colours)) {
        palette <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7",
            "#882255", "#44AA99", "#332288", "#AA4499",
            "#117733", "#661100", "#56B4E9", "#E69F00")
        colours <- if (length(precursor_names) <= length(palette)) {
            palette[seq_along(precursor_names)]
        } else {
            grDevices::hcl.colors(length(precursor_names), palette = "Dark 3")
        }
    } else {
        if (!is.character(colours) || anyNA(colours) ||
            any(!nzchar(colours)))
            stop("colours must contain valid color strings", call. = FALSE)
        if (!is.null(names(colours))) {
            if (anyDuplicated(names(colours)) ||
                !all(precursor_names %in% names(colours)))
                stop("Named colours must cover every plotted precursor",
                    call. = FALSE)
            colours <- colours[precursor_names]
        } else if (length(colours) != length(precursor_names)) {
            stop("colours must have one value per plotted precursor", call. = FALSE)
        }
        tryCatch(grDevices::col2rgb(colours), error = function(error) {
            stop("colours must contain valid color strings", call. = FALSE)
        })
    }
    colours <- stats::setNames(unname(colours), precursor_names)
    if (style == "BPG") {
        return(lattice::xyplot(log2_quantity ~ sample_index | protein,
            data = plot_data, groups = plot_data$precursor, type = "b",
            layout = c(length(labels), 1L),
            scales = list(x = list(at = seq_along(samples),
                labels = sample_labels, rot = 45, cex = 0.8, font = 2,
                alternating = FALSE, relation = "free", tck = c(0.5, 0)),
                y = list(relation = "free", cex = 0.8, font = 2,
                    tck = c(0.5, 0))),
            xlim = c(0.5, length(samples) + 0.5),
            xlab = NULL, ylab = "Log2 precursor quantity",
            auto.key = FALSE,
            key = list(
                points = list(col = c(NA, unname(colours)),
                    pch = c(NA, rep(19, length(colours)))),
                text = list(lab = c("precursor", names(colours)),
                    cex = c(1, rep(0.8, length(colours))), font = 2),
                columns = 1, space = "right"),
            par.settings = list(
                superpose.line = list(col = unname(colours), lwd = 1),
                superpose.symbol = list(col = unname(colours), pch = 19, cex = 0.8),
                strip.background = list(col = "grey85"),
                axis.line = list(col = "black", lwd = 1),
                par.ylab.text = list(font = 2),
                par.xlab.text = list(font = 2)),
            panel = function(x, y, groups, subscripts, ...) {
                panel_data <- plot_data[subscripts, , drop = FALSE]
                for (precursor in unique(panel_data$precursor)) {
                    rows <- panel_data[panel_data$precursor == precursor,
                        , drop = FALSE]
                    rows <- rows[order(rows$sample_index), , drop = FALSE]
                    lattice::panel.lines(rows$sample_index, rows$log2_quantity,
                        col = colours[[precursor]], lwd = 1)
                    lattice::panel.points(rows$sample_index, rows$log2_quantity,
                        col = colours[[precursor]], pch = 19, cex = 0.8)
                }
                limits <- lattice::current.panel.limits()
                lattice::panel.rect(limits$xlim[[1L]], limits$ylim[[1L]],
                    limits$xlim[[2L]], limits$ylim[[2L]],
                    col = NA, border = "black", lwd = 0.6)
            }))
    }
    ggplot2::ggplot(plot_data,
        ggplot2::aes(x = .data$sample, y = .data$log2_quantity,
            colour = .data$precursor, group = .data$precursor)) +
        ggplot2::geom_line(linewidth = 0.6) +
        ggplot2::geom_point(size = 2) +
        ggplot2::facet_wrap(ggplot2::vars(.data$protein), scales = "free_y") +
        ggplot2::scale_colour_manual(values = colours, name = "precursor") +
        ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(angle = 45)) +
        ggplot2::labs(x = NULL, y = "Log2 precursor quantity") +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::theme(panel.grid = ggplot2::element_blank(),
            panel.border = ggplot2::element_rect(colour = "black", fill = NA,
                linewidth = 0.6),
            axis.title = ggplot2::element_text(face = "bold"),
            axis.text = ggplot2::element_text(colour = "black", face = "bold"),
            strip.text = ggplot2::element_text(face = "bold"),
            legend.title = ggplot2::element_text(face = "bold"),
            legend.text = ggplot2::element_text(face = "bold"))
}
