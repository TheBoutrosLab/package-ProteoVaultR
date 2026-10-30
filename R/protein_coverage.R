parse_openms_modifications <- function(modified_sequence, peptide_sequence) {
    characters <- strsplit(modified_sequence, split = "", fixed = TRUE)[[1L]]
    residues <- character()
    positions <- integer()
    names <- character()
    locations <- character()
    index <- 1L
    terminal_marker <- FALSE
    while (index <= length(characters)) {
        character <- characters[[index]]
        if (character %in% LETTERS) {
            residues <- c(residues, character)
            terminal_marker <- FALSE
            index <- index + 1L
            next
        }
        if (character %in% c(".", "-")) {
            terminal_marker <- TRUE
            index <- index + 1L
            next
        }
        if (!character %in% c("(", "[")) {
            stop(
                "Unsupported OpenMS modified sequence: ", modified_sequence,
                call. = FALSE
            )
        }
        closing <- if (character == "(") ")" else "]"
        depth <- 1L
        start <- index + 1L
        index <- start
        while (index <= length(characters) && depth > 0L) {
            if (characters[[index]] == character) depth <- depth + 1L
            if (characters[[index]] == closing) depth <- depth - 1L
            index <- index + 1L
        }
        if (depth != 0L || index <= start + 1L) {
            stop(
                "Malformed OpenMS modification: ", modified_sequence,
                call. = FALSE
            )
        }
        name <- paste(
            characters[seq.int(from = start, to = index - 2L)],
            collapse = ""
        )
        names <- c(names, name)
        if (!length(residues)) {
            positions <- c(positions, NA_integer_)
            locations <- c(locations, "N-term")
        } else if (terminal_marker) {
            positions <- c(positions, NA_integer_)
            locations <- c(locations, "C-term")
        } else {
            positions <- c(positions, length(residues))
            locations <- c(locations, "residue")
        }
        terminal_marker <- FALSE
    }
    if (paste(residues, collapse = "") != peptide_sequence) {
        stop(
            "Modified and unmodified peptide sequences disagree: ",
            modified_sequence,
            call. = FALSE
        )
    }
    data.frame(
        peptide_position = positions,
        modification = names,
        location = locations
    )
}


#' Plot validated peptide coverage across a protein sequence
#'
#' The sequence is wrapped into fixed-width rows. Each row shows the amino acids
#' above one coverage strip per sample. A peptide matching more than one
#' position in the same protein is shown as a possible location at every match.
#' Peptides shared with other proteins are included for each referenced protein;
#' this plot does not establish protein-specific identification.
#'
#' @inheritParams close_vault
#' @param accession Exact protein accession in the vault.
#' @param samples Sample names to show, in display order. `NULL` shows every
#'   sample in the vault.
#' @param row_width Number of residues per wrapped row.
#' @param show_ptms Whether to mark residue-level OpenMS modifications from
#'   validated PSMs.
#' @param ptm_names Optional modification names to show, such as `Oxidation`.
#'   `NULL` shows all observed residue-level modifications except those in
#'   `ptm_denylist`.
#' @param ptm_denylist Modification names to hide. By default,
#'   `Carbamidomethyl` is omitted because it is usually a fixed modification.
#'   Use `character()` to show it.
#' @return A ggplot object. The subtitle reports the union of residues covered
#'   by peptides with a single position in this protein.
#' @export
plot_protein_coverage <- function(
    vault,
    accession,
    samples = NULL,
    row_width = 50L,
    show_ptms = TRUE,
    ptm_names = NULL,
    ptm_denylist = "Carbamidomethyl"
) {
    if (!is.character(accession) || length(accession) != 1L ||
        is.na(accession) || !nzchar(accession)) {
        stop("accession must be one nonempty string", call. = FALSE)
    }
    if (!is.numeric(row_width) || length(row_width) != 1L ||
        is.na(row_width) || !is.finite(row_width) ||
        row_width < 1L || row_width != as.integer(row_width)) {
        stop("row_width must be one positive integer", call. = FALSE)
    }
    row_width <- as.integer(row_width)
    if (!is.logical(show_ptms) || length(show_ptms) != 1L ||
        is.na(show_ptms)) {
        stop("show_ptms must be TRUE or FALSE", call. = FALSE)
    }
    if (!is.null(ptm_names) &&
        (!is.character(ptm_names) || anyNA(ptm_names) ||
            any(!nzchar(ptm_names)) || anyDuplicated(ptm_names))) {
        stop("ptm_names must contain distinct nonempty names", call. = FALSE)
    }
    if (!is.character(ptm_denylist) || anyNA(ptm_denylist) ||
        any(!nzchar(ptm_denylist)) || anyDuplicated(ptm_denylist)) {
        stop("ptm_denylist must contain distinct nonempty names", call. = FALSE)
    }
    connection <- vault_connection(vault = vault)
    protein <- DBI::dbGetQuery(
        conn = connection,
        statement = "SELECT sequence FROM proteins WHERE accession = ?",
        params = list(accession)
    )
    if (nrow(protein) != 1L) {
        stop("Expected exactly one protein with accession ", accession, call. = FALSE)
    }
    if (is.na(protein$sequence[[1L]]) || !nzchar(protein$sequence[[1L]])) {
        stop(
            "Protein sequence is missing for ", accession,
            "; add the search FASTA with proteovault add-fasta",
            call. = FALSE
        )
    }
    protein_sequence <- toupper(protein$sequence[[1L]])
    protein_length <- nchar(protein_sequence)

    available_samples <- DBI::dbGetQuery(
        conn = connection,
        statement = "SELECT sample_name FROM samples ORDER BY sample_name"
    )$sample_name
    if (is.null(samples)) {
        samples <- available_samples
    } else if (!is.character(samples) || !length(samples) ||
        anyNA(samples) || any(!nzchar(samples)) || anyDuplicated(samples)) {
        stop("samples must contain distinct nonempty sample names", call. = FALSE)
    }
    if (!length(samples)) {
        stop("The vault has no samples to plot", call. = FALSE)
    }
    missing_samples <- setdiff(samples, available_samples)
    if (length(missing_samples)) {
        stop(
            "Unknown samples: ", paste(missing_samples, collapse = ", "),
            call. = FALSE
        )
    }

    evidence <- DBI::dbGetQuery(
        conn = connection,
        statement = paste(
            "SELECT DISTINCT s.sample_name, pep.sequence AS peptide_sequence,",
            "form.modified_sequence, form.notation",
            "FROM proteins prot",
            "JOIN psm_proteins pp ON pp.protein_id = prot.protein_id",
            "JOIN psms psm ON psm.psm_id = pp.psm_id",
            "JOIN processing_stages stage ON stage.stage_id = psm.stage_id",
            "JOIN spectra spectrum ON spectrum.spectrum_id = psm.spectrum_id",
            "JOIN runs run ON run.run_id = spectrum.run_id",
            "JOIN samples s ON s.sample_id = run.sample_id",
            "JOIN precursors precursor ON precursor.precursor_id = psm.precursor_id",
            "JOIN peptidoforms form ON form.peptidoform_id = precursor.peptidoform_id",
            "JOIN peptides pep ON pep.peptide_id = form.peptide_id",
            "WHERE prot.accession = ? AND stage.name = 'validation'"
        ),
        params = list(accession)
    )
    evidence <- evidence[evidence$sample_name %in% samples, , drop = FALSE]
    observed <- matrix(
        data = FALSE,
        nrow = protein_length,
        ncol = length(samples),
        dimnames = list(NULL, samples)
    )
    possible <- observed
    ptm_markers <- list()
    unplaced_ptms <- character()
    for (index in seq_len(nrow(evidence))) {
        peptide <- toupper(evidence$peptide_sequence[[index]])
        if (is.na(peptide) || !nzchar(peptide)) {
            stop("A validated PSM has an empty peptide sequence", call. = FALSE)
        }
        peptide_length <- nchar(peptide)
        candidate_starts <- seq_len(max(0L, protein_length - peptide_length + 1L))
        starts <- candidate_starts[
            substring(
                text = protein_sequence,
                first = candidate_starts,
                last = candidate_starts + peptide_length - 1L
            ) == peptide
        ]
        if (!length(starts)) {
            stop(
                "Validated peptide ", peptide,
                " does not match the stored sequence of ", accession,
                call. = FALSE
            )
        }
        positions <- unlist(
            lapply(
                X = starts,
                FUN = function(start) {
                    seq.int(from = start, length.out = peptide_length)
                }
            ),
            use.names = FALSE
        )
        sample_index <- match(evidence$sample_name[[index]], samples)
        if (length(starts) == 1L) {
            observed[positions, sample_index] <- TRUE
        } else {
            possible[positions, sample_index] <- TRUE
        }
        if (show_ptms && evidence$notation[[index]] == "OpenMS") {
            modifications <- parse_openms_modifications(
                modified_sequence = evidence$modified_sequence[[index]],
                peptide_sequence = evidence$peptide_sequence[[index]]
            )
            if (!is.null(ptm_names)) {
                modifications <- modifications[
                    modifications$modification %in% ptm_names,
                    ,
                    drop = FALSE
                ]
            }
            modifications <- modifications[
                !modifications$modification %in% ptm_denylist,
                ,
                drop = FALSE
            ]
            unplaced_ptms <- c(
                unplaced_ptms,
                modifications$modification[modifications$location != "residue"]
            )
            modifications <- modifications[
                modifications$location == "residue",
                ,
                drop = FALSE
            ]
            for (site in seq_len(nrow(modifications))) {
                ptm_markers[[length(ptm_markers) + 1L]] <- data.frame(
                    position = starts +
                        modifications$peptide_position[[site]] - 1L,
                    sample_name = evidence$sample_name[[index]],
                    modification = modifications$modification[[site]],
                    placement = if (length(starts) == 1L) {
                        "One position"
                    } else {
                        "Multiple positions"
                    }
                )
            }
        }
    }
    if (length(unplaced_ptms)) {
        warning(
            "Terminal modifications are not shown on residue tracks: ",
            paste(sort(unique(unplaced_ptms)), collapse = ", "),
            call. = FALSE
        )
    }
    residue_positions <- seq_len(protein_length)
    coverage <- data.frame(
        position = rep(residue_positions, times = length(samples)),
        sample_name = rep(samples, each = protein_length),
        status = ifelse(
            as.vector(observed),
            "Observed",
            ifelse(as.vector(possible), "Possible location", "Not observed")
        )
    )
    sequence_row <- data.frame(
        position = residue_positions,
        sample_name = "Sequence",
        status = NA_character_
    )
    coverage <- rbind(coverage, sequence_row)
    coverage$block_number <- (coverage$position - 1L) %/% row_width + 1L
    coverage$column <- (coverage$position - 1L) %% row_width + 1L
    block_start <- (coverage$block_number - 1L) * row_width + 1L
    block_end <- pmin(coverage$block_number * row_width, protein_length)
    coverage$block <- sprintf("%d\u2013%d", block_start, block_end)
    block_levels <- unique(coverage$block[coverage$sample_name == "Sequence"])
    coverage$block <- factor(coverage$block, levels = block_levels)
    coverage$residue <- ""
    sequence_indices <- coverage$sample_name == "Sequence"
    coverage$residue[sequence_indices] <- substring(
        text = protein_sequence,
        first = residue_positions,
        last = residue_positions
    )
    sample_coverage <- colSums(observed) / protein_length * 100
    sample_labels <- stats::setNames(
        sprintf("%s (%.1f%%)", samples, sample_coverage),
        samples
    )
    coverage$sample_label <- ifelse(
        sequence_indices,
        "Sequence",
        sample_labels[coverage$sample_name]
    )
    coverage$sample_label <- factor(
        coverage$sample_label,
        levels = c(rev(unname(sample_labels)), "Sequence")
    )
    coverage$status <- factor(
        coverage$status,
        levels = c("Not observed", "Possible location", "Observed")
    )
    union_coverage <- sum(rowSums(observed) > 0L)
    subtitle <- sprintf(
        "%d/%d residues covered by uniquely placed peptides (%.1f%%)",
        union_coverage,
        protein_length,
        100 * union_coverage / protein_length
    )

    plot <- ggplot2::ggplot(
        data = coverage,
        mapping = ggplot2::aes(x = column, y = sample_label, fill = status)
    ) +
        ggplot2::geom_tile(
            colour = "#DDDDDD",
            linewidth = 0.15,
            width = 1,
            height = 0.9
        ) +
        ggplot2::geom_text(
            data = coverage[sequence_indices, , drop = FALSE],
            mapping = ggplot2::aes(label = residue),
            family = "mono",
            fontface = "bold",
            size = 2.5
        ) +
        ggplot2::facet_grid(
            rows = ggplot2::vars(block),
            drop = FALSE,
            switch = "y"
        ) +
        ggplot2::scale_fill_manual(
            values = c(
                "Not observed" = "#F5F5F5",
                "Possible location" = "#E8B365",
                "Observed" = "#4E9DA6"
            ),
            drop = TRUE,
            na.value = "white",
            na.translate = FALSE,
            guide = ggplot2::guide_legend(order = 1, ncol = 1)
        ) +
        ggplot2::scale_x_continuous(
            breaks = c(
                1L,
                if (row_width >= 10L) {
                    seq.int(from = 10L, to = row_width, by = 10L)
                } else {
                    integer()
                }
            ),
            limits = c(0.5, row_width + 0.5),
            expand = ggplot2::expansion(mult = 0)
        ) +
        ggplot2::scale_y_discrete(
            drop = FALSE,
            expand = ggplot2::expansion(add = 0.5)
        ) +
        ggplot2::labs(
            title = accession,
            subtitle = subtitle,
            x = "Position within sequence row",
            y = NULL,
            fill = "Coverage"
        ) +
        ggplot2::theme_bw(base_size = 11) +
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
            strip.background = ggplot2::element_rect(
                fill = "white",
                colour = "black",
                linewidth = 0.5
            ),
            strip.text = ggplot2::element_text(face = "bold"),
            legend.position = "right",
            legend.direction = "vertical",
            legend.box = "vertical",
            legend.box.just = "left"
        )
    if (length(ptm_markers)) {
        ptms <- unique(do.call(rbind, ptm_markers))
        ptms$block_number <- (ptms$position - 1L) %/% row_width + 1L
        ptms$column <- (ptms$position - 1L) %% row_width + 1L
        ptms$block <- factor(
            block_levels[ptms$block_number],
            levels = block_levels
        )
        ptms$sample_label <- factor(
            sample_labels[ptms$sample_name],
            levels = levels(coverage$sample_label)
        )
        ptms$placement <- factor(
            ptms$placement,
            levels = c("One position", "Multiple positions")
        )
        ambiguous_ptms <- any(ptms$placement == "Multiple positions")
        names <- sort(unique(ptms$modification))
        colours <- stats::setNames(
            grDevices::hcl.colors(length(names), palette = "Dark 3"),
            names
        )
        plot <- plot +
            ggplot2::geom_point(
                data = ptms,
                mapping = ggplot2::aes(
                    x = column,
                    y = sample_label,
                    colour = modification,
                    shape = placement
                ),
                inherit.aes = FALSE,
                size = 2.2,
                stroke = 0.7
            ) +
            ggplot2::scale_colour_manual(
                values = colours,
                name = "PTM",
                guide = ggplot2::guide_legend(order = 2, ncol = 1)
            ) +
            ggplot2::scale_shape_manual(
                values = c(
                    "One position" = 16,
                    "Multiple positions" = 1
                ),
                name = NULL,
                breaks = "Multiple positions",
                labels = "Ambiguous peptide position",
                drop = TRUE,
                guide = if (ambiguous_ptms) {
                    ggplot2::guide_legend(order = 3, ncol = 1)
                } else {
                    "none"
                }
            )
    }
    plot
}
