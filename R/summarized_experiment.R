#' Export quantification as a SummarizedExperiment
#'
#' Protein-group and precursor measurements are stored in the vault. DDA
#' ProteinQuantifier peptide rows are exported directly; DIA peptide quantities
#' are derived by aggregating precursors with the same stripped sequence.
#' For DIA-NN groups, the first accession in
#' `source_group_id` is treated as the leading protein. It is matched to a
#' stored group member to provide `leading_protein_id` and
#' `leading_protein_accession`. A feature linked to several groups retains
#' their leaders as list columns; the singular fields are missing then.
#' Each UniProt mapping type adds `leading_<type>` and `member_<type>` character
#' columns. Distinct values are joined with semicolons; missing values are `NA`.
#'
#' @inheritParams matrix_protein_abundance
#' @param level `"protein_group"`, `"precursor"`, or `"peptide"`.
#' @param mapping_types Optional UniProt mapping types to add as paired leading
#'   and member columns in rowData. `NULL` retains all available types. Requires an imported
#'   UniProt mapping source when specified.
#' @param mapping_source_id Optional UniProt mapping source identifier.
#' @param duplicates Resolution of repeated feature/sample measurements:
#'   `"error"` (default outside TMT), `"best"`, `"sum"`, or `"mean"`.
#'   TMT precursor exports default to `"sum"`; only `"sum"` and `"error"`
#'   are supported for TMT to avoid mixing reporter channels independently.
#' @param peptide_aggregation Required when deriving DIA peptides from
#'   precursors: `"sum"`, `"max"`, or `"mean"`. Native DDA peptide rows do not
#'   use this option.
#' @return A [SummarizedExperiment::SummarizedExperiment] with a `quantity`
#'   assay, feature annotations in rowData, and sample metadata in colData.
#' @export
as_summarized_experiment <- function(
    vault,
    level = c("protein_group", "precursor", "peptide"),
    method = NULL,
    source_id = NULL,
    samples = NULL,
    duplicates = c("error", "best", "sum", "mean"),
    mapping_types = NULL,
    mapping_source_id = NULL,
    peptide_aggregation = NULL
) {
    level <- match.arg(level)
    duplicates_missing <- missing(duplicates)
    duplicates <- match.arg(duplicates)
    if (!is.null(mapping_types) &&
        (!is.character(mapping_types) || anyNA(mapping_types) ||
            any(!nzchar(mapping_types)))) {
        stop("mapping_types must contain nonempty strings", call. = FALSE)
    }
    native_peptide <- level == "peptide" &&
        (identical(method, "ProteinQuantifier") ||
            (is.null(method) && any(vault$sources$format == "openms_peptide_quant")))
    if (level == "peptide" && !native_peptide) {
        if (is.null(peptide_aggregation)) {
            stop("peptide_aggregation is required for derived peptide quantities", call. = FALSE)
        }
        peptide_aggregation <- match.arg(peptide_aggregation, c("sum", "max", "mean"))
    } else if (!is.null(peptide_aggregation)) {
        stop("peptide_aggregation applies only to derived peptide quantities", call. = FALSE)
    }
    if (!"uniprot-idmapping" %in% vault_info(vault)$layers &&
        (!is.null(mapping_types) || !is.null(mapping_source_id))) {
        stop("Protein mappings require an imported UniProt idmapping source", call. = FALSE)
    }
    quant_level <- if (level == "protein_group") "protein_group" else if (
        native_peptide) "peptide" else "precursor"
    method <- resolve_quant_method(vault, quant_level, method)
    if (duplicates_missing && quant_level == "precursor" && startsWith(method, "TMT.")) {
        duplicates <- "sum"
    }
    source_id <- resolve_quant_source(vault, quant_level, source_id, method)
    if (level == "protein_group") {
        assay <- matrix_protein_abundance(vault, source_id, method, samples,
            row_id = "group_id", duplicates = duplicates)
        features <- protein_group_export_data(vault, source_id, rownames(assay),
            mapping_types, mapping_source_id)
    } else if (native_peptide) {
        assay <- matrix_peptide_abundance(vault, source_id, method, samples,
            duplicates)
        features <- native_peptide_export_data(vault, rownames(assay),
            mapping_types, mapping_source_id)
    } else {
        precursor_assay <- matrix_precursor_abundance(vault, source_id, method,
            samples, duplicates)
        precursors <- precursor_export_data(vault, rownames(precursor_assay))
        if (level == "precursor") {
            assay <- precursor_assay
            features <- precursors
        } else {
            groups <- split(seq_len(nrow(precursors)), precursors$peptide_id)
            assay <- matrix(NA_real_, nrow = length(groups),
                ncol = ncol(precursor_assay),
                dimnames = list(names(groups), colnames(precursor_assay)))
            choice <- if (peptide_aggregation == "max") "best" else peptide_aggregation
            for (index in seq_along(groups)) {
                assay[index, ] <- apply(precursor_assay[groups[[index]], , drop = FALSE],
                    2L, aggregate_duplicate_quantities, duplicates = choice)
            }
            features <- data.frame(
                peptide_id = rownames(assay),
                sequence = vapply(groups, function(indices) {
                    precursors$sequence[indices[[1L]]]
                }, character(1)),
                stringsAsFactors = FALSE
            )
            features$precursor_ids <- I(lapply(groups, function(indices) {
                precursors$precursor_id[indices]
            }))
        }
        features <- add_precursor_group_annotations(vault, features, level,
            source_id, mapping_types, mapping_source_id)
    }
    col_data <- sample_export_data(vault, colnames(assay))
    features <- as.data.frame(features)
    col_data <- as.data.frame(col_data)
    rownames(features) <- rownames(assay)
    rownames(col_data) <- colnames(assay)
    SummarizedExperiment::SummarizedExperiment(
        assays = list(quantity = assay),
        rowData = S4Vectors::DataFrame(features),
        colData = S4Vectors::DataFrame(col_data),
        metadata = list(level = level, method = method,
            source_id = as.character(source_id),
            peptide_aggregation = peptide_aggregation)
    )
}

sample_export_data <- function(vault, samples) {
    stored <- dplyr::collect(tbl_vault(vault, "samples"))
    selected <- stored[match(samples, stored$sample_name), , drop = FALSE]
    names(selected)[names(selected) == "sample_id"] <- "vault_sample_id"
    runs <- dplyr::collect(tbl_vault(vault, "runs"))
    selected$run_names <- I(lapply(selected$vault_sample_id, function(id) {
        runs$run_name[!is.na(runs$sample_id) &
            as.character(runs$sample_id) == as.character(id)]
    }))
    selected$file_names <- I(lapply(selected$vault_sample_id, function(id) {
        runs$file_name[!is.na(runs$sample_id) &
            as.character(runs$sample_id) == as.character(id)]
    }))
    tmt <- DBI::dbGetQuery(vault_connection(vault), paste(
        "SELECT s.sample_name, p.plex_name, t.reporter_label,",
        "t.channel_description, coalesce(t.channel_id = p.reference_channel_id, false) AS is_reference,",
        "r.run_name, r.file_name FROM samples s",
        "JOIN quant_channels c USING(sample_id)",
        "JOIN tmt_channels t USING(channel_id)",
        "JOIN tmt_plexes p USING(plex_id)",
        "LEFT JOIN tmt_runs tr USING(plex_id)",
        "LEFT JOIN runs r ON tr.run_id = r.run_id",
        "ORDER BY s.sample_name, r.run_name"
    ))
    if (nrow(tmt)) {
        index <- match(samples, tmt$sample_name)
        selected$tmt_plex <- tmt$plex_name[index]
        selected$tmt_reporter <- tmt$reporter_label[index]
        selected$tmt_channel_description <- tmt$channel_description[index]
        selected$tmt_reference <- tmt$is_reference[index]
        for (position in which(!is.na(index))) {
            rows <- tmt[tmt$sample_name == samples[[position]], , drop = FALSE]
            selected$run_names[[position]] <- unique(rows$run_name[!is.na(rows$run_name)])
            selected$file_names[[position]] <- unique(rows$file_name[!is.na(rows$file_name)])
        }
    }
    metadata <- sample_metadata(vault)
    if (!is.null(metadata)) {
        metadata <- metadata[match(samples, metadata$sample_name), , drop = FALSE]
        extra <- setdiff(names(metadata), "sample_name")
        collisions <- intersect(extra, names(selected))
        if (length(collisions)) {
            stop("Sample metadata conflicts with vault columns: ",
                paste(collisions, collapse = ", "), call. = FALSE)
        }
        for (name in extra) selected[[name]] <- metadata[[name]]
    }
    selected
}

native_peptide_export_data <- function(vault, ids, mapping_types, mapping_source_id) {
    features <- tbl_vault(vault, "peptide_quant_features") %>%
        dplyr::inner_join(tbl_vault(vault, "peptides"), by = "peptide_id") %>%
        dplyr::collect()
    features$feature_id <- as.character(features$feature_id)
    features$peptide_id <- as.character(features$peptide_id)
    features <- features[match(ids, features$feature_id), , drop = FALSE]
    names(features)[names(features) == "sequence"] <- "peptide"
    protein_source <- resolve_quant_source(vault, "protein_group",
        method = "ProteinQuantifier")
    groups <- dplyr::collect(tbl_protein_groups(vault,
        source_id = protein_source, method = "ProteinQuantifier"))
    groups$group_id <- as.character(groups$group_id)
    group_index <- match(features$source_protein, groups$source_group_id)
    features$protein_group_id <- groups$group_id[group_index]
    for (column in c("member_protein_ids", "member_accessions", "member_names")) {
        features[[column]] <- I(rep(list(character()), nrow(features)))
    }
    features$leading_protein_id <- rep(NA_character_, nrow(features))
    features$leading_protein_accession <- rep(NA_character_, nrow(features))
    matched_ids <- unique(stats::na.omit(features$protein_group_id))
    if (length(matched_ids)) {
        annotations <- protein_group_export_data(vault, protein_source,
            matched_ids, mapping_types, mapping_source_id)
        index <- match(features$protein_group_id, annotations$group_id)
        for (column in c("member_protein_ids", "member_accessions", "member_names")) {
            features[[column]] <- I(lapply(index, function(i) {
                if (is.na(i)) character() else annotations[[column]][[i]]
            }))
        }
        mapped <- grep("^(leading_|member_)", names(annotations), value = TRUE)
        mapped <- setdiff(mapped, c("member_protein_ids", "member_accessions",
            "member_names"))
        for (column in mapped) features[[column]] <- annotations[[column]][index]
    } else {
        features <- add_mapping_columns(vault, features, mapping_types,
            mapping_source_id)
    }
    features
}

precursor_export_data <- function(vault, ids) {
    features <- tbl_vault(vault, "precursors") %>%
        dplyr::inner_join(tbl_vault(vault, "peptidoforms"), by = "peptidoform_id") %>%
        dplyr::inner_join(tbl_vault(vault, "peptides"), by = "peptide_id") %>%
        dplyr::select(dplyr::all_of(c("precursor_id", "peptide_id", "sequence",
            "modified_sequence", "notation", "charge"))) %>%
        dplyr::collect()
    features$precursor_id <- as.character(features$precursor_id)
    features$peptide_id <- as.character(features$peptide_id)
    features[match(ids, features$precursor_id), , drop = FALSE]
}

protein_group_export_data <- function(vault, source_id, ids,
    mapping_types, mapping_source_id) {
    groups <- dplyr::collect(tbl_protein_groups(vault, source_id = source_id))
    groups$group_id <- as.character(groups$group_id)
    groups <- groups[match(ids, groups$group_id), , drop = FALSE]
    members <- dplyr::collect(tbl_group_members(vault, source_id = source_id))
    members$group_id <- as.character(members$group_id)
    members$protein_id <- as.character(members$protein_id)
    member_rows <- split(members, members$group_id)
    for (column in c("protein_id", "accession", "description")) {
        output <- switch(column, protein_id = "member_protein_ids",
            accession = "member_accessions", description = "member_names")
        groups[[output]] <- I(lapply(ids, function(id) {
            rows <- member_rows[[id]]
            if (is.null(rows)) character() else unique(rows[[column]])
        }))
    }
    leaders <- lapply(seq_along(ids), function(index) {
        first_group_member(groups$source_group_id[[index]], member_rows[[ids[[index]]]])
    })
    groups$leading_protein_id <- vapply(leaders, `[[`, character(1), "protein_id")
    groups$leading_protein_accession <- vapply(leaders, `[[`, character(1), "accession")
    add_mapping_columns(vault, groups, mapping_types, mapping_source_id)
}

add_precursor_group_annotations <- function(vault, features, level, source_id,
    mapping_types, mapping_source_id) {
    links <- dplyr::collect(tbl_vault(vault, "precursor_groups"))
    links$precursor_id <- as.character(links$precursor_id)
    links$group_id <- as.character(links$group_id)
    group_info <- dplyr::collect(tbl_protein_groups(vault, source_id = source_id))
    group_info$group_id <- as.character(group_info$group_id)
    links <- links[links$group_id %in% group_info$group_id, , drop = FALSE]
    group_names <- stats::setNames(group_info$source_group_id, group_info$group_id)
    link_groups <- split(links$group_id, links$precursor_id)
    members <- dplyr::collect(tbl_group_members(vault, source_id = source_id))
    members$group_id <- as.character(members$group_id)
    members$protein_id <- as.character(members$protein_id)
    member_groups <- split(members, members$group_id)
    feature_precursors <- if (level == "precursor") {
        lapply(features$precursor_id, function(id) id)
    } else {
        features$precursor_ids
    }
    feature_groups <- lapply(feature_precursors, function(ids) {
        unique(unlist(link_groups[ids], use.names = FALSE))
    })
    features$protein_group_ids <- I(feature_groups)
    features$protein_group_names <- I(lapply(feature_groups, function(ids) {
        unname(group_names[ids])
    }))
    feature_members <- lapply(feature_groups, function(ids) {
        rows <- member_groups[ids]
        rows <- rows[!vapply(rows, is.null, logical(1))]
        if (length(rows)) do.call(rbind, rows) else NULL
    })
    for (column in c("protein_id", "accession", "description")) {
        output <- switch(column, protein_id = "member_protein_ids",
            accession = "member_accessions", description = "member_names")
        features[[output]] <- I(lapply(feature_members, function(rows) {
            if (is.null(rows)) character() else unique(rows[[column]])
        }))
    }
    group_leaders <- lapply(seq_len(nrow(group_info)), function(index) {
        id <- group_info$group_id[[index]]
        first_group_member(group_info$source_group_id[[index]], member_groups[[id]])
    })
    names(group_leaders) <- group_info$group_id
    features$leading_protein_ids <- I(lapply(feature_groups, function(ids) {
        vapply(group_leaders[ids], `[[`, character(1), "protein_id")
    }))
    features$leading_protein_accessions <- I(lapply(feature_groups, function(ids) {
        vapply(group_leaders[ids], `[[`, character(1), "accession")
    }))
    features$leading_protein_id <- vapply(features$leading_protein_ids,
        function(ids) if (length(ids) == 1L) unname(ids) else NA_character_,
        character(1))
    features$leading_protein_accession <- vapply(features$leading_protein_accessions,
        function(ids) if (length(ids) == 1L) unname(ids) else NA_character_,
        character(1))
    add_mapping_columns(vault, features, mapping_types, mapping_source_id)
}

add_mapping_columns <- function(vault, features, mapping_types, mapping_source_id) {
    if (!"uniprot-idmapping" %in% vault_info(vault)$layers) return(features)
    annotations <- dplyr::collect(tbl_protein_identifiers(vault,
        source_id = mapping_source_id, mapping_type = mapping_types))
    annotations$protein_id <- as.character(annotations$protein_id)
    types <- if (is.null(mapping_types)) unique(annotations$mapping_type) else unique(mapping_types)
    columns <- paste0(rep(c("leading_", "member_"), each = length(types)),
        rep(tolower(gsub("[^[:alnum:]]+", "_", types)), times = 2L))
    if (anyDuplicated(columns) || any(columns %in% names(features))) {
        stop("Mapping types produce duplicate rowData column names", call. = FALSE)
    }
    for (type in types) {
        suffix <- tolower(gsub("[^[:alnum:]]+", "_", type))
        mapped <- split(annotations$mapped_id[annotations$mapping_type == type],
            annotations$protein_id[annotations$mapping_type == type])
        collapse_ids <- function(ids) {
            values <- unlist(mapped[unique(ids[!is.na(ids)])], use.names = FALSE)
            values <- sort(unique(values[!is.na(values) & nzchar(values)]))
            if (length(values)) paste(values, collapse = ";") else NA_character_
        }
        features[[paste0("leading_", suffix)]] <- vapply(features$leading_protein_id,
            collapse_ids, character(1))
        features[[paste0("member_", suffix)]] <- vapply(features$member_protein_ids,
            collapse_ids, character(1))
    }
    features
}

first_group_member <- function(group_name, members) {
    missing <- c(protein_id = NA_character_, accession = NA_character_)
    if (is.null(members) || is.na(group_name) || !nzchar(group_name)) return(missing)
    first_accession <- trimws(strsplit(group_name, "[;/]")[[1L]][[1L]])
    index <- match(first_accession, members$accession)
    if (is.na(index)) return(missing)
    c(protein_id = members$protein_id[[index]],
        accession = members$accession[[index]])
}
