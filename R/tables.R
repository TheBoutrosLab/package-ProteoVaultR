#' List tables and views in a vault
#'
#' Returns the names of physical tables and views in the vault's `main` schema.
#' The result reflects the open database, including optional imported layers.
#' Pass one of the returned names to [tbl_vault()] to query it.
#'
#' @inheritParams close_vault
#' @return A sorted character vector of table and view names.
#' @export
list_vault_tables <- function(vault) {
    DBI::dbGetQuery(
        conn = vault_connection(vault),
        statement = paste(
            "SELECT table_name FROM information_schema.tables",
            "WHERE table_schema = 'main'",
            "ORDER BY table_name"
        )
    )$table_name
}


#' Access a saved table or view lazily
#'
#' This escape hatch exposes the schema without loading an entire table into R.
#' Use dplyr filters and projections before collecting. Optional PIN and raw PSM
#' feature tables require their corresponding source layer.
#'
#' @inheritParams close_vault
#' @param name One table or view name in the database.
#' @return A lazy dbplyr table. The vault must remain open while it is used.
#' @export
tbl_vault <- function(vault, name) {
    connection <- vault_connection(vault = vault)
    if (length(name) != 1L || is.na(name)) {
        stop("name must identify one table or view", call. = FALSE)
    }
    layer <- switch(name,
        pin_features = "msbooster",
        psm_features = "raw-search",
        psm_metadata = "raw-search",
        qc_metrics = "qc",
        protein_identifiers = "uniprot-idmapping",
        NULL
    )
    if (!is.null(layer)) require_layer(vault = vault, layer = layer)
    readable <- c(
        psms = "psms_readable",
        protein_groups = "protein_groups_readable",
        psm_scores = "psm_scores_readable",
        protein_quant = "protein_quant_readable",
        precursor_quant = "precursor_quant_readable",
        peptide_quant = "peptide_quant_readable"
    )
    if (name %in% names(readable)) name <- readable[[name]]
    dbplyr::tbl_sql(
        subclass = "duckdb_connection",
        src = dbplyr::src_dbi(connection),
        from = dbplyr::in_schema(schema = "main", table = name)
    )
}


#' Query PSMs with their acquisition and peptide context
#'
#' Selects validated PSM observations by default. Requesting `search` requires
#' the optional raw-search layer. Joining protein memberships is deliberately
#' left separate so a shared peptide does not duplicate PSM rows.
#'
#' @inheritParams close_vault
#' @param stage One or both of `validation` and `search`.
#' @return A lazy dbplyr query containing PSM, spectrum, run, sample or TMT
#'   plex, peptide, modification, and charge information. TMT acquisitions use
#'   `plex_name` as `sample_name` because one run contains several reporters.
#'   RT remains in seconds.
#' @export
tbl_psms <- function(vault, stage = "validation") {
    stage <- match.arg(
        arg = stage,
        choices = c("validation", "search"),
        several.ok = TRUE
    )
    if ("search" %in% stage) require_layer(vault = vault, layer = "raw-search")
    spectrum <- tbl_vault(vault = vault, name = "spectra") %>%
        dplyr::select(dplyr::all_of(c(
            "spectrum_id", "run_id", "native_id", "rt_seconds", "ms_level", "precursor_mz"
        )))
    tmt_run_plex <- tbl_vault(vault, "tmt_runs") %>%
        dplyr::inner_join(tbl_vault(vault, "tmt_plexes"), by = "plex_id") %>%
        dplyr::select(dplyr::all_of(c("run_id", "plex_name")))
    tbl_vault(vault = vault, name = "psms") %>%
        dplyr::filter(.data$stage %in% .env$stage) %>%
        dplyr::inner_join(y = spectrum, by = "spectrum_id") %>%
        dplyr::inner_join(
            y = tbl_vault(vault = vault, name = "runs"),
            by = "run_id"
        ) %>%
        dplyr::left_join(
            y = tbl_vault(vault = vault, name = "samples"),
            by = "sample_id"
        ) %>%
        dplyr::left_join(y = tmt_run_plex, by = "run_id") %>%
        dplyr::mutate(sample_name = dplyr::coalesce(
            .data$sample_name, .data$plex_name)) %>%
        dplyr::select(-dplyr::all_of("plex_name")) %>%
        dplyr::inner_join(
            y = tbl_vault(vault = vault, name = "precursors"),
            by = "precursor_id"
        ) %>%
        dplyr::inner_join(
            y = tbl_vault(vault = vault, name = "peptidoforms"),
            by = "peptidoform_id"
        ) %>%
        dplyr::inner_join(
            y = tbl_vault(vault = vault, name = "peptides"),
            by = "peptide_id"
        )
}

#' Query PSM scores without losing their meaning
#'
#' Scores remain in long form with their original names and higher-better flag.
#' Multiple scores for one PSM are distinct observations, not additional PSMs.
#'
#' @inheritParams tbl_psms
#' @return A lazy dbplyr query with score name/value/orientation and PSM context.
#' @export
tbl_psm_scores <- function(vault, stage = "validation") {
    context <- tbl_psms(vault = vault, stage = stage) %>%
        dplyr::select(dplyr::all_of(c(
            "psm_id", "sample_name", "stage", "source_id", "search_engine",
            "native_id", "sequence", "modified_sequence", "charge"
        )))
    tbl_vault(vault = vault, name = "psm_scores") %>%
        dplyr::inner_join(y = context, by = "psm_id")
}


#' Query source-specific protein groups
#'
#' Quantitative DIA-NN groups and protein inference groups remain separate.
#' OpenMS protein_group and indistinguishable_proteins definitions are also
#' distinct. Filter method, type, and source before comparing definitions.
#'
#' @inheritParams close_vault
#' @param source_id Optional source identifier; NULL retains all sources.
#' @param method Optional inference/grouping method as stored in the database.
#' @param group_type Optional group definition, such as `quantification`,
#'   `protein_group`, or `indistinguishable_proteins`.
#' @return A lazy dbplyr query with group score and score_type. Multiple sources
#'   are retained when not filtered; no equivalence between them is inferred.
#' @export
tbl_protein_groups <- function(
    vault,
    source_id = NULL,
    method = NULL,
    group_type = NULL
) {
    query <- tbl_vault(vault = vault, name = "protein_group_summary")
    if (!is.null(source_id)) {
        query <- dplyr::filter(query, .data$source_id %in% .env$source_id)
    }
    if (!is.null(method)) {
        query <- dplyr::filter(query, .data$method %in% .env$method)
    }
    if (!is.null(group_type)) {
        query <- dplyr::filter(query, .data$group_type %in% .env$group_type)
    }
    query
}

#' Query the proteins belonging to each group
#'
#' Returns one row per group-member association, including original accession
#' and decoy/contaminant flags. Does not select representatives or propagate a
#' group quantity to every member protein.
#'
#' @inheritParams tbl_protein_groups
#' @return A lazy dbplyr query joining groups, memberships, and protein records.
#' @export
tbl_group_members <- function(
    vault,
    source_id = NULL,
    method = NULL,
    group_type = NULL
) {
    groups <- tbl_protein_groups(
        vault = vault,
        source_id = source_id,
        method = method,
        group_type = group_type
    )
    tbl_vault(vault = vault, name = "protein_group_members") %>%
        dplyr::inner_join(y = groups, by = "group_id") %>%
        dplyr::inner_join(
            y = tbl_vault(vault = vault, name = "proteins"),
            by = "protein_id"
        )
}


#' Query protein identifier annotations
#'
#' Returns every mapped identifier in long form with the original protein
#' accession and mapping-source provenance. A protein can have multiple values
#' for the same mapping type. The query does not pick a representative value.
#'
#' @inheritParams close_vault
#' @param source_id Optional mapping archive source identifier; NULL retains all.
#' @param mapping_type Optional UniProt mapping type, such as `Gene_Name`,
#'   `GeneID`, or `Ensembl`.
#' @return A lazy dbplyr query with protein_id, accession, source_id,
#'   source_accession, mapping_type, and mapped_id.
#' @export
tbl_protein_identifiers <- function(vault, source_id = NULL, mapping_type = NULL) {
    identifiers <- tbl_vault(vault = vault, name = "protein_identifiers")
    if (!is.null(source_id)) {
        identifiers <- dplyr::filter(identifiers, .data$source_id %in% .env$source_id)
    }
    if (!is.null(mapping_type)) {
        identifiers <- dplyr::filter(identifiers, .data$mapping_type %in% .env$mapping_type)
    }
    identifiers %>%
        dplyr::inner_join(
            y = tbl_vault(vault = vault, name = "proteins") %>%
                dplyr::select(dplyr::all_of(c("protein_id", "accession"))),
            by = "protein_id"
        )
}


#' Query protein-group quantities
#'
#' Returns source-specific group quantities with original sample names and the
#' recorded grouping/quantification methods. Does not normalize or refilter.
#'
#' @inheritParams resolve_quant_source
#' @return A lazy dbplyr query. Missing values remain SQL NULL until collected.
#' @export
tbl_protein_quant <- function(vault, source_id = NULL, method = NULL) {
    method <- resolve_quant_method(vault, "protein_group", method)
    source_id <- resolve_quant_source(
        vault = vault,
        level = "protein_group",
        source_id = source_id,
        method = method
    )
    tbl_vault(vault = vault, name = "protein_quant_long") %>%
        dplyr::filter(
            .data$source_id == .env$source_id,
            .data$quantification_method == .env$method
        )
}

#' Query precursor quantities
#'
#' Precursors are charged modified peptide forms. This function does not combine
#' charge states, modifications, or repeated TMT observations. `quant_id` and
#' source element/map IDs retain individual TMT scans or features. `rt_seconds`
#' is the TMT centroid RT or DIA-NN report RT; matrix rows have no RT. A
#' `spectrum_id` is present only for scan measurements matched to accepted PSMs.
#' Report confidence is returned where supplied.
#'
#' @inheritParams resolve_quant_source
#' @return A lazy dbplyr query with source notation and quantification method.
#' @export
tbl_precursor_quant <- function(vault, source_id = NULL, method = NULL) {
    method <- resolve_quant_method(vault, "precursor", method)
    source_id <- resolve_quant_source(
        vault = vault,
        level = "precursor",
        source_id = source_id,
        method = method
    )
    tbl_vault(vault = vault, name = "precursor_quant_long") %>%
        dplyr::filter(.data$source_id == .env$source_id, .data$method == .env$method)
}


#' Query native peptide quantities
#'
#' Returns one ProteinQuantifier source row per peptide feature and sample.
#' `feature_id` preserves source rows, including their modification, charge,
#' protein, and fraction annotations. A charge of zero means an aggregate row.
#'
#' @inheritParams resolve_quant_source
#' @return A lazy dbplyr query with source peptide and sample context.
#' @export
tbl_peptide_quant <- function(vault, source_id = NULL, method = NULL) {
    method <- resolve_quant_method(vault, "peptide", method)
    source_id <- resolve_quant_source(vault, "peptide", source_id, method)
    tbl_vault(vault, "peptide_quant_long") %>%
        dplyr::filter(.data$source_id == .env$source_id,
            .data$method == .env$method)
}


#' Query published cohort QC metrics
#'
#' Returns the final merged DIA or DDA metric table without reshaping or
#' recomputing values. Cohort rows have an empty `sample_name`; unavailable
#' metrics retain their status and have missing values. The `value` column is
#' text because some published metrics are categorical.
#'
#' @inheritParams close_vault
#' @return A lazy dbplyr table with the published metric dimensions, value,
#'   status, and source file.
#' @export
tbl_qc_metrics <- function(vault) {
    require_layer(vault = vault, layer = "qc")
    tbl_vault(vault = vault, name = "qc_metrics")
}
