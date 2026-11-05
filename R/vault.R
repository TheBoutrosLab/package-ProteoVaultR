#' Open a ProteoVault database
#'
#' Opens an existing database read-only and checks schema compatibility. The
#' handle contains the connection and small metadata tables, not the full data.
#' Lazy queries remain valid only while the connection is open. Close and reopen
#' after Python adds a layer or replaces the database file.
#'
#' @param path Path to an existing ProteoVault DuckDB file.
#' @param config Named list of DuckDB configuration options, such as
#'   `list(memory_limit = "2GB", threads = "2")`.
#' @return An S3 `proteovault` connection handle.
#' @details DuckDB storage compatibility and ProteoVault schema compatibility
#'   are separate. A file using DuckDB 2.0 storage needs a compatible R engine.
#'   Nanosecond timestamps and other BIGINT columns are retained as integer64.
#' @examples
#' \dontrun{
#' vault <- open_vault(path = "results.proteovault.duckdb")
#' vault_info(vault)
#' close_vault(vault)
#' }
#' @export
open_vault <- function(path, config = list()) {
    path <- normalizePath(path = path, mustWork = TRUE)
    connection <- tryCatch(
        DBI::dbConnect(
            drv = duckdb::duckdb(shared_home = FALSE),
            dbdir = path,
            read_only = TRUE,
            config = config,
            bigint = "integer64"
        ),
        error = function(error) {
            stop(
                "Cannot open ProteoVault database. Check that the R DuckDB engine ",
                "can read its storage format; DuckDB 2.0 files need a compatible ",
                "reader. Original error: ",
                conditionMessage(error),
                call. = FALSE
            )
        }
    )
    complete <- FALSE
    on.exit(
        if (!complete) DBI::dbDisconnect(conn = connection),
        add = TRUE
    )
    metadata <- check_schema(connection = connection)
    vault <- structure(
        list(
            connection = connection,
            path = path,
            metadata = metadata,
            analysis = tibble::as_tibble(DBI::dbGetQuery(
                conn = connection,
                statement = "SELECT * FROM analysis"
            )),
            sources = tibble::as_tibble(DBI::dbGetQuery(
                conn = connection,
                statement = "SELECT * FROM source_files ORDER BY relative_path"
            )),
            sample_metadata = NULL
        ),
        class = "proteovault"
    )
    if (!bit64::is.integer64(vault$sources$mtime_ns)) {
        stop("Source timestamps must retain their 64-bit representation", call. = FALSE)
    }
    complete <- TRUE
    vault
}

#' Attach sample metadata to a vault handle
#'
#' Metadata is stored only in the returned R handle. The DuckDB database is
#' never modified. Match rows using the original `sample_name`; an optional
#' `sample_id` column supplies shorter labels for quantity distribution plots.
#' Samples without a metadata row retain their original plot label.
#'
#' @inheritParams close_vault
#' @param metadata A data frame with a unique, nonmissing `sample_name` column
#'   matching names in the vault. Other columns are retained in memory.
#' @return A vault handle with attached sample metadata. Assign the result.
#' @export
set_sample_metadata <- function(vault, metadata) {
    connection <- vault_connection(vault)
    if (!is.data.frame(metadata) || !"sample_name" %in% names(metadata)) {
        stop("metadata must be a data frame with a sample_name column", call. = FALSE)
    }
    sample_names <- metadata$sample_name
    if (!is.character(sample_names) || anyNA(sample_names) || any(!nzchar(sample_names)) || anyDuplicated(sample_names)) {
        stop("metadata sample_name must contain unique, nonempty strings", call. = FALSE)
    }
    known <- DBI::dbGetQuery(connection, "SELECT sample_name FROM samples")$sample_name
    if (!all(sample_names %in% known)) {
        stop("metadata sample_name contains names absent from the vault", call. = FALSE)
    }
    if ("sample_id" %in% names(metadata)) {
        labels <- metadata$sample_id
        if (!is.character(labels) || anyNA(labels) || any(!nzchar(labels)) || anyDuplicated(labels)) {
            stop("metadata sample_id must contain unique, nonempty strings", call. = FALSE)
        }
        displayed <- c(labels, setdiff(known, sample_names))
        if (anyDuplicated(displayed)) {
            stop("metadata sample_id conflicts with an unlabeled sample_name", call. = FALSE)
        }
    }
    vault$sample_metadata <- tibble::as_tibble(metadata)
    vault
}

#' Read sample metadata attached to a vault handle
#'
#' @inheritParams close_vault
#' @return The in-memory metadata table, or NULL when none is attached.
#' @export
sample_metadata <- function(vault) {
    vault_connection(vault)
    vault$sample_metadata
}

vault_sample_labels <- function(vault) {
    metadata <- sample_metadata(vault)
    if (is.null(metadata) || !"sample_id" %in% names(metadata)) return(NULL)
    stats::setNames(metadata$sample_id, metadata$sample_name)
}

#' Close a ProteoVault connection
#'
#' Closing an already closed handle is harmless. Collected R results remain
#' usable, while lazy queries require the connection to stay open.
#'
#' @param vault A handle returned by [open_vault()].
#' @return The handle, invisibly.
#' @export
close_vault <- function(vault) {
    if (!inherits(vault, "proteovault")) {
        stop("vault must be a proteovault handle", call. = FALSE)
    }
    if (DBI::dbIsValid(dbObj = vault$connection)) {
        DBI::dbDisconnect(conn = vault$connection)
    }
    invisible(vault)
}

#' Access the DBI connection owned by a vault
#'
#' Provides a connection for direct read-only DBI queries and checks that the
#' handle is still open. Does not repeat schema or full-data validation.
#'
#' @inheritParams close_vault
#' @return A valid DuckDB DBI connection.
#' @export
vault_connection <- function(vault) {
    if (!inherits(vault, "proteovault")) {
        stop("vault must be a proteovault handle", call. = FALSE)
    }
    if (!DBI::dbIsValid(dbObj = vault$connection)) {
        stop("Vault connection is closed; reopen it with open_vault()", call. = FALSE)
    }
    vault$connection
}

#' Print a ProteoVault handle
#'
#' Displays small metadata without counting or collecting fact tables.
#'
#' @param x A `proteovault` handle.
#' @param ... Additional arguments, unused.
#' @return The handle, invisibly.
#' @export
print.proteovault <- function(x, ...) {
    state <- if (DBI::dbIsValid(dbObj = x$connection)) "open" else "closed"
    cat("<ProteoVault: ", state, ">\n", sep = "")
    cat("Dataset: ", paste(x$analysis$dataset_id, collapse = ", "), "\n", sep = "")
    cat("Schema: ", x$metadata[["schema_version"]], "\n", sep = "")
    if (state == "open") {
        cat("Layers: ", paste(vault_info(vault = x)$layers, collapse = ", "), "\n", sep = "")
    }
    cat("Path: ", x$path, "\n", sep = "")
    invisible(x)
}


#' Describe the supported database schema
#'
#' @return A named list of required tables/views and the columns used by the
#'   reader. The Python package owns the schema definition.
#' @export
vault_schema <- function() {
    list(
        metadata = c("key", "value"),
        analysis = c("dataset_id", "pipeline", "pipeline_version", "root_path"),
        source_files = c("source_id", "relative_path", "format", "stage"),
        samples = c("sample_id", "sample_name"),
        runs = c("run_id", "sample_id", "run_name", "file_name"),
        spectra = c("spectrum_id", "run_id", "native_id", "rt_seconds", "ms_level", "precursor_mz"),
        peptides = c("peptide_id", "sequence"),
        peptidoforms = c("peptidoform_id", "peptide_id", "modified_sequence", "notation"),
        precursors = c("precursor_id", "peptidoform_id", "charge"),
        psms = c("psm_id", "spectrum_id", "precursor_id", "stage_id", "source_id", "search_engine"),
        psms_readable = c("psm_id", "spectrum_id", "precursor_id", "stage", "source_id", "search_engine"),
        psm_scores = c("psm_id", "score_id", "value", "higher_better"),
        psm_scores_readable = c("psm_id", "name", "value", "higher_better"),
        psm_features = c("psm_id", "name", "value"),
        psm_metadata = c("psm_id", "name", "value"),
        proteins = c("protein_id", "accession", "is_decoy", "is_contaminant"),
        protein_identifiers = c("protein_id", "source_id", "source_accession", "mapping_type", "mapped_id"),
        protein_groups = c(
            "group_id", "source_id", "source_group_id", "method_id",
            "group_type", "score", "score_type"
        ),
        protein_groups_readable = c(
            "group_id", "source_id", "source_group_id", "method",
            "group_type", "score", "score_type"
        ),
        protein_group_members = c("group_id", "protein_id"),
        protein_quant = c("group_id", "channel_id", "source_id", "method_id", "quantity"),
        protein_quant_readable = c("group_id", "channel_id", "source_id", "method", "quantity"),
        precursor_quant = c("precursor_id", "channel_id", "source_id", "method_id", "quantity"),
        precursor_quant_readable = c("precursor_id", "channel_id", "source_id", "method", "quantity"),
        quant_methods = c("method_id", "name"),
        grouping_methods = c("method_id", "name"),
        processing_stages = c("stage_id", "name"),
        score_types = c("score_id", "name"),
        score_aliases = c("score_id", "source_name"),
        protein_quant_long = c(
            "sample_name", "group_id", "source_group_id", "grouping_method",
            "quantification_method", "quantity", "source_id"
        ),
        precursor_quant_long = c(
            "sample_name", "precursor_id", "modified_sequence", "notation", "charge",
            "method", "quantity", "q_value", "pep", "source_id"
        ),
        protein_group_summary = c("group_id", "method", "group_type", "member_count"),
        sample_protein_counts = c("sample_name", "quantified_protein_groups"),
        pin_features = c("source_id", "spectrum_reference", "scan_number", "feature_name", "value"),
        qc_metrics = c(
            "source_id", "row_number", "dataset_id", "layer", "sample_name",
            "metric_name", "group", "value", "status", "source_file"
        )
    )
}
#' Check ProteoVault schema compatibility
#'
#' Checks the format identifier, schema version, and required columns. This is a
#' metadata check, not a replacement for Python's full database validation.
#'
#' @param connection A valid DBI connection to a ProteoVault database.
#' @return A named character vector of database metadata.
#' @export
check_schema <- function(connection) {
    columns <- DBI::dbGetQuery(
        conn = connection,
        statement = paste(
            "SELECT table_name, column_name FROM information_schema.columns",
            "WHERE table_schema = 'main'"
        )
    )
    expected <- vault_schema()
    for (table in names(expected)) {
        actual <- columns$column_name[columns$table_name == table]
        missing <- setdiff(x = expected[[table]], y = actual)
        if (length(missing)) {
            stop("Missing schema columns in ", table, ": ", paste(missing, collapse = ", "), call. = FALSE)
        }
    }
    rows <- DBI::dbGetQuery(
        conn = connection,
        statement = "SELECT key, value FROM metadata"
    )
    metadata <- stats::setNames(
        object = rows$value,
        nm = rows$key
    )
    if (!identical(unname(metadata["format"]), "ProteoVault")) {
        stop("Database is not identified as ProteoVault", call. = FALSE)
    }
    if (!identical(unname(metadata["schema_version"]), "1.1.0")) {
        stop("Unsupported ProteoVault schema version: ", metadata["schema_version"], call. = FALSE)
    }
    metadata
}


#' Inspect database metadata and imported layers
#'
#' Returns only small metadata and source-inventory tables. Does not count or
#' collect PSM or quantity tables. An imported source with zero rows still marks
#' its layer as present.
#'
#' @inheritParams close_vault
#' @return A list containing path, metadata, analysis, sources, layers, and the
#'   R DuckDB package version.
#' @export
vault_info <- function(vault) {
    vault_connection(vault = vault)
    stage_layers <- c(
        quantitation = "core",
        `quant-details` = "quant-details",
        validation = "core",
        inference = "core",
        search = "raw-search",
        MSBooster = "msbooster",
        qc = "qc",
        reference = "fasta",
        annotation = "uniprot-idmapping"
    )
    layers <- unname(stage_layers[vault$sources$stage])
    list(
        path = vault$path,
        metadata = vault$metadata,
        analysis = vault$analysis,
        sources = vault$sources,
        layers = sort(unique(layers[!is.na(layers)])),
        reader_duckdb_version = as.character(utils::packageVersion("duckdb"))
    )
}

#' Require an imported optional layer
#'
#' @inheritParams close_vault
#' @param layer Layer name, such as `raw-search` or `msbooster`.
#' @return The handle, invisibly, or an actionable error if the layer is absent.
#' @export
require_layer <- function(vault, layer) {
    if (!layer %in% vault_info(vault = vault)$layers) {
        if (identical(layer, "uniprot-idmapping")) {
            stop(
                "Layer 'uniprot-idmapping' is not imported. Add it with: ",
                "proteovault add-uniprot-idmapping ", shQuote(vault$path),
                " --mapping /path/to/UP000005640_9606.idmapping.gz",
                call. = FALSE
            )
        }
        stop(
            "Layer '", layer, "' is not imported. Add it with: proteovault add-layer ",
            shQuote(vault$path), " --layer ", layer,
            call. = FALSE
        )
    }
    invisible(vault)
}
