#' Resolve a quantification source
#'
#' Resolves source IDs from saved measurement methods and the source inventory.
#' Only the small source-ID result is collected. Known empty matrix/report and
#' ProteinQuantifier formats remain selectable from provenance. Ambiguous
#' sources must be selected explicitly, rather than combined or chosen by order.
#'
#' @inheritParams close_vault
#' @param level `protein_group`, `precursor`, or `peptide`.
#' @param source_id Optional source identifier from `vault_info(vault)$sources`.
#' @param method Quantification method. `NULL` selects `matrix` when present,
#'   otherwise the sole available method (such as `ProteinQuantifier`). Specify
#'   a method when several non-matrix methods are available.
#' @return One source identifier.
#' @export
resolve_quant_source <- function(vault, level, source_id = NULL, method = NULL) {
    vault_connection(vault = vault)
    level <- match.arg(arg = level, choices = c("protein_group", "precursor", "peptide"))
    method <- resolve_quant_method(vault, level, method)
    table <- switch(level, protein_group = "protein_quant",
        precursor = "precursor_quant", peptide = "peptide_quant")
    table <- paste0(table, "_readable")
    recorded <- DBI::dbGetQuery(
        conn = vault_connection(vault = vault),
        statement = paste0("SELECT DISTINCT source_id FROM ", table, " WHERE method = ?"),
        params = list(method)
    )
    format <- if (method == "matrix") {
        if (level == "protein_group") "diann_pg_matrix" else "diann_pr_matrix"
    } else if (method == "ProteinQuantifier") {
        if (level == "protein_group") {
            "openms_protein_quant"
        } else if (level == "peptide") {
            "openms_peptide_quant"
        } else {
            character()
        }
    } else if (level == "precursor" && method %in% c(
        "Precursor.Quantity", "Precursor.Normalised", "Precursor.Translated"
    )) {
        "diann_report"
    } else {
        character()
    }
    candidates <- vault$sources[
        vault$sources$source_id %in% recorded$source_id | vault$sources$format %in% format, ,
        drop = FALSE
    ]
    sources <- candidates
    if (!is.null(source_id)) {
        if (length(source_id) != 1L || is.na(source_id)) {
            stop("source_id must select one source", call. = FALSE)
        }
        sources <- sources[sources$source_id == source_id, , drop = FALSE]
    }
    if (nrow(sources) != 1L) {
        choices <- paste(
            paste(candidates$source_id, candidates$relative_path, sep = " : "),
            collapse = "; "
        )
        stop("Expected one compatible ", level, " source for ", method,
            "; select source_id explicitly. Candidates: ", choices, call. = FALSE)
    }
    sources$source_id[[1L]]
}

resolve_quant_method <- function(vault, level, method = NULL) {
    if (!is.null(method)) {
        if (!is.character(method) || length(method) != 1L ||
            is.na(method) || !nzchar(method)) {
            stop("method must select one recorded measurement method", call. = FALSE)
        }
        return(method)
    }
    table <- switch(level, protein_group = "protein_quant_readable",
        precursor = "precursor_quant_readable", peptide = "peptide_quant_readable")
    available <- DBI::dbGetQuery(vault_connection(vault),
        paste0("SELECT DISTINCT method FROM ", table))$method
    if ("matrix" %in% available) return("matrix")
    if (length(available) == 1L) return(available[[1L]])
    if (length(available) == 0L) {
        formats <- vault$sources$format
        if (level == "protein_group" && "diann_pg_matrix" %in% formats) return("matrix")
        if (level == "precursor" && "diann_pr_matrix" %in% formats) return("matrix")
        if (level == "protein_group" && "openms_protein_quant" %in% formats) return("ProteinQuantifier")
        if (level == "peptide" && "openms_peptide_quant" %in% formats) return("ProteinQuantifier")
    }
    stop("Select a quantification method for ", level, "; available: ",
        paste(available, collapse = ", "), call. = FALSE)
}
