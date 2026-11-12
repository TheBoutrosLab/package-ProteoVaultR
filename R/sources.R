#' Resolve a quantification source
#'
#' Resolves source IDs from the saved source-method inventory. Methods from
#' empty quantification files remain selectable. Ambiguous sources must be
#' selected explicitly, rather than combined or chosen by order.
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
    recorded <- DBI::dbGetQuery(
        conn = vault_connection(vault = vault),
        statement = paste(
            "SELECT q.source_id FROM source_quant_methods q",
            "JOIN quant_methods m USING(method_id)",
            "WHERE q.level = ? AND m.name = ?"
        ),
        params = list(level, method)
    )
    candidates <- vault$sources[vault$sources$source_id %in% recorded$source_id, , drop = FALSE]
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
    available <- DBI::dbGetQuery(vault_connection(vault),
        paste("SELECT DISTINCT m.name FROM source_quant_methods q",
            "JOIN quant_methods m USING(method_id) WHERE q.level = ?"),
        params = list(level))$name
    if ("matrix" %in% available) return("matrix")
    if (length(available) == 1L) return(available[[1L]])
    stop("Select a quantification method for ", level, "; available: ",
        paste(available, collapse = ", "), call. = FALSE)
}
