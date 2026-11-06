#' Resolve a quantification source
#'
#' Resolves source IDs from saved measurement methods and the source inventory.
#' Only the small source-ID result is collected. Known empty matrix/report
#' formats remain selectable from provenance. Ambiguous
#' sources must be selected explicitly, rather than combined or chosen by order.
#'
#' @inheritParams close_vault
#' @param level `protein_group` or `precursor`.
#' @param source_id Optional source identifier from `vault_info(vault)$sources`.
#' @param method Quantification method. `matrix` selects a saved matrix;
#'   precursor reports also support `Precursor.Quantity`, `Precursor.Normalised`,
#'   and `Precursor.Translated`.
#' @return One source identifier.
#' @export
resolve_quant_source <- function(vault, level, source_id = NULL, method = "matrix") {
    vault_connection(vault = vault)
    level <- match.arg(arg = level, choices = c("protein_group", "precursor"))
    if (length(method) != 1L || is.na(method) || !nzchar(method)) {
        stop("method must select one recorded measurement method", call. = FALSE)
    }
    table <- if (level == "protein_group") "protein_quant" else "precursor_quant"
    table <- paste0(table, "_readable")
    recorded <- DBI::dbGetQuery(
        conn = vault_connection(vault = vault),
        statement = paste0("SELECT DISTINCT source_id FROM ", table, " WHERE method = ?"),
        params = list(method)
    )
    format <- if (method == "matrix") {
        if (level == "protein_group") "diann_pg_matrix" else "diann_pr_matrix"
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
            paste(candidates$source_id,
                candidates$relative_path,
                sep = " : "
            ),
            collapse = "; "
        )
        stop(
            "Expected one compatible ", level, " source for ", method,
            "; select source_id explicitly. Candidates: ", choices,
            call. = FALSE
        )
    }
    sources$source_id[[1L]]
}
