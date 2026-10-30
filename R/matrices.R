#' Materialize a quantitative matrix
#'
#' Filters to one source and method before retrieving feature/sample quantities.
#' Rejects duplicate feature/sample records instead of summing them. Missing
#' combinations become NA; explicit zeros remain zero. This operation allocates
#' an R matrix and is not lazy.
#'
#' @param query A lazy quantitative query with sample_name and quantity columns.
#' @param feature Name of the feature identifier column.
#' @param samples Optional character vector selecting and ordering sample names.
#' @return A numeric feature-by-sample matrix with original identifiers.
#' @export
matrix_quant <- function(query, feature, samples = NULL) {
    available <- query %>%
        dplyr::distinct(.data$sample_name) %>%
        dplyr::collect()
    if (is.null(samples)) samples <- sort(available$sample_name)
    if (anyNA(samples) || anyDuplicated(samples) || !all(samples %in% available$sample_name)) {
        stop("samples must select unique sample names present in this quantification source", call. = FALSE)
    }
    selected <- query %>%
        dplyr::filter(.data$sample_name %in% .env$samples) %>%
        dplyr::select(dplyr::all_of(c(feature, "sample_name", "quantity")))
    duplicates <- selected %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(c(feature, "sample_name")))) %>%
        dplyr::summarise(records = dplyr::n(), .groups = "drop") %>%
        dplyr::filter(.data$records > 1L) %>%
        utils::head(n = 1L) %>%
        dplyr::collect()
    if (nrow(duplicates)) {
        stop("Duplicate feature/sample quantities; select a single source and method", call. = FALSE)
    }
    values <- dplyr::collect(selected)
    rows <- sort(unique(values[[feature]]))
    result <- matrix(
        data = NA_real_,
        nrow = length(rows),
        ncol = length(samples),
        dimnames = list(as.character(rows), samples)
    )
    if (nrow(values)) {
        result[cbind(match(values[[feature]], rows), match(values$sample_name, samples))] <- values$quantity
    }
    result
}

#' Create a protein-group abundance matrix
#'
#' Row names are database group IDs; access annotations with [tbl_protein_groups()]
#' and [tbl_group_members()]. Quantities are not transformed or normalized.
#'
#' @inheritParams tbl_protein_quant
#' @inheritParams matrix_quant
#' @return A numeric group-by-sample matrix, collected into R memory.
#' @export
matrix_protein_abundance <- function(
    vault,
    source_id = NULL,
    method = "matrix",
    samples = NULL
) {
    matrix_quant(
        query = tbl_protein_quant(
            vault = vault,
            source_id = source_id,
            method = method
        ),
        feature = "group_id",
        samples = samples
    )
}

#' Create a charged-peptidoform abundance matrix
#'
#' Does not collapse modifications or charge states into peptide quantities.
#'
#' @inheritParams tbl_precursor_quant
#' @inheritParams matrix_quant
#' @return A numeric precursor-by-sample matrix, collected into R memory.
#' @export
matrix_precursor_abundance <- function(
    vault,
    source_id = NULL,
    method = "matrix",
    samples = NULL
) {
    matrix_quant(
        query = tbl_precursor_quant(
            vault = vault,
            source_id = source_id,
            method = method
        ),
        feature = "precursor_id",
        samples = samples
    )
}
