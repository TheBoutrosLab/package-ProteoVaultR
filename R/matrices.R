#' Materialize a quantitative matrix
#'
#' Retrieves one source and method before constructing a feature-by-sample
#' matrix. Missing combinations remain NA and explicit zeros remain zero.
#' Duplicate resolution takes place within each feature/sample pair.
#'
#' @param query A lazy quantitative query with sample_name and quantity columns.
#' @param feature Name of the feature identifier column.
#' @param samples Optional character vector selecting and ordering sample names.
#' @param duplicates How to resolve repeated feature/sample measurements:
#'   `"best"` (default) keeps the largest nonmissing quantity, `"sum"` sums
#'   nonmissing quantities, `"mean"` averages them, and `"error"` rejects
#'   duplicates. An all-missing pair remains NA for each numeric choice.
#' @return A numeric feature-by-sample matrix with the selected identifiers.
#' @export
matrix_quant <- function(
    query,
    feature,
    samples = NULL,
    duplicates = c("best", "error", "sum", "mean")
) {
    duplicates <- match.arg(arg = duplicates)
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
    if (duplicates == "error") {
        repeated <- selected %>%
            dplyr::group_by(dplyr::across(dplyr::all_of(c(feature, "sample_name")))) %>%
            dplyr::summarise(records = dplyr::n(), .groups = "drop") %>%
            dplyr::filter(.data$records > 1L) %>%
            utils::head(n = 1L) %>%
            dplyr::collect()
        if (nrow(repeated)) {
            stop("Duplicate feature/sample quantities", call. = FALSE)
        }
    }
    values <- dplyr::collect(selected)
    if (duplicates != "error") {
        duplicate_ids <- unique(values[[feature]][duplicated(values[c(feature, "sample_name")])])
        message(sprintf(
            "Collapsed %d duplicated %s ID(s) using %s",
            length(duplicate_ids), feature, duplicates
        ))
    }
    if (nrow(values) && duplicates != "error") {
        values <- values %>%
            dplyr::group_by(dplyr::across(dplyr::all_of(c(feature, "sample_name")))) %>%
            dplyr::summarise(
                quantity = aggregate_duplicate_quantities(
                    values = .data$quantity,
                    duplicates = duplicates
                ),
                .groups = "drop"
            )
    }
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

aggregate_duplicate_quantities <- function(values, duplicates) {
    if (all(is.na(values))) return(NA_real_)
    values <- values[!is.na(values)]
    switch(
        duplicates,
        best = max(values),
        sum = sum(values),
        mean = mean(values)
    )
}

#' Create a protein-group abundance matrix
#'
#' By default row names use the original source-group labels. Use
#' `row_id = "group_id"` for vault-internal numeric group IDs. Quantities are
#' not transformed or normalized. Repeated source-group labels are resolved
#' separately for each sample according to `duplicates`.
#'
#' @inheritParams tbl_protein_quant
#' @inheritParams matrix_quant
#' @param row_id Row identifier: `"source_group_id"` (default) or `"group_id"`.
#' @return A numeric group-by-sample matrix, collected into R memory.
#' @export
matrix_protein_abundance <- function(
    vault,
    source_id = NULL,
    method = NULL,
    samples = NULL,
    row_id = c("source_group_id", "group_id"),
    duplicates = c("best", "error", "sum", "mean")
) {
    row_id <- match.arg(arg = row_id)
    matrix_quant(
        query = tbl_protein_quant(
            vault = vault,
            source_id = source_id,
            method = method
        ),
        feature = row_id,
        samples = samples,
        duplicates = duplicates
    )
}

#' Create a charged-peptidoform abundance matrix
#'
#' Row names are internal precursor IDs. This function does not collapse
#' modifications or charge states into peptide quantities. TMT can have several
#' scan or feature observations for one precursor and channel. TMT defaults to
#' summing observations across every reporter channel. `"best"` and `"mean"`
#' are rejected for TMT because they can combine channels from different
#' observations or use different denominators. Use [tbl_precursor_quant()] to
#' inspect each original observation and its retention time.
#'
#' @inheritParams tbl_precursor_quant
#' @inheritParams matrix_quant
#' @param duplicates For TMT, `"sum"` is the default and `"error"` is the
#'   only other choice. For other methods, `"best"` is the default; `"sum"`,
#'   `"mean"`, and `"error"` are also available.
#' @return A numeric precursor-by-sample matrix, collected into R memory.
#' @export
matrix_precursor_abundance <- function(
    vault,
    source_id = NULL,
    method = NULL,
    samples = NULL,
    duplicates = c("best", "error", "sum", "mean")
) {
    method <- resolve_quant_method(vault, "precursor", method)
    if (startsWith(method, "TMT.")) {
        if (missing(duplicates)) duplicates <- "sum"
        duplicates <- match.arg(duplicates)
        if (!duplicates %in% c("sum", "error")) {
            stop("TMT precursor duplicates must be 'sum' or 'error'", call. = FALSE)
        }
    }
    matrix_quant(
        query = tbl_precursor_quant(
            vault = vault,
            source_id = source_id,
            method = method
        ),
        feature = "precursor_id",
        samples = samples,
        duplicates = duplicates
    )
}


#' Create a native peptide abundance matrix
#'
#' Rows use source feature IDs, retaining separate modified peptide, charge,
#' protein, and fraction records. Use [tbl_peptide_quant()] to inspect them.
#'
#' @inheritParams tbl_peptide_quant
#' @inheritParams matrix_quant
#' @return A numeric peptide-feature-by-sample matrix.
#' @export
matrix_peptide_abundance <- function(
    vault,
    source_id = NULL,
    method = NULL,
    samples = NULL,
    duplicates = c("best", "error", "sum", "mean")
) {
    matrix_quant(
        query = tbl_peptide_quant(vault, source_id, method),
        feature = "feature_id",
        samples = samples,
        duplicates = duplicates
    )
}
