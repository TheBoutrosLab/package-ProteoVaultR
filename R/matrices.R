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
#' selecting the observation with the largest sum of reporter intensities in
#' each plex. All channels in a plex come from the same observation. The
#' `"best_across_plexes"` choice selects one linked feature for all plexes and
#' requires FeatureBased data; an unobserved plex remains NA. These choices
#' require intensity-scale quantities, not log2 ratios. The vault must record
#' the TMT quantity scale during import. Use
#' [tbl_precursor_quant()] to inspect the original observations.
#'
#' @inheritParams tbl_precursor_quant
#' @inheritParams matrix_quant
#' @param duplicates For TMT, `"best_per_plex"` (default) selects the highest
#'   total reporter intensity per precursor and plex, `"sum"` adds observations,
#'   `"best_across_plexes"` selects one linked feature by its total across
#'   plexes, and `"error"` rejects duplicates. For other methods, `"best"` is
#'   the default; `"sum"`, `"mean"`, and `"error"` are also available.
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
        if (missing(duplicates)) duplicates <- "best_per_plex"
        duplicates <- match.arg(duplicates,
            c("best_per_plex", "sum", "best_across_plexes", "error"))
        if (duplicates %in% c("best_per_plex", "best_across_plexes")) {
            return(matrix_tmt_precursor_abundance(
                vault, source_id, method, samples, duplicates))
        }
        if (duplicates == "sum") validate_tmt_intensity_scale(vault, source_id, method)
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

matrix_tmt_precursor_abundance <- function(vault, source_id, method, samples, duplicates) {
    if (duplicates == "best_across_plexes" && method != "TMT.FeatureBased") {
        stop("best_across_plexes requires TMT.FeatureBased linked features", call. = FALSE)
    }
    validate_tmt_intensity_scale(vault, source_id, method)
    query <- tbl_precursor_quant(vault, source_id, method)
    available <- dplyr::collect(dplyr::distinct(query, .data$sample_name))$sample_name
    if (is.null(samples)) samples <- sort(available)
    if (anyNA(samples) || anyDuplicated(samples) || !all(samples %in% available)) {
        stop("samples must select unique sample names present in this quantification source",
            call. = FALSE)
    }
    values <- query %>%
        dplyr::select(dplyr::all_of(c("precursor_id", "sample_name",
            "source_element_id", "quantity"))) %>%
        dplyr::collect()
    if (anyNA(values$source_element_id) || any(!nzchar(values$source_element_id))) {
        stop("TMT observation selection requires source_element_id", call. = FALSE)
    }
    if (any(is.finite(values$quantity) & values$quantity < 0, na.rm = TRUE)) {
        stop("TMT observation selection requires nonnegative intensity-scale quantities",
            call. = FALSE)
    }
    sample_plex <- DBI::dbGetQuery(vault_connection(vault), paste(
        "SELECT s.sample_name, p.plex_name FROM samples s",
        "JOIN quant_channels c USING(sample_id)",
        "JOIN tmt_channels t USING(channel_id)",
        "JOIN tmt_plexes p USING(plex_id)"
    ))
    values$plex <- sample_plex$plex_name[match(values$sample_name, sample_plex$sample_name)]
    if (anyNA(values$plex)) {
        stop("TMT observation selection requires a plex for every selected sample",
            call. = FALSE)
    }
    observations <- values %>%
        dplyr::group_by(.data$precursor_id, .data$plex, .data$source_element_id) %>%
        dplyr::summarise(score = if (all(is.na(.data$quantity))) NA_real_ else
            sum(.data$quantity, na.rm = TRUE), .groups = "drop")
    if (duplicates == "best_across_plexes") {
        scores <- observations %>%
            dplyr::group_by(.data$precursor_id, .data$source_element_id) %>%
            dplyr::summarise(score = if (all(is.na(.data$score))) NA_real_ else
                sum(.data$score, na.rm = TRUE), .groups = "drop")
        keys <- c("precursor_id", "source_element_id")
    } else {
        scores <- observations
        keys <- c("precursor_id", "plex", "source_element_id")
    }
    scores$rank_score <- ifelse(is.na(scores$score), -Inf, scores$score)
    group_keys <- setdiff(keys, "source_element_id")
    chosen <- scores %>%
        dplyr::arrange(dplyr::desc(.data$rank_score), .data$source_element_id) %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(group_keys))) %>%
        dplyr::slice_head(n = 1L) %>%
        dplyr::ungroup() %>%
        dplyr::select(dplyr::all_of(keys))
    selected <- dplyr::inner_join(values, chosen, by = keys) %>%
        dplyr::filter(.data$sample_name %in% .env$samples)
    repeated <- duplicated(selected[c("precursor_id", "sample_name")])
    if (any(repeated)) {
        stop("Selected TMT observation has duplicate precursor/sample quantities",
            call. = FALSE)
    }
    rows <- sort(unique(values$precursor_id[values$sample_name %in% samples]))
    result <- matrix(NA_real_, length(rows), length(samples),
        dimnames = list(as.character(rows), samples))
    if (nrow(selected)) {
        result[cbind(match(selected$precursor_id, rows),
            match(selected$sample_name, samples))] <- selected$quantity
    }
    result
}

tmt_quantity_scale <- function(vault, source_id, method) {
    source_id <- resolve_quant_source(vault, "precursor", source_id, method)
    recorded <- DBI::dbGetQuery(vault_connection(vault),
        "SELECT value FROM metadata WHERE key = ?",
        params = list(paste0("tmt_quantity_scale:", as.character(source_id))))
    if (nrow(recorded) != 1L) {
        stop("TMT quantity scale was not recorded; rebuild the vault", call. = FALSE)
    }
    recorded$value[[1L]]
}

validate_tmt_intensity_scale <- function(vault, source_id, method) {
    scale <- tmt_quantity_scale(vault, source_id, method)
    if (scale == "log2_ratio") {
        stop("TMT aggregation requires intensity-scale quantities, not log2 ratios",
            call. = FALSE)
    }
    if (method == "TMT.FeatureBased" && scale != "quasi_intensity") {
        warning("FeatureBased TMT uses intensity-scale values without quasi-intensity normalization",
            call. = FALSE)
    }
    invisible(scale)
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
