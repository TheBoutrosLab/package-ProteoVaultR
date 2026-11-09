#' Count quantified features per sample
#'
#' Quantified means quantity greater than zero. Missing and zero quantities do
#' not contribute. Features are counted once per sample even if several channels
#' refer to that sample. Samples represented only by missing/zero rows retain a count
#' of zero. The aggregation runs in DuckDB; only the small summary is collected.
#'
#' @inheritParams resolve_quant_source
#' @return A tibble with sample_name, quantified_count, source_id, method, and
#'   level. Count values use the numeric representation returned by the driver.
#' @export
summarize_quantified_counts <- function(
    vault,
    level = "protein_group",
    source_id = NULL,
    method = NULL
) {
    level <- match.arg(arg = level, choices = c("protein_group", "precursor", "peptide"))
    method <- resolve_quant_method(vault, level, method)
    source_id <- resolve_quant_source(
        vault = vault,
        level = level,
        source_id = source_id,
        method = method
    )
    query <- if (level == "protein_group") {
        tbl_protein_quant(
            vault = vault,
            source_id = source_id,
            method = method
        )
    } else if (level == "precursor") {
        tbl_precursor_quant(
            vault = vault,
            source_id = source_id,
            method = method
        )
    } else {
        tbl_peptide_quant(vault, source_id, method)
    }
    feature <- switch(level, protein_group = "group_id", precursor = "precursor_id",
        peptide = "feature_id")
    query %>%
        dplyr::group_by(.data$sample_name) %>%
        dplyr::summarise(
            quantified_count = dplyr::n_distinct(
                dplyr::if_else(
                    condition = !is.na(.data$quantity) & .data$quantity > 0,
                    true = as.character(.data[[feature]]),
                    false = NA_character_
                ),
                na.rm = TRUE
            ),
            .groups = "drop"
        ) %>%
        dplyr::arrange(.data$sample_name) %>%
        dplyr::collect() %>%
        dplyr::mutate(
            source_id = source_id,
            method = method,
            level = level
        )
}


#' Measure feature coverage across samples
#'
#' Counts the distinct samples with positive quantity for each feature. Retains
#' features with zero positive samples and reports the total number of samples
#' represented in the selected measurement. All aggregation runs in DuckDB.
#'
#' @inheritParams summarize_quantified_counts
#' @return A tibble containing feature_id, n_samples, total_samples, source_id,
#'   method, and level. Use n_samples greater than zero for the cohort union and
#'   n_samples equal to total_samples for the intersection.
#' @export
summarize_quantification_coverage <- function(
    vault,
    level = "protein_group",
    source_id = NULL,
    method = NULL
) {
    level <- match.arg(arg = level, choices = c("protein_group", "precursor", "peptide"))
    method <- resolve_quant_method(vault, level, method)
    source_id <- resolve_quant_source(
        vault = vault,
        level = level,
        source_id = source_id,
        method = method
    )
    query <- if (level == "protein_group") {
        tbl_protein_quant(vault = vault, source_id = source_id, method = method)
    } else if (level == "precursor") {
        tbl_precursor_quant(vault = vault, source_id = source_id, method = method)
    } else {
        tbl_peptide_quant(vault, source_id, method)
    }
    feature <- switch(level, protein_group = "group_id", precursor = "precursor_id",
        peptide = "feature_id")
    total <- query %>%
        dplyr::summarise(total_samples = dplyr::n_distinct(.data$sample_name)) %>%
        dplyr::collect()
    query %>%
        dplyr::group_by(feature_id = .data[[feature]]) %>%
        dplyr::summarise(
            n_samples = dplyr::n_distinct(
                dplyr::if_else(
                    condition = !is.na(.data$quantity) & .data$quantity > 0,
                    true = .data$sample_name,
                    false = NA_character_
                ),
                na.rm = TRUE
            ),
            .groups = "drop"
        ) %>%
        dplyr::collect() %>%
        dplyr::mutate(
            total_samples = total$total_samples[[1L]],
            source_id = source_id,
            method = method,
            level = level
        )
}
