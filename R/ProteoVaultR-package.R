#' Read and analyze linked ProteoVault results
#'
#' Open a read-only ProteoVault database with [open_vault()], query lazy result
#' tables, and explicitly collect results or construct quantitative matrices.
#' Core PSM access selects validated observations by default. Source identifiers,
#' inference group definitions, and quantification methods remain explicit.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data .env
#' @keywords internal
"_PACKAGE"

# Columns used in ggplot2's nonstandard evaluation.
utils::globalVariables(c(
    "block", "charge", "column", "group", "log2_quantity",
    "modification", "p05", "p25", "p50", "p75", "p95",
    "placement", "quantity_type", "residue", "sample_label",
    "sample_name", "status", "value"
))
