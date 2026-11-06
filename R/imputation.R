#' Conditionally impute wholly missing sample bins
#'
#' Uses a per-sample limit of detection (LOD), calculated as the smallest
#' positive observed feature intensity multiplied by `lod_factor`. For each
#' feature, a wholly missing bin is filled only when another bin has at least
#' `min_observed` values. Bins with insufficient observations are masked.
#' The input assay is retained and a new assay contains the imputed values.
#'
#' @param se A SummarizedExperiment.
#' @param bin Categorical `colData` column assigning every sample to a bin.
#' @param min_observed Minimum observed values in a bin to permit imputation.
#' @param lod_factor Positive multiplier below one, default 0.5.
#' @param assay Input assay name, default `"quantity"`.
#' @param output_assay Output assay name, default `"quantity_imputed"`.
#' @return The updated SummarizedExperiment. `rowData` contains per-bin counts
#'   and statuses; `colData` contains per-sample minima and LODs.
#' @export
impute_conditional_lod <- function(se, bin, min_observed = 2L,
    lod_factor = 0.5, assay = "quantity", output_assay = "quantity_imputed") {
    if (!methods::is(se, "SummarizedExperiment"))
        stop("se must be a SummarizedExperiment", call. = FALSE)
    if (!is.character(bin) || length(bin) != 1L || is.na(bin) ||
        !bin %in% names(SummarizedExperiment::colData(se)))
        stop("bin must name a colData column", call. = FALSE)
    if (!is.numeric(min_observed) || length(min_observed) != 1L ||
        is.na(min_observed) || min_observed < 1 || min_observed != as.integer(min_observed))
        stop("min_observed must be a positive integer", call. = FALSE)
    if (!is.numeric(lod_factor) || length(lod_factor) != 1L ||
        !is.finite(lod_factor) || lod_factor <= 0 || lod_factor >= 1)
        stop("lod_factor must be between zero and one", call. = FALSE)
    if (!assay %in% SummarizedExperiment::assayNames(se) ||
        output_assay %in% SummarizedExperiment::assayNames(se))
        stop("Input assay must exist and output assay must be new", call. = FALSE)
    values <- SummarizedExperiment::assay(se, assay)
    if (!is.matrix(values) || !is.numeric(values) ||
        any(values < 0 | is.infinite(values), na.rm = TRUE))
        stop("assay must contain nonnegative finite intensities", call. = FALSE)
    groups <- as.character(SummarizedExperiment::colData(se)[[bin]])
    if (anyNA(groups) || any(!nzchar(groups)))
        stop("Every sample needs a bin", call. = FALSE)
    bins <- unique(groups)
    suffix <- make.names(bins, unique = TRUE)
    names(suffix) <- bins
    observed <- values
    observed[!is.na(observed) & observed == 0] <- NA_real_
    minima <- apply(observed, 2L, function(column) {
        positive <- column[!is.na(column)]
        if (length(positive)) min(positive) else NA_real_
    })
    lod <- minima * lod_factor
    result <- observed
    counts <- vapply(bins, function(group) {
        rowSums(!is.na(observed[, groups == group, drop = FALSE]))
    }, numeric(nrow(observed)))
    if (is.null(dim(counts))) counts <- matrix(counts, ncol = 1L)
    eligible <- rowSums(counts >= min_observed) > 0L
    row_data <- SummarizedExperiment::rowData(se)
    if (any(c("imputed_any", paste0("n_observed_", suffix),
        paste0("imputation_status_", suffix)) %in% names(row_data)))
        stop("Imputation rowData columns already exist", call. = FALSE)
    row_data$imputed_any <- rep(FALSE, nrow(se))
    for (index in seq_along(bins)) {
        columns <- which(groups == bins[[index]])
        status <- ifelse(counts[, index] >= min_observed, "observed",
            ifelse(counts[, index] == 0L & eligible, "lod_imputed", "insufficient"))
        result[status == "insufficient", columns] <- NA_real_
        fill <- which(status == "lod_imputed")
        if (length(fill) && anyNA(lod[columns]))
            stop("Cannot calculate LOD for a sample needing imputation", call. = FALSE)
        if (length(fill)) {
            result[fill, columns] <- matrix(rep(lod[columns], each = length(fill)),
                nrow = length(fill))
            row_data$imputed_any[fill] <- TRUE
        }
        row_data[[paste0("n_observed_", suffix[[index]])]] <- counts[, index]
        row_data[[paste0("imputation_status_", suffix[[index]])]] <- status
    }
    col_data <- SummarizedExperiment::colData(se)
    if (any(c("observed_minimum", "imputation_lod") %in% names(col_data)))
        stop("Imputation colData columns already exist", call. = FALSE)
    col_data$observed_minimum <- minima
    col_data$imputation_lod <- lod
    SummarizedExperiment::rowData(se) <- row_data
    SummarizedExperiment::colData(se) <- col_data
    SummarizedExperiment::assay(se, output_assay) <- result
    S4Vectors::metadata(se)$conditional_imputation <- list(bin = bin,
        min_observed = min_observed, lod_method = "sample_min",
        lod_factor = lod_factor, assay = assay, output_assay = output_assay)
    se
}
