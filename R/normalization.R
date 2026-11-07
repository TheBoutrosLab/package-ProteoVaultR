#' Normalize using geNorm-style stable reference proteins
#'
#' Selects abundant proteins observed in every sample without using sample
#' labels. Repeatedly removes the candidate with the largest mean pairwise SD
#' of log2 ratios, then centers each sample by the median abundance of the
#' selected references. The normalized assay remains on the input intensity
#' scale. Supply `exclude` before inspecting target results.
#'
#' @param se A [SummarizedExperiment::SummarizedExperiment].
#' @param assay Input assay containing positive, untransformed intensities.
#' @param output_assay Name of the new normalized assay.
#' @param exclude Feature IDs (row names) or `source_group_id` values excluded
#'   from reference selection. Every supplied value must match a feature.
#' @param candidate_count Number of most abundant complete proteins considered.
#' @param reference_count Number of references to select.
#' @return The updated SummarizedExperiment, with a normalized assay, per-sample
#'   log2 offsets in `colData`, selection flags in `rowData`, and settings in
#'   `metadata`.
#' @export
normalize_stable_reference_genorm <- function(se, assay = "quantity",
    output_assay = "quantity_genorm_normalized", exclude = character(),
    candidate_count = 1000L, reference_count = 200L) {
    normalize_stable_reference(se, "genorm", assay, output_assay, exclude,
        candidate_count, reference_count)
}

#' Normalize using PCA-selected stable reference proteins
#'
#' Candidate proteins are selected without sample labels. PCA identifies a
#' shared run pattern; references have a slope near the median slope and low
#' residual variation around that pattern. Sample offsets are calculated from
#' the selected references, as for [normalize_stable_reference_genorm()].
#'
#' @inheritParams normalize_stable_reference_genorm
#' @param output_assay Name of the new normalized assay.
#' @param slope_window Inclusive lower and upper bounds on PC1 slope relative
#'   to the median; default `c(0.5, 1.5)`.
#' @export
normalize_stable_reference_pca <- function(se, assay = "quantity",
    output_assay = "quantity_pca_normalized", exclude = character(),
    candidate_count = 1000L, reference_count = 200L,
    slope_window = c(0.5, 1.5)) {
    if (!is.numeric(slope_window) || length(slope_window) != 2L ||
        any(!is.finite(slope_window)) || slope_window[[1L]] <= 0 ||
        slope_window[[2L]] < slope_window[[1L]]) {
        stop("slope_window must contain positive ordered bounds", call. = FALSE)
    }
    normalize_stable_reference(se, "pca", assay, output_assay, exclude,
        candidate_count, reference_count, slope_window)
}

normalize_stable_reference <- function(se, method, assay, output_assay, exclude,
    candidate_count, reference_count, slope_window = NULL) {
    if (!methods::is(se, "SummarizedExperiment"))
        stop("se must be a SummarizedExperiment", call. = FALSE)
    if (!is.character(assay) || length(assay) != 1L || is.na(assay) ||
        !assay %in% SummarizedExperiment::assayNames(se))
        stop("assay must name an existing assay", call. = FALSE)
    if (!is.character(output_assay) || length(output_assay) != 1L ||
        is.na(output_assay) || !nzchar(output_assay) ||
        output_assay %in% SummarizedExperiment::assayNames(se))
        stop("output_assay must name a new assay", call. = FALSE)
    if (!is.character(exclude) || anyNA(exclude))
        stop("exclude must contain feature IDs", call. = FALSE)
    valid_count <- function(value) {
        is.numeric(value) && length(value) == 1L &&
            is.finite(value) && value >= 1L && value <= .Machine$integer.max &&
            value == as.integer(value)
    }
    if (!valid_count(candidate_count) || !valid_count(reference_count) ||
        candidate_count < reference_count)
        stop("candidate_count and reference_count must be positive integers, ",
            "with candidate_count >= reference_count", call. = FALSE)
    values <- SummarizedExperiment::assay(se, assay)
    if (!is.matrix(values) || !is.numeric(values) || ncol(values) < 2L ||
        is.null(rownames(values)) || anyDuplicated(rownames(values)) ||
        any(values < 0 | is.infinite(values), na.rm = TRUE))
        stop("assay must have distinct feature IDs, at least two samples, ",
            "and nonnegative finite intensities", call. = FALSE)
    source_groups <- if ("source_group_id" %in%
        names(SummarizedExperiment::rowData(se))) {
        as.character(SummarizedExperiment::rowData(se)$source_group_id)
    } else {
        rep(NA_character_, nrow(se))
    }
    known <- c(rownames(values), source_groups)
    if (any(!exclude %in% known))
        stop("Unmatched exclude values: ",
            paste(exclude[!exclude %in% known], collapse = ", "), call. = FALSE)
    excluded <- rownames(values) %in% exclude | source_groups %in% exclude
    excluded[is.na(excluded)] <- FALSE
    positive <- values
    positive[!is.na(positive) & positive == 0] <- NA_real_
    logged <- log2(positive)
    complete <- rowSums(is.finite(logged)) == ncol(logged) &
        !excluded
    candidates <- logged[complete, , drop = FALSE]
    if (nrow(candidates) < reference_count)
        stop("Fewer complete, nonexcluded candidates than reference_count", call. = FALSE)
    ordered <- order(-rowMeans(candidates), rownames(candidates))
    candidates <- candidates[utils::head(ordered, candidate_count), , drop = FALSE]
    centered <- sweep(candidates, 1L, rowMeans(candidates), "-")
    diagnostics <- list()
    if (method == "genorm") {
        row_ss <- rowSums(centered^2)
        distances <- sqrt(pmax(outer(row_ss, row_ss, "+") -
            2 * tcrossprod(centered), 0) / (ncol(centered) - 1L))
        remaining <- seq_len(nrow(centered))
        while (length(remaining) > reference_count) {
            instability <- rowSums(distances[remaining, remaining, drop = FALSE]) /
                (length(remaining) - 1L)
            remaining <- remaining[-which.max(instability)]
        }
        selected <- rownames(centered)[remaining]
    } else {
        pc <- stats::prcomp(t(centered), center = FALSE, scale. = FALSE)
        scores <- pc$x[, 1L]
        if (sum(scores^2) <= .Machine$double.eps)
            stop("Candidates have no shared PCA variation", call. = FALSE)
        slopes <- as.vector(centered %*% scores / sum(scores^2))
        if (stats::median(slopes) < 0) {
            slopes <- -slopes
            scores <- -scores
        }
        slope_median <- stats::median(slopes)
        if (abs(slope_median) <= .Machine$double.eps)
            stop("Median PCA slope is zero", call. = FALSE)
        relative <- slopes / slope_median
        residuals <- centered - tcrossprod(slopes, scores)
        residual_sd <- sqrt(rowSums(residuals^2) / (ncol(centered) - 1L))
        eligible <- which(relative >= slope_window[[1L]] &
            relative <= slope_window[[2L]])
        if (length(eligible) < reference_count)
            stop("Fewer candidates pass the PCA slope window than reference_count", call. = FALSE)
        selected <- rownames(centered)[eligible[order(residual_sd[eligible],
            rownames(centered)[eligible])[seq_len(reference_count)]]]
        diagnostics <- list(pc1_variance_fraction = pc$sdev[[1L]]^2 /
            sum(pc$sdev^2), slope_window = slope_window)
    }
    reference_medians <- apply(logged[selected, , drop = FALSE], 2L, stats::median)
    offsets <- reference_medians - stats::median(reference_medians)
    normalized <- 2^sweep(logged, 2L, offsets, "-")
    if (any(is.infinite(normalized), na.rm = TRUE))
        stop("Normalization produced infinite intensities", call. = FALSE)
    prefix <- paste0("reference_", method)
    row_names <- paste0(prefix, c("_candidate", "_selected"))
    col_name <- paste0(prefix, "_log2_offset")
    if (any(row_names %in% names(SummarizedExperiment::rowData(se))) ||
        col_name %in% names(SummarizedExperiment::colData(se)))
        stop("Reference annotations already exist", call. = FALSE)
    row_data <- SummarizedExperiment::rowData(se)
    row_data[[row_names[[1L]]]] <- rownames(se) %in% rownames(candidates)
    row_data[[row_names[[2L]]]] <- rownames(se) %in% selected
    col_data <- SummarizedExperiment::colData(se)
    col_data[[col_name]] <- offsets
    SummarizedExperiment::rowData(se) <- row_data
    SummarizedExperiment::colData(se) <- col_data
    SummarizedExperiment::assay(se, output_assay) <- normalized
    info <- S4Vectors::metadata(se)
    info[[paste0("reference_normalization_", method)]] <- c(list(
        method = method, input_assay = assay, output_assay = output_assay,
        exclude = exclude, candidate_count = candidate_count,
        reference_count = reference_count, candidates = rownames(candidates),
        selected = selected), diagnostics)
    S4Vectors::metadata(se) <- info
    se
}
