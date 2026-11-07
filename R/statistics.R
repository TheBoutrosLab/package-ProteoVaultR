#' Fit limma contrasts to ProteoVault quantities
#'
#' Fits a linear model to one assay of a `SummarizedExperiment`, usually from
#' [as_summarized_experiment()]. The design formula uses columns of `colData`.
#' Design column names are converted with [base::make.names()] so treatment
#' labels containing hyphens can be used in contrast expressions. The returned
#' `design` matrix shows the exact names to use.
#'
#' The default expects positive intensity values and applies log2 without a
#' pseudocount. Zeros become missing values before transformation; negative or
#' infinite intensities are rejected. Missing values are passed to limma without
#' imputation. If the assay is already log2 transformed, set
#' `log2_transform = FALSE`; zero values are then retained by default.
#' Normalization and feature filtering remain the caller's decisions.
#'
#' @param se A [SummarizedExperiment::SummarizedExperiment] with a numeric
#'   feature-by-sample assay and sample information in `colData`.
#' @param design A model formula, such as `~ 0 + condition + batch`, or a numeric
#'   design matrix with one row per sample in assay column order. The model must
#'   have full column rank and residual degrees of freedom.
#' @param contrasts A named character vector of limma contrast expressions,
#'   such as `c(Drug_vs_DMSO = "conditionDrug - conditionDMSO")`.
#'   Expressions refer to the columns of the returned design matrix.
#' @param assay Name of the assay to analyze; defaults to `"quantity"`.
#' @param log2_transform Whether to apply log2 to positive intensities.
#' @param zero_as_missing Whether to replace zero assay values with `NA` before
#'   fitting. Defaults to `log2_transform`, so zeros of an already log2 assay
#'   are retained.
#' @param trend Passed to [limma::eBayes()]; use `TRUE` to model a mean-variance
#'   trend if justified by diagnostic plots.
#' @param robust Passed to [limma::eBayes()] for robust empirical Bayes
#'   variance moderation.
#' @return A list of class `proteovault_limma` with `fit` (the moderated limma
#'   fit), `results` (one row per feature and contrast, with BH-adjusted
#'   p-values and `rowData` annotations), `design`, `contrast_matrix`, and
#'   `settings`. Results are not filtered by p-value or fold change.
#' @details For a group-means design, write `~ 0 + condition` and a contrast
#'   such as `"conditionTreatment - conditionDMSO"`. A positive `logFC` means
#'   greater abundance in Treatment. If a condition is named `SBI-9456`, its
#'   design column becomes `conditionSBI.9456`. Each selected contrast is fitted
#'   against the same model; other treatment groups may remain in the design.
#'   The assay should contain comparable quantities from one source and method.
#'   Count-data workflows requiring `voom()` are outside this function.
#' @examples
#' \dontrun{
#' se <- as_summarized_experiment(vault, level = "protein_group")
#' fit <- fit_limma(
#'   se,
#'   design = ~ 0 + condition,
#'   contrasts = c(SBI_vs_DMSO = "conditionSBI.9456 - conditionDMSO")
#' )
#' head(fit$results)
#' colnames(fit$design)
#' }
#' @export
fit_limma <- function(
    se,
    design,
    contrasts,
    assay = "quantity",
    log2_transform = TRUE,
    zero_as_missing = log2_transform,
    trend = FALSE,
    robust = FALSE
) {
    if (!methods::is(se, "SummarizedExperiment")) {
        stop("se must be a SummarizedExperiment", call. = FALSE)
    }
    if (!is.character(assay) || length(assay) != 1L ||
        is.na(assay) || !assay %in% SummarizedExperiment::assayNames(se)) {
        stop("assay must name an assay in se", call. = FALSE)
    }
    for (name in c("log2_transform", "zero_as_missing", "trend", "robust")) {
        value <- get(name)
        if (!is.logical(value) || length(value) != 1L || is.na(value)) {
            stop(name, " must be TRUE or FALSE", call. = FALSE)
        }
    }
    values <- SummarizedExperiment::assay(se, assay)
    if (!is.matrix(values) || !is.numeric(values) ||
        !nrow(values) || !ncol(values)) {
        stop("assay must be a nonempty numeric matrix", call. = FALSE)
    }
    if (any(is.infinite(values))) {
        stop("assay contains infinite values", call. = FALSE)
    }
    if (log2_transform && any(values < 0, na.rm = TRUE)) {
        stop("negative intensities cannot be log2 transformed", call. = FALSE)
    }
    if (zero_as_missing) values[!is.na(values) & values == 0] <- NA_real_
    if (log2_transform) values <- log2(values)

    design_matrix <- build_limma_design(se, design)
    if (!is.character(contrasts) || !length(contrasts) ||
        is.null(names(contrasts)) || anyNA(contrasts) ||
        any(!nzchar(contrasts)) || anyNA(names(contrasts)) ||
        any(!nzchar(names(contrasts))) || anyDuplicated(names(contrasts))) {
        stop("contrasts must be a named, nonempty character vector", call. = FALSE)
    }
    contrast_matrix <- tryCatch(
        limma::makeContrasts(contrasts = unname(contrasts), levels = design_matrix),
        error = function(error) {
            stop("Cannot evaluate contrasts against design columns ",
                paste(colnames(design_matrix), collapse = ", "), ": ",
                conditionMessage(error), call. = FALSE)
        }
    )
    colnames(contrast_matrix) <- names(contrasts)
    fitted <- limma::lmFit(values, design_matrix)
    fitted <- limma::contrasts.fit(fitted, contrast_matrix)
    fitted <- limma::eBayes(fitted, trend = trend, robust = robust)

    annotation <- as.data.frame(SummarizedExperiment::rowData(se))
    annotation <- annotation[, setdiff(names(annotation),
        c("feature_id", "contrast", "logFC", "AveExpr", "t", "P.Value",
            "adj.P.Val", "B")), drop = FALSE]
    output <- lapply(seq_along(contrasts), function(index) {
        table <- limma::topTable(fitted, coef = index, number = Inf,
            adjust.method = "BH", sort.by = "none")
        indices <- match(rownames(table), rownames(se))
        cbind(
            data.frame(feature_id = rownames(table),
                contrast = names(contrasts)[[index]], stringsAsFactors = FALSE),
            table[, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B"),
                drop = FALSE],
            annotation[indices, , drop = FALSE]
        )
    })
    result <- list(
        fit = fitted,
        results = tibble::as_tibble(do.call(rbind, output)),
        design = design_matrix,
        contrast_matrix = contrast_matrix,
        settings = list(assay = assay, log2_transform = log2_transform,
            zero_as_missing = zero_as_missing, trend = trend, robust = robust)
    )
    class(result) <- "proteovault_limma"
    result
}

build_limma_design <- function(se, design) {
    sample_data <- as.data.frame(SummarizedExperiment::colData(se))
    if (inherits(design, "formula")) {
        design_matrix <- tryCatch({
            frame <- stats::model.frame(design, data = sample_data,
                na.action = stats::na.fail)
            stats::model.matrix(design, data = frame)
        }, error = function(error) {
            stop("Cannot build design from colData: ", conditionMessage(error),
                call. = FALSE)
        })
    } else if (is.matrix(design) && is.numeric(design)) {
        design_matrix <- design
    } else {
        stop("design must be a formula or numeric matrix", call. = FALSE)
    }
    if (nrow(design_matrix) != ncol(se) || !ncol(design_matrix) ||
        anyNA(design_matrix) || any(!is.finite(design_matrix))) {
        stop("design must have one finite, complete row per sample", call. = FALSE)
    }
    if (!is.null(rownames(design_matrix)) &&
        !identical(rownames(design_matrix), colnames(se))) {
        stop("design row names must match assay sample names in order", call. = FALSE)
    }
    if (is.null(colnames(design_matrix)) || any(!nzchar(colnames(design_matrix)))) {
        stop("design matrix needs nonempty column names", call. = FALSE)
    }
    colnames(design_matrix) <- make.names(colnames(design_matrix), unique = TRUE)
    rownames(design_matrix) <- colnames(se)
    if (qr(design_matrix)$rank != ncol(design_matrix)) {
        stop("design matrix is not full rank", call. = FALSE)
    }
    if (nrow(design_matrix) <= ncol(design_matrix)) {
        stop("design needs residual degrees of freedom", call. = FALSE)
    }
    design_matrix
}

#' List pairwise limma contrasts for a sample group
#'
#' Builds all pairwise comparisons between levels of a categorical `colData`
#' variable in a design formula. Each expression is ready for the `contrasts`
#' argument of [fit_limma()]. The direction is later level minus earlier level.
#' Other additive design terms are held constant. Interactions are not supported.
#'
#' @param se A [SummarizedExperiment::SummarizedExperiment].
#' @param design A model formula or a single string containing a formula, such
#'   as `"~ 0 + condition + batch"`.
#' @param group Name of the categorical `colData` variable to compare. May be
#'   omitted when the design contains exactly one categorical variable.
#' @return A named character vector of limma contrast expressions.
#' @examples
#' \dontrun{
#' list_limma_contrasts(se, "~ 0 + condition")
#' }
#' @export
list_limma_contrasts <- function(se, design, group = NULL) {
    if (!methods::is(se, "SummarizedExperiment")) {
        stop("se must be a SummarizedExperiment", call. = FALSE)
    }
    if (is.character(design) && length(design) == 1L && !is.na(design)) {
        design <- tryCatch(stats::as.formula(design), error = function(error) {
            stop("design must be a valid formula", call. = FALSE)
        })
    }
    if (!inherits(design, "formula")) {
        stop("design must be a formula or one formula string", call. = FALSE)
    }
    terms <- stats::terms(design)
    if (any(attr(terms, "order") > 1L)) {
        stop("Interactions are not supported by list_limma_contrasts", call. = FALSE)
    }
    sample_data <- as.data.frame(SummarizedExperiment::colData(se))
    design_matrix <- build_limma_design(se, design)
    variables <- all.vars(design)
    candidates <- variables[vapply(sample_data[variables], function(value) {
        is.factor(value) || is.character(value)
    }, logical(1))]
    if (is.null(group)) {
        if (length(candidates) != 1L) {
            stop("Specify group; categorical design variables are: ",
                paste(candidates, collapse = ", "), call. = FALSE)
        }
        group <- candidates
    }
    if (!is.character(group) || length(group) != 1L || is.na(group) ||
        !group %in% candidates) {
        stop("group must name a categorical variable in the design", call. = FALSE)
    }
    frame <- stats::model.frame(design, data = sample_data, na.action = stats::na.fail)
    for (column in names(frame)) {
        if (is.character(frame[[column]])) frame[[column]] <- factor(frame[[column]])
    }
    levels <- levels(as.factor(frame[[group]]))
    if (length(levels) < 2L) {
        stop("group needs at least two observed levels", call. = FALSE)
    }
    profiles <- lapply(levels, function(level) {
        row <- frame[1L, , drop = FALSE]
        row[[group]] <- factor(level, levels = levels)
        profile <- stats::model.matrix(design, data = row)
        colnames(profile) <- make.names(colnames(profile), unique = TRUE)
        as.numeric(profile[1L, colnames(design_matrix)])
    })
    pairs <- utils::combn(seq_along(levels), 2L)
    expressions <- apply(pairs, 2L, function(indices) {
        difference <- profiles[[indices[[2L]]]] - profiles[[indices[[1L]]]]
        positive <- colnames(design_matrix)[difference == 1]
        negative <- colnames(design_matrix)[difference == -1]
        if (any(!difference %in% c(-1, 0, 1))) {
            stop("Cannot express group comparison using simple design coefficients",
                call. = FALSE)
        }
        negative_terms <- if (length(negative)) paste0("- ", negative) else character()
        expression <- paste(c(positive, negative_terms), collapse = " + ")
        expression <- sub("\\+ - ", "- ", expression)
        expression
    })
    names(expressions) <- make.unique(vapply(seq_len(ncol(pairs)), function(index) {
        paste0(make.names(levels[[pairs[2L, index]]]), "_vs_",
            make.names(levels[[pairs[1L, index]]]))
    }, character(1)))
    expressions
}
