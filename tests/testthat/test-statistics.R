make_limma_test_se <- function() {
    values <- rbind(
        P1 = c(100, 110, 90, 200, 210, 190, 120, 125, 115),
        P2 = c(80, 85, 75, 82, 88, 78, 150, 160, 140),
        P3 = c(60, 65, 55, 0, 63, 57, 61, 66, 56)
    )
    colnames(values) <- paste0("S", seq_len(ncol(values)))
    SummarizedExperiment::SummarizedExperiment(
        assays = list(quantity = values),
        rowData = S4Vectors::DataFrame(gene = c("A", "B", "C")),
        colData = S4Vectors::DataFrame(
            condition = factor(rep(c("DMSO", "SBI-9456", "Other"), each = 3),
                levels = c("DMSO", "SBI-9456", "Other"))
        )
    )
}

test_that("limma tests selected treatments against DMSO and retains annotations", {
    se <- make_limma_test_se()
    result <- fit_limma(se, design = ~ 0 + condition,
        contrasts = c(SBI_vs_DMSO = "conditionSBI.9456 - conditionDMSO"))
    expect_s3_class(result, "proteovault_limma")
    expect_equal(colnames(result$design),
        c("conditionDMSO", "conditionSBI.9456", "conditionOther"))
    expect_equal(colnames(result$contrast_matrix), "SBI_vs_DMSO")
    expect_equal(nrow(result$results), 3L)
    expect_equal(result$results$gene, c("A", "B", "C"))
    expect_true(result$results$logFC[result$results$feature_id == "P1"] > 0.8)
    expect_true(all(result$results$adj.P.Val >= 0 &
        result$results$adj.P.Val <= 1, na.rm = TRUE))
    expect_equal(result$settings$zero_as_missing, TRUE)
    expect_equal(SummarizedExperiment::assay(se)["P3", "S4"], 0)
})

test_that("limma accepts several contrasts and an already log2 assay", {
    se <- make_limma_test_se()
    values <- log2(SummarizedExperiment::assay(se) + 1)
    SummarizedExperiment::assay(se) <- values
    result <- fit_limma(se, design = ~ 0 + condition,
        contrasts = c(
            SBI_vs_DMSO = "conditionSBI.9456 - conditionDMSO",
            Other_vs_DMSO = "conditionOther - conditionDMSO"
        ), log2_transform = FALSE)
    expect_equal(nrow(result$results), 6L)
    expect_setequal(result$results$contrast,
        c("SBI_vs_DMSO", "Other_vs_DMSO"))
    expect_false(result$settings$zero_as_missing)
})

test_that("limma rejects ambiguous designs, invalid contrasts, and intensities", {
    se <- make_limma_test_se()
    expect_error(fit_limma(se, ~ 0 + condition,
        c(test = "missing - conditionDMSO")), "design columns")
    expect_error(fit_limma(se, ~ 0 + condition,
        "conditionSBI.9456 - conditionDMSO"), "named")
    values <- SummarizedExperiment::assay(se)
    values[1, 1] <- -1
    SummarizedExperiment::assay(se) <- values
    expect_error(fit_limma(se, ~ 0 + condition,
        c(test = "conditionSBI.9456 - conditionDMSO")), "negative")
    values[1, 1] <- 100
    SummarizedExperiment::assay(se) <- values
    matrix_design <- stats::model.matrix(~ 0 + condition,
        data = as.data.frame(SummarizedExperiment::colData(se)))
    rownames(matrix_design) <- rev(colnames(se))
    expect_error(fit_limma(se, matrix_design,
        c(test = "conditionSBI.9456 - conditionDMSO")), "row names")
    rownames(matrix_design) <- colnames(se)
    colnames(matrix_design) <- c("DMSO", "SBI", "Other")
    matrix_design[, 3] <- matrix_design[, 2]
    expect_error(fit_limma(se, matrix_design,
        c(test = "SBI - DMSO")), "full rank")
    se$condition[1] <- NA
    expect_error(fit_limma(se, ~ 0 + condition,
        c(test = "conditionSBI.9456 - conditionDMSO")), "Cannot build design")
})

test_that("limma reports an all-missing feature without inventing an effect", {
    se <- make_limma_test_se()
    values <- SummarizedExperiment::assay(se)
    values[3, ] <- 0
    SummarizedExperiment::assay(se) <- values
    result <- fit_limma(se, ~ 0 + condition,
        c(test = "conditionSBI.9456 - conditionDMSO"))
    expect_true(is.na(result$results$logFC[result$results$feature_id == "P3"]))
})

test_that("available contrasts match no-intercept and intercept designs", {
    se <- make_limma_test_se()
    contrasts <- list_limma_contrasts(se, "~ 0 + condition")
    expect_equal(contrasts, c(
        SBI.9456_vs_DMSO = "conditionSBI.9456 - conditionDMSO",
        Other_vs_DMSO = "conditionOther - conditionDMSO",
        Other_vs_SBI.9456 = "conditionOther - conditionSBI.9456"
    ))
    expect_equal(list_limma_contrasts(se, ~ condition), c(
        SBI.9456_vs_DMSO = "conditionSBI.9456",
        Other_vs_DMSO = "conditionOther",
        Other_vs_SBI.9456 = "conditionOther - conditionSBI.9456"
    ))
    fit <- fit_limma(se, ~ 0 + condition, contrasts)
    expect_equal(colnames(fit$contrast_matrix), names(contrasts))
})

test_that("available contrasts select a grouping variable in additive designs", {
    se <- make_limma_test_se()
    se$batch <- factor(rep(c("A", "B", "C"), 3))
    expect_error(list_limma_contrasts(se, ~ 0 + condition + batch),
        "Specify group")
    contrasts <- list_limma_contrasts(se, ~ 0 + condition + batch,
        group = "condition")
    expect_equal(contrasts[[1L]], "conditionSBI.9456 - conditionDMSO")
    expect_error(list_limma_contrasts(se, ~ condition * batch,
        group = "condition"), "Interactions")
})
