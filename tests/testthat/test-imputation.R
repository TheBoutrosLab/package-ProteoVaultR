test_that("conditional LOD fills only wholly missing bins", {
    values <- rbind(P1 = c(10, 12, NA, NA),
        P2 = c(8, NA, 9, 11), P3 = c(NA, NA, 4, NA))
    colnames(values) <- paste0("S", 1:4)
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(quantity = values),
        colData = S4Vectors::DataFrame(group = c("A", "A", "B", "B")))
    output <- impute_conditional_lod(se, "group", min_observed = 2)
    expect_equal(SummarizedExperiment::assay(output, "quantity"), values)
    expect_equal(SummarizedExperiment::assay(output, "quantity_imputed")["P1", ],
        c(S1 = 10, S2 = 12, S3 = 2, S4 = 5.5))
    expect_true(SummarizedExperiment::rowData(output)$imputed_any[[1L]])
    expect_equal(SummarizedExperiment::rowData(output)$imputation_status_B[[1L]],
        "lod_imputed")
    expect_true(all(is.na(SummarizedExperiment::assay(output, "quantity_imputed")["P3", ])))
})
