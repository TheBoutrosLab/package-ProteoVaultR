make_reference_test_se <- function() {
    base <- seq(12, 19)
    shifts <- c(-1, -0.5, 0, 0.5, 1, 1.5)
    logged <- outer(base, shifts, "+")
    logged <- logged + outer(seq_along(base), c(0, 0.01, -0.01, 0, 0.01, -0.01),
        function(a, b) a * b / 10)
    logged <- rbind(logged, Target = c(15, 15, 15, 17, 17, 17))
    rownames(logged) <- c(paste0("R", seq_along(base)), "Target")
    colnames(logged) <- paste0("S", seq_len(ncol(logged)))
    SummarizedExperiment::SummarizedExperiment(
        assays = list(quantity = 2^logged),
        rowData = S4Vectors::DataFrame(source_group_id = c(
            paste0("sp|", seq_along(base), "|REF"), "sp|P48426|PI42A_HUMAN")),
        colData = S4Vectors::DataFrame(group = rep(c("A", "B"), each = 3)))
}

test_that("both stable-reference methods preserve input and coexist", {
    se <- make_reference_test_se()
    ge <- normalize_stable_reference_genorm(se,
        exclude = "sp|P48426|PI42A_HUMAN",
        candidate_count = 8, reference_count = 4)
    both <- normalize_stable_reference_pca(ge, exclude = "Target",
        candidate_count = 8, reference_count = 4)
    expect_equal(SummarizedExperiment::assay(both, "quantity"),
        SummarizedExperiment::assay(se, "quantity"))
    expect_true(all(c("quantity_genorm_normalized", "quantity_pca_normalized") %in%
        SummarizedExperiment::assayNames(both)))
    expect_false(SummarizedExperiment::rowData(both)$reference_genorm_candidate[[9L]])
    expect_false(SummarizedExperiment::rowData(both)$reference_pca_selected[[9L]])
    changed_labels <- se
    changed_labels$group <- rev(changed_labels$group)
    blind <- normalize_stable_reference_genorm(changed_labels, exclude = "Target",
        candidate_count = 8, reference_count = 4)
    expect_equal(SummarizedExperiment::assay(ge, "quantity_genorm_normalized"),
        SummarizedExperiment::assay(blind, "quantity_genorm_normalized"))
    for (method in c("genorm", "pca")) {
        info <- S4Vectors::metadata(both)[[paste0("reference_normalization_", method)]]
        expect_length(info$selected, 4L)
        offsets <- SummarizedExperiment::colData(both)[[paste0("reference_", method,
            "_log2_offset")]]
        expect_equal(stats::median(offsets), 0)
        expect_equal(SummarizedExperiment::assay(both,
            paste0("quantity_", method, "_normalized")),
            sweep(SummarizedExperiment::assay(se, "quantity"), 2L, 2^offsets, "/"))
    }
})

test_that("stable-reference methods reject insufficient candidates", {
    se <- make_reference_test_se()
    expect_error(normalize_stable_reference_genorm(se, exclude = "Target",
        reference_count = 9), "Fewer complete")
    expect_error(normalize_stable_reference_pca(se, candidate_count = 8,
        reference_count = 4, slope_window = c(2, 3)), "slope window")
    expect_error(normalize_stable_reference_genorm(se, exclude = "not_a_protein",
        candidate_count = 8, reference_count = 4), "Unmatched exclude")
})

test_that("reference diagnostics show before and after profiles for both methods", {
    se <- make_reference_test_se()
    se <- normalize_stable_reference_genorm(se, exclude = "Target",
        candidate_count = 8, reference_count = 4)
    se <- normalize_stable_reference_pca(se, exclude = "Target",
        candidate_count = 8, reference_count = 4)
    for (method in c("genorm", "pca")) {
        bpg <- plot_stable_reference_diagnostics(se, method)
        ggplot <- plot_stable_reference_diagnostics(se, method,
            style = "ggplot2")
        expect_s3_class(bpg, "trellis")
        expect_s3_class(ggplot, "ggplot")
        expect_equal(levels(ggplot$layers[[1L]]$data$stage),
            c("Before", "After"))
        expect_equal(sum(ggplot$layers[[2L]]$data$selected), 4 * ncol(se) * 2)
        expect_silent(ggplot2::ggplot_build(ggplot))
        grDevices::pdf(file = tempfile(fileext = ".pdf"))
        expect_silent(print(bpg))
        grDevices::dev.off()
    }
    expect_error(plot_stable_reference_diagnostics(make_reference_test_se()),
        "Run the selected")
})
