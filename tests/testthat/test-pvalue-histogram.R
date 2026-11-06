make_histogram_test_fit <- function() {
    structure(list(results = data.frame(
        contrast = "Drug_vs_DMSO",
        P.Value = c(0.01, 0.05, 0, NA_real_)
    )), class = "proteovault_limma")
}

test_that("p-value histogram uses raw p-values in both styles", {
    fit <- make_histogram_test_fit()
    bpg <- plot_limma_pvalue_histogram(fit, bins = 10)
    ggplot <- plot_limma_pvalue_histogram(fit, bins = 10, style = "ggplot2")
    expect_s3_class(bpg, "trellis")
    expect_s3_class(ggplot, "ggplot")
    expect_equal(ggplot$data$p_value, c(0.01, 0.05, 0))
    expect_equal(ggplot$labels$x, "p-value")
    expect_equal(ggplot$labels$y, "Count")
    expect_equal(bpg$x.scales$tck, c(0.5, 0))
    expect_equal(bpg$y.scales$tck, c(0.5, 0))
    expect_equal(sum(ggplot2::ggplot_build(ggplot)$data[[1L]]$count), 3)
    grDevices::pdf(file = tempfile(fileext = ".pdf"))
    on.exit(grDevices::dev.off(), add = TRUE)
    expect_silent(print(bpg))
})

test_that("p-value histogram validates its contrast and values", {
    fit <- make_histogram_test_fit()
    expect_error(plot_limma_pvalue_histogram(fit, bins = 0), "bins")
    expect_error(plot_limma_pvalue_histogram(fit, bins = 2.5), "bins")
    expect_error(plot_limma_pvalue_histogram(fit, contrast = "missing"), "contrast")
    fit$results$P.Value[1:3] <- NA_real_
    expect_error(plot_limma_pvalue_histogram(fit), "No finite p-values")
    fit$results$P.Value[[1L]] <- 1.2
    expect_error(plot_limma_pvalue_histogram(fit), "between 0 and 1")
})
