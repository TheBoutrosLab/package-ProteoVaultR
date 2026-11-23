make_volcano_test_fit <- function() {
    rows <- data.frame(
        feature_id = c("g1", "g2", "g3", "g4"),
        contrast = "Drug_vs_DMSO",
        logFC = c(1.5, -0.8, 2, NA_real_),
        P.Value = c(0.01, 0.05, 0, NA_real_),
        adj.P.Val = c(0.03, 0.1, 0, NA_real_),
        source_group_id = c("sp|O75083|WDR1_HUMAN", "P2", "P3", "P4")
    )
    rows$member_accessions <- I(list("O75083", "P2", "P3", "P4"))
    rows$leading_protein_accession <- c("O75083", "P2", "P3", "P4")
    rows$leading_gene_name <- c("WDR1", "GENE2", "GENE3", "GENE4")
    rows$member_gene_name <- rows$leading_gene_name
    structure(list(results = rows), class = "proteovault_limma")
}

test_that("volcano plots label transformed probabilities in both styles", {
    fit <- make_volcano_test_fit()
    ggplot <- plot_limma_volcano(fit, style = "ggplot2",
        proteins = c("O75083", "GENE2"))
    bpg <- plot_limma_volcano(fit, proteins = "WDR1")
    expect_s3_class(ggplot, "ggplot")
    expect_s3_class(bpg, "trellis")
    expect_equal(bpg$x.scales$tck, c(0.5, 0))
    expect_equal(bpg$y.scales$tck, c(0.5, 0))
    expect_equal(bpg$y.scales$labels[1:3], c("0", "1", "2"))
    key <- bpg$legend$right$args$key
    expect_equal(key$title, as.expression(bquote(bold("q-value"))))
    expect_equal(key$cex.title, bpg$xlab$cex)
    expect_equal(key$text$cex, bpg$x.scales$cex[[1L]])
    expect_equal(key$text$font, 2)
    expect_equal(ggplot$scales$get_scales("colour")$name, "q-value")
    expect_equal(ggplot$theme$legend.title$size, ggplot$theme$axis.title$size)
    expect_equal(ggplot$theme$legend.text$size, ggplot$theme$axis.text$size)
    expect_equal(ggplot$theme$legend.title$face, "bold")
    expect_equal(ggplot$theme$legend.text$face, "bold")
    expect_equal(levels(ggplot$data$significant), c("> 0.05",
        "≤ 0.05; FC < 0", "≤ 0.05; FC ≥ 0"))
    expect_equal(as.character(ggplot$data$significant[[1L]]),
        "≤ 0.05; FC ≥ 0")
    expect_equal(nrow(ggplot$data), 3L)
    expect_equal(ggplot$data$label, c("O75083", "GENE2", ""))
    expect_equal(ggplot$layers[[2L]]$data$feature_id, c("g1", "g2"))
    expect_equal(ggplot$layers[[3L]]$aes_params$size, 1.8)
    expect_null(ggplot$layers[[3L]]$aes_params$colour)
    expect_equal(ggplot2::ggplot_build(ggplot)$data[[3L]]$colour,
        c("#B24A59", "#777777"))
    labeled <- plot_limma_volcano(fit, style = "ggplot2",
        proteins = c("O75083", "P2"), label_by = "leading_gene_name")
    labeled_bpg <- plot_limma_volcano(fit,
        proteins = c("O75083", "P2"), label_by = "leading_gene_name")
    expect_equal(labeled$data$label, c("WDR1", "GENE2", ""))
    expect_s3_class(labeled_bpg, "trellis")
    expect_true(all(labeled$layers[[2L]]$data$label_y >
        labeled$layers[[2L]]$data$negative_log10))
    expect_error(plot_limma_volcano(fit, proteins = "O75083",
        label_by = "missing"), "label_by")
    expect_equal(ggplot$data$negative_log10[[1L]], -log10(0.03))
    expect_equal(ggplot$data$negative_log10[[3L]], -log10(0.003))
    expect_equal(ggplot$labels$y, "-log10(q-value)")
    expect_s3_class(ggplot$theme$panel.grid, "element_blank")
    expect_true(all(c("0", "1", "2") %in%
        ggplot2::ggplot_build(ggplot)$layout$panel_params[[1]]$y$get_labels()))
    expect_silent(ggplot2::ggplot_build(ggplot))
    grDevices::pdf(file = tempfile(fileext = ".pdf"))
    on.exit(grDevices::dev.off(), add = TRUE)
    expect_silent(print(bpg))
})

test_that("volcano selects p-values and validates contrasts and protein labels", {
    fit <- make_volcano_test_fit()
    plot <- plot_limma_volcano(fit, significance = "p_value", style = "ggplot2")
    expect_equal(plot$labels$y, "-log10(p-value)")
    expect_equal(plot$scales$get_scales("colour")$name, "p-value")
    expect_equal(as.character(plot$data$significant[[2L]]),
        "≤ 0.05; FC < 0")
    expect_equal(plot$data$negative_log10[[1L]], 2)
    expect_error(plot_limma_volcano(fit, proteins = "missing"),
        "No plotted feature matches")
    expect_equal(plot_limma_volcano(fit, proteins = "sp|O75083|WDR1_HUMAN",
        style = "ggplot2")$data$label[[1L]], "sp|O75083|WDR1_HUMAN")
    expect_error(plot_limma_volcano(fit, contrast = "missing"),
        "contrast must name")
    fit$results <- rbind(fit$results,
        transform(fit$results, contrast = "Other_vs_DMSO"))
    expect_error(plot_limma_volcano(fit), "Select one contrast")
    expect_s3_class(plot_limma_volcano(fit, contrast = "Other_vs_DMSO"),
        "trellis")
})

test_that("p-value axis can color by q-value cutoff in both styles", {
    fit <- make_volcano_test_fit()
    plot <- plot_limma_volcano(fit, significance = "p_value",
        color_by = "q_value", cutoff = 0.04, style = "ggplot2")
    bpg <- plot_limma_volcano(fit, significance = "p_value",
        color_by = "q_value", cutoff = 0.04)
    expect_equal(plot$labels$y, "-log10(p-value)")
    expect_equal(plot$scales$get_scales("colour")$name, "q-value")
    expect_equal(plot$data$negative_log10[[2L]], -log10(0.05))
    expect_equal(as.character(plot$data$significant),
        c("≤ 0.04; FC ≥ 0", "> 0.04", "≤ 0.04; FC ≥ 0"))
    labeled <- plot_limma_volcano(fit, significance = "p_value",
        color_by = "q_value", cutoff = 0.04, style = "ggplot2",
        proteins = c("O75083", "P2"))
    expect_equal(ggplot2::ggplot_build(labeled)$data[[3L]]$colour,
        c("#B24A59", "#777777"))
    expect_equal(bpg$legend$right$args$key$title,
        as.expression(bquote(bold("q-value"))))
    expect_error(plot_limma_volcano(fit, significance = "q_value",
        color_by = "p_value"), "color_by must be q_value")
})

test_that("ggplot legend draws both cutoff symbols when no feature passes", {
    fit <- make_volcano_test_fit()
    fit$results$adj.P.Val[1:3] <- 0.2
    plot <- plot_limma_volcano(fit, style = "ggplot2")
    grDevices::pdf(file = tempfile(fileext = ".pdf"))
    on.exit(grDevices::dev.off(), add = TRUE)
    grob <- ggplot2::ggplotGrob(plot)
    guide <- grob$grobs[[which(grob$layout$name == "guide-box-right")]]$grobs[[1L]]
    keys <- guide$grobs[grep("^key-", guide$layout$name)]
    colours <- vapply(keys, function(key) key$children[[2L]]$gp$col, character(1))
    expect_equal(colours, c("#777777B2", "#4C78A8B2", "#B24A59B2"))
})
