test_that("protein precursor plot maps groups and quantification method", {
    path <- make_test_vault()
    edit_test_vault(path, "INSERT INTO precursor_groups VALUES (51,91),(52,92)")
    edit_test_vault(path,
        "UPDATE precursor_quant SET quantity=15 WHERE precursor_id=51 AND channel_id=112 AND source_id=102")
    vault <- open_vault(path)
    on.exit(close_vault(vault), add = TRUE)
    vault <- set_sample_metadata(vault, data.frame(
        sample_name = c("S1", "S2"), sample_id = c("Control", "Drug")))
    proteins <- c(First = "P1;P2", Third = "P3")
    ggplot <- plot_protein_precursor_quantity_lineplot(vault, proteins, method = "matrix",
        style = "ggplot2")
    bpg <- plot_protein_precursor_quantity_lineplot(vault, proteins, method = "matrix")
    expect_s3_class(ggplot, "ggplot")
    expect_s3_class(bpg, "trellis")
    expect_setequal(as.character(unique(ggplot$data$protein)), c("First", "Third"))
    expect_equal(unique(ggplot$data$precursor[ggplot$data$protein == "First"]),
        "PEPTIDE z2")
    expect_equal(levels(ggplot$data$sample), c("Control", "Drug"))
    expect_equal(ggplot$data$log2_quantity[ggplot$data$sample_name == "S2"],
        log2(c(15, 20)))
    expect_silent(ggplot2::ggplot_build(ggplot))
    grDevices::pdf(file = tempfile(fileext = ".pdf"))
    on.exit(grDevices::dev.off(), add = TRUE)
    expect_silent(print(bpg))
    custom <- plot_protein_precursor_quantity_lineplot(vault, proteins,
        method = "matrix", style = "ggplot2",
        colours = c("OTHER z3" = "#123456", "PEPTIDE z2" = "#ABCDEF"))
    expect_equal(custom$scales$get_scales("colour")$palette(2),
        c("OTHER z3" = "#123456", "PEPTIDE z2" = "#ABCDEF"))
    expect_error(plot_protein_precursor_quantity_lineplot(vault, proteins,
        method = "matrix", colours = "red"), "one value")
    report <- plot_protein_precursor_quantity_lineplot(vault, c(First = "P1;P2"),
        method = "Precursor.Quantity", style = "ggplot2")
    expect_equal(report$data$quantity, 5)
    expect_error(plot_protein_precursor_quantity_lineplot(vault, proteins,
        method = "Precursor.Quantity"), "No positive precursor quantities")
})
