test_that("protein abundance plot selects one group and uses sample labels", {
    vault <- open_vault(make_test_vault())
    on.exit(close_vault(vault), add = TRUE)
    vault <- set_sample_metadata(vault, data.frame(
        sample_name = c("S1", "S2"), sample_id = c("Control", "Treatment"),
        group = c("DMSO", "SBI-9456")
    ))
    plot <- plot_protein_abundance(vault, "P1;P2", samples = c("S2", "S1"))
    ggplot <- plot_protein_abundance(vault, "P1;P2",
        samples = c("S2", "S1"), style = "ggplot2")
    expect_s3_class(plot, "trellis")
    expect_s3_class(ggplot, "ggplot")
    expect_equal(ggplot$layers[[1L]]$aes_params$colour, "black")
    expect_equal(ggplot$layers[[1L]]$aes_params$linewidth, 0.7)
    expect_equal(levels(ggplot$data$sample_name), c("S2", "S1"))
    expect_equal(ggplot$data$quantity, c(100, NA_real_))
    expect_equal(
        levels(ggplot$data$sample), c("Control", "Treatment")
    )
    grDevices::pdf(file = tempfile(fileext = ".pdf"))
    on.exit(grDevices::dev.off(), add = TRUE)
    expect_silent(print(plot))
    colored <- plot_protein_abundance(vault, "P1;P2", color_by = "group")
    colored_ggplot <- plot_protein_abundance(vault, "P1;P2",
        color_by = "group", style = "ggplot2")
    expect_s3_class(colored, "trellis")
    expect_s3_class(colored_ggplot, "ggplot")
    expect_equal(as.character(colored_ggplot$data$group), c("DMSO", "SBI-9456"))
    expect_silent(print(colored))
    expect_error(plot_protein_abundance(vault, "P1;P2", color_by = "missing"),
        "color_by")
    expect_error(plot_protein_abundance(vault, "missing"), "No quantities")
    expect_error(plot_protein_abundance(vault, "P1;P2", samples = "unknown"), "samples")
})

test_that("BPG colors quantified bars by sample metadata", {
    path <- make_test_vault()
    edit_test_vault(path,
        "UPDATE protein_quant SET quantity=120 WHERE group_id=91 AND channel_id=112")
    vault <- open_vault(path)
    on.exit(close_vault(vault), add = TRUE)
    vault <- set_sample_metadata(vault, data.frame(
        sample_name = c("S1", "S2"),
        group = c("DMSO", "SBI-9456")
    ))
    plot <- plot_protein_abundance(vault, "P1;P2", color_by = "group")
    expect_s3_class(plot, "trellis")
    grDevices::pdf(file = tempfile(fileext = ".pdf"))
    on.exit(grDevices::dev.off(), add = TRUE)
    expect_silent(print(plot))
})
