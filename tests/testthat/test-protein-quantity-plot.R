test_that("protein-group quantity distributions use the normalized matrix", {
    path <- make_test_vault()
    edit_test_vault(
        path = path,
        statement = paste(
            "UPDATE protein_quant SET quantity=110",
            "WHERE group_id=91 AND channel_id=112 AND source_id=101"
        )
    )
    vault <- open_vault(path = path)
    on.exit(close_vault(vault = vault), add = TRUE)

    boxplot <- plot_protein_quantity_distribution(
        vault = vault,
        samples = c("S2", "S1")
    )
    violin <- plot_protein_quantity_distribution(
        vault = vault,
        plot_type = "violin"
    )

    expect_s3_class(boxplot, "ggplot")
    expect_s3_class(violin, "ggplot")
    expect_equal(nrow(boxplot$data), 3L)
    expect_setequal(as.character(unique(boxplot$data$group_id)), c("91", "92"))
    expect_equal(levels(boxplot$data$sample_name), c("S2", "S1"))
    expect_equal(boxplot$data$log2_quantity, log2(boxplot$data$quantity))
    raw_plot <- plot_protein_quantity_distribution(vault, transform = "none")
    expect_equal(nrow(raw_plot$data), 4L)
    expect_true(0 %in% raw_plot$data$quantity)
    expect_equal(raw_plot$data$plot_quantity, raw_plot$data$quantity)
    expect_equal(raw_plot$labels$y, "protein-group quantity")
    expect_equal(nrow(violin$data), 3L)
    expect_s3_class(violin$layers[[1L]]$geom, "GeomViolin")
    expect_silent(ggplot2::ggplot_build(boxplot))
    expect_error(
        plot_protein_quantity_distribution(
            vault = vault,
            plot_type = "scatter"
        ),
        "boxplot.*violin"
    )
})

test_that("protein-group distribution selects a report quantification method", {
    path <- make_test_vault()
    edit_test_vault(path, "INSERT INTO quant_methods VALUES (3,'PG.MaxLFQ')")
    edit_test_vault(path, paste(
        "INSERT INTO protein_groups VALUES",
        "(95,103,'P1;P2',1,'quantification',NULL,NULL)"
    ))
    edit_test_vault(path, "INSERT INTO protein_quant VALUES (95,111,103,3,500)")
    vault <- open_vault(path)
    on.exit(close_vault(vault), add = TRUE)

    plot <- plot_protein_quantity_distribution(vault, method = "PG.MaxLFQ")
    expect_equal(plot$data$quantity, 500)
    expect_equal(plot$data$quantity_type, factor("PG.MaxLFQ", levels = "PG.MaxLFQ"))
})
