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
