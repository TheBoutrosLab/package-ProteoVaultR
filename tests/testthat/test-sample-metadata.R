test_that("sample metadata labels plots without changing vault identifiers", {
    path <- make_test_vault()
    original <- open_vault(path)
    on.exit(close_vault(original), add = TRUE)
    vault <- set_sample_metadata(
        original,
        data.frame(
            sample_name = c("S1", "S2"),
            sample_id = c("Control", "Treatment"),
            batch = c("A", "B")
        )
    )

    expect_null(sample_metadata(original))
    expect_equal(sample_metadata(vault)$batch, c("A", "B"))
    expect_equal(DBI::dbGetQuery(vault_connection(vault),
        "SELECT sample_name FROM samples ORDER BY sample_name")$sample_name,
        c("S1", "S2", "S3"))

    protein <- plot_protein_quantity_distribution(vault, samples = c("S2", "S1"))
    precursor <- plot_precursor_quantity_distribution(vault)
    expect_equal(levels(protein$data$sample_name), c("S2", "S1"))
    expect_equal(as.character(ggplot2::ggplot_build(protein)$layout$panel_params[[1]]$x$get_labels()),
        c("Treatment", "Control"))
    expect_equal(as.character(ggplot2::ggplot_build(precursor)$layout$panel_params[[1]]$x$get_labels()),
        c("Control", "Treatment"))
    expect_equal(colnames(matrix_protein_abundance(vault)), c("S1", "S2"))

    reopened <- open_vault(path)
    on.exit(close_vault(reopened), add = TRUE)
    expect_null(sample_metadata(reopened))
})

test_that("sample metadata validates join keys and labels", {
    vault <- open_vault(make_test_vault())
    on.exit(close_vault(vault), add = TRUE)
    expect_error(set_sample_metadata(vault, data.frame(sample_id = "short")),
        "sample_name column")
    expect_error(set_sample_metadata(vault, data.frame(sample_name = c("S1", "S1"))),
        "unique")
    expect_error(set_sample_metadata(vault, data.frame(sample_name = "missing")),
        "absent")
    expect_error(set_sample_metadata(vault,
        data.frame(sample_name = "S1", sample_id = "S2")), "conflicts")
    expect_error(set_sample_metadata(vault,
        data.frame(sample_name = c("S1", "S2"), sample_id = c("same", "same"))),
        "unique")
})
