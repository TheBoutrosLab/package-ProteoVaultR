test_that("protein coverage wraps sequence and includes undetected samples", {
    path <- make_test_vault()
    vault <- open_vault(path = path)
    on.exit(close_vault(vault = vault), add = TRUE)

    plot <- plot_protein_coverage(
        vault = vault,
        accession = "P1",
        row_width = 4L
    )
    expect_s3_class(plot, "ggplot")
    built <- ggplot2::ggplot_build(plot)
    expect_equal(
        as.character(built$plot$scales$get_scales("fill")$get_breaks()),
        c("Not observed", "Observed")
    )
    rows <- plot$data
    expect_equal(levels(rows$block), c("1–4", "5–7"))
    expect_equal(
        as.character(rows$status[rows$sample_name == "S1"]),
        rep("Observed", 7)
    )
    expect_equal(
        as.character(rows$status[rows$sample_name == "S2"]),
        rep("Not observed", 7)
    )
    expect_equal(
        rows$residue[rows$sample_name == "Sequence"],
        strsplit("PEPTIDE", split = "", fixed = TRUE)[[1L]]
    )
    expect_match(plot$labels$subtitle, "7/7 residues")
})

test_that("repeated peptide positions are marked as possible", {
    path <- make_test_vault()
    edit_test_vault(
        path = path,
        statement = "UPDATE proteins SET sequence='PEPTIDEPEPTIDE' WHERE accession='P1'"
    )
    vault <- open_vault(path = path)
    on.exit(close_vault(vault = vault), add = TRUE)

    plot <- plot_protein_coverage(
        vault = vault,
        accession = "P1",
        samples = "S1",
        row_width = 7L
    )
    expect_true(all(
        as.character(plot$data$status[plot$data$sample_name == "S1"]) ==
            "Possible location"
    ))
    expect_match(plot$labels$subtitle, "0/14 residues")
    built <- ggplot2::ggplot_build(plot)
    expect_true(
        "Possible location" %in%
            as.character(built$plot$scales$get_scales("fill")$get_breaks())
    )
})

test_that("protein coverage requires a matching full sequence", {
    path <- make_test_vault()
    vault <- open_vault(path = path)
    on.exit(close_vault(vault = vault), add = TRUE)

    expect_error(
        plot_protein_coverage(vault = vault, accession = "P2"),
        "add the search FASTA"
    )
    expect_error(
        plot_protein_coverage(vault = vault, accession = "P1", samples = "missing"),
        "Unknown samples"
    )
})


test_that("overlapping peptide matches remain possible locations", {
    path <- make_test_vault()
    edit_test_vault(
        path = path,
        statement = "UPDATE proteins SET sequence='AAAAA' WHERE accession='P1'"
    )
    edit_test_vault(
        path = path,
        statement = "UPDATE peptides SET sequence='AAA' WHERE peptide_id=31"
    )
    vault <- open_vault(path = path)
    on.exit(close_vault(vault = vault), add = TRUE)

    plot <- plot_protein_coverage(
        vault = vault,
        accession = "P1",
        samples = "S1",
        row_width = 5L
    )
    expect_equal(
        as.character(plot$data$status[plot$data$sample_name == "S1"]),
        rep("Possible location", 5)
    )
})
