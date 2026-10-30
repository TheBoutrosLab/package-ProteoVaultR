test_that("protein coverage marks OpenMS residue modifications", {
    path <- make_test_vault()
    edit_test_vault(
        path = path,
        statement = paste(
            "UPDATE peptidoforms SET",
            "modified_sequence='PEPT(Phospho)IDE', notation='OpenMS'",
            "WHERE peptidoform_id=41"
        )
    )
    vault <- open_vault(path = path)
    on.exit(close_vault(vault = vault), add = TRUE)

    plot <- plot_protein_coverage(
        vault = vault,
        accession = "P1"
    )
    expect_equal(length(plot$layers), 3L)
    markers <- plot$layers[[3L]]$data
    expect_equal(markers$position, 4L)
    expect_equal(markers$modification, "Phospho")
    expect_equal(markers$sample_name, "S1")
    expect_equal(as.character(markers$placement), "One position")
    expect_silent(ggplot2::ggplot_build(plot))
    expect_equal(plot$theme$legend.position, "right")
    expect_equal(plot$theme$legend.direction, "vertical")
    expect_equal(plot$theme$legend.box, "vertical")
    expect_equal(plot$scales$get_scales("shape")$guide, "none")

    filtered <- plot_protein_coverage(
        vault = vault,
        accession = "P1",
        ptm_names = "Oxidation"
    )
    expect_equal(length(filtered$layers), 2L)
    denylisted <- plot_protein_coverage(
        vault = vault,
        accession = "P1",
        ptm_denylist = "Phospho"
    )
    expect_equal(length(denylisted$layers), 2L)
    hidden <- plot_protein_coverage(
        vault = vault,
        accession = "P1",
        show_ptms = FALSE
    )
    expect_equal(length(hidden$layers), 2L)
})

test_that("repeated peptides retain ambiguous PTM positions", {
    path <- make_test_vault()
    edit_test_vault(
        path = path,
        statement = "UPDATE proteins SET sequence='PEPTIDEPEPTIDE' WHERE accession='P1'"
    )
    edit_test_vault(
        path = path,
        statement = paste(
            "UPDATE peptidoforms SET",
            "modified_sequence='PEPT(Phospho)IDE', notation='OpenMS'",
            "WHERE peptidoform_id=41"
        )
    )
    vault <- open_vault(path = path)
    on.exit(close_vault(vault = vault), add = TRUE)

    plot <- plot_protein_coverage(
        vault = vault,
        accession = "P1",
        samples = "S1"
    )
    markers <- plot$layers[[3L]]$data
    expect_equal(markers$position, c(4L, 11L))
    expect_true(all(as.character(markers$placement) == "Multiple positions"))
    expect_equal(
        plot$scales$get_scales("shape")$breaks,
        "Multiple positions"
    )
})

test_that("terminal OpenMS modifications are not assigned to residues", {
    path <- make_test_vault()
    edit_test_vault(
        path = path,
        statement = paste(
            "UPDATE peptidoforms SET",
            "modified_sequence='.(Acetyl)PEPTIDE', notation='OpenMS'",
            "WHERE peptidoform_id=41"
        )
    )
    vault <- open_vault(path = path)
    on.exit(close_vault(vault = vault), add = TRUE)

    expect_warning(
        plot <- plot_protein_coverage(vault = vault, accession = "P1"),
        "Terminal modifications are not shown"
    )
    expect_equal(length(plot$layers), 2L)
})

test_that("PTM denylist hides carbamidomethyl by default and is configurable", {
    path <- make_test_vault()
    edit_test_vault(
        path = path,
        statement = "UPDATE proteins SET sequence='PEPTCDE' WHERE accession='P1'"
    )
    edit_test_vault(
        path = path,
        statement = "UPDATE peptides SET sequence='PEPTCDE' WHERE peptide_id=31"
    )
    edit_test_vault(
        path = path,
        statement = paste(
            "UPDATE peptidoforms SET",
            "modified_sequence='PEPTC(Carbamidomethyl)DE', notation='OpenMS'",
            "WHERE peptidoform_id=41"
        )
    )
    vault <- open_vault(path = path)
    on.exit(close_vault(vault = vault), add = TRUE)

    default_plot <- plot_protein_coverage(vault = vault, accession = "P1")
    expect_equal(length(default_plot$layers), 2L)
    included <- plot_protein_coverage(
        vault = vault,
        accession = "P1",
        ptm_denylist = character()
    )
    expect_equal(included$layers[[3L]]$data$modification, "Carbamidomethyl")
    excluded <- plot_protein_coverage(
        vault = vault,
        accession = "P1",
        ptm_denylist = c("Carbamidomethyl", "Oxidation")
    )
    expect_equal(length(excluded$layers), 2L)
    expect_error(
        plot_protein_coverage(
            vault = vault,
            accession = "P1",
            ptm_denylist = NA_character_
        ),
        "ptm_denylist"
    )
})
