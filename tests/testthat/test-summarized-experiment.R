test_that("protein-group export keeps group quantities and one-to-many annotations", {
    path <- make_test_vault()
    edit_test_vault(path,
        "UPDATE protein_groups SET source_group_id = 'P2;P1' WHERE group_id = 91")
    edit_test_vault(path, paste0(
        "INSERT INTO source_files VALUES ",
        "(109,'/sources/mapping','mapping','uniprot_idmapping','annotation',",
        "10,1720000000000000001,'imported',NULL)"
    ))
    edit_test_vault(path, paste0(
        "INSERT INTO protein_identifiers VALUES ",
        "(61,109,'P1','Gene_Name','GENE1'),",
        "(61,109,'P1','Gene_Name','GENE1_ALT'),",
        "(62,109,'P2','Gene_Name','GENE2')"
    ))
    vault <- open_vault(path)
    on.exit(close_vault(vault), add = TRUE)
    vault <- set_sample_metadata(vault, data.frame(
        sample_name = c("S1", "S2"), sample_id = c("Control", "Treated"),
        condition = c("DMSO", "Drug")
    ))
    result <- as_summarized_experiment(vault, samples = c("S2", "S1"),
        mapping_types = "Gene_Name")
    expect_s4_class(result, "SummarizedExperiment")
    expect_equal(dim(result), c(2L, 2L))
    expect_equal(colnames(result), c("S2", "S1"))
    expect_equal(SummarizedExperiment::assay(result, "quantity")["91", "S1"], 100)
    expect_true(is.na(SummarizedExperiment::assay(result, "quantity")["91", "S2"]))
    rows <- SummarizedExperiment::rowData(result)
    expect_equal(sort(rows$member_accessions[[1L]]), c("P1", "P2"))
    expect_false("mapped_identifiers" %in% names(rows))
    expect_equal(rows$leading_gene_name[[1L]], "GENE2")
    expect_equal(rows$member_gene_name[[1L]], "GENE1;GENE1_ALT;GENE2")
    expect_equal(rows$leading_protein_id[[1L]], "62")
    expect_equal(rows$leading_protein_accession[[1L]], "P2")
    columns <- SummarizedExperiment::colData(result)
    expect_equal(as.character(columns$sample_id), c("Treated", "Control"))
    expect_equal(as.character(columns$condition), c("Drug", "DMSO"))
    expect_equal(columns$run_names[[2L]], "S1")
})

test_that("precursor and peptide exports preserve their distinct feature levels", {
    path <- make_test_vault()
    edit_test_vault(path,
        "INSERT INTO protein_groups VALUES (96,102,'P1;P2',1,'quantification',NULL,NULL)")
    edit_test_vault(path,
        "INSERT INTO protein_group_members VALUES (96,61),(96,62)")
    edit_test_vault(path, "INSERT INTO precursor_groups VALUES (51,96)")
    edit_test_vault(path, paste0(
        "INSERT INTO source_files VALUES ",
        "(109,'/sources/mapping','mapping','uniprot_idmapping','annotation',",
        "10,1720000000000000001,'imported',NULL)"))
    edit_test_vault(path, paste0(
        "INSERT INTO protein_identifiers VALUES ",
        "(61,109,'P1','Gene_Name','GENE1'),",
        "(61,109,'P1','Gene_Name','GENE1_ALT'),",
        "(62,109,'P2','Gene_Name','GENE2')"))
    vault <- open_vault(path)
    on.exit(close_vault(vault), add = TRUE)
    precursor <- as_summarized_experiment(vault, level = "precursor",
        mapping_types = "Gene_Name")
    expect_equal(rownames(precursor), c("51", "52"))
    expect_equal(as.character(SummarizedExperiment::rowData(precursor)$sequence),
        c("PEPTIDE", "OTHER"))
    expect_equal(SummarizedExperiment::assay(precursor)["52", "S2"], 20)
    expect_equal(SummarizedExperiment::rowData(precursor)$protein_group_names[[1L]],
        "P1;P2")
    expect_equal(sort(SummarizedExperiment::rowData(precursor)$member_accessions[[1L]]),
        c("P1", "P2"))
    expect_equal(SummarizedExperiment::rowData(precursor)$leading_protein_id[[1L]], "61")
    expect_equal(SummarizedExperiment::rowData(precursor)$leading_gene_name[[1L]],
        "GENE1;GENE1_ALT")
    expect_equal(SummarizedExperiment::rowData(precursor)$member_gene_name[[1L]],
        "GENE1;GENE1_ALT;GENE2")
    expect_error(as_summarized_experiment(vault, level = "peptide"),
        "peptide_aggregation is required")
    peptide <- as_summarized_experiment(vault, level = "peptide",
        peptide_aggregation = "sum")
    expect_equal(SummarizedExperiment::assay(peptide)["31", "S1"], 10)
    expect_equal(SummarizedExperiment::assay(peptide)["32", "S2"], 20)
    expect_equal(SummarizedExperiment::rowData(peptide)$protein_group_names[[1L]],
        "P1;P2")
    expect_equal(SummarizedExperiment::rowData(peptide)$leading_protein_id[[1L]], "61")
    expect_equal(SummarizedExperiment::rowData(peptide)$member_gene_name[[1L]],
        "GENE1;GENE1_ALT;GENE2")
    expect_equal(S4Vectors::metadata(peptide)$peptide_aggregation, "sum")
    expect_true(is.na(SummarizedExperiment::rowData(precursor)$member_gene_name[[2L]]))
})

test_that("peptide aggregation combines distinct precursors after sample resolution", {
    path <- make_test_vault()
    edit_test_vault(path,
        "INSERT INTO peptidoforms VALUES (43,31,'PEPTIDE[Oxidation]','DIA-NN')")
    edit_test_vault(path, "INSERT INTO precursors VALUES (53,43,2)")
    edit_test_vault(path,
        "INSERT INTO precursor_quant VALUES (53,111,102,1,4,NULL,NULL)")
    vault <- open_vault(path)
    on.exit(close_vault(vault), add = TRUE)

    summed <- as_summarized_experiment(vault, level = "peptide",
        peptide_aggregation = "sum")
    highest <- as_summarized_experiment(vault, level = "peptide",
        peptide_aggregation = "max")
    averaged <- as_summarized_experiment(vault, level = "peptide",
        peptide_aggregation = "mean")
    expect_equal(SummarizedExperiment::assay(summed)["31", "S1"], 14)
    expect_equal(SummarizedExperiment::assay(highest)["31", "S1"], 10)
    expect_equal(SummarizedExperiment::assay(averaged)["31", "S1"], 7)
    expect_equal(length(SummarizedExperiment::rowData(summed)$precursor_ids[[1L]]), 2L)
})
