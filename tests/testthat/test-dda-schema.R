test_that("DDA peptide and TMT schema are queryable", {
    path <- make_test_vault()
    edit_test_vault(path, paste(
        "UPDATE analysis SET pipeline='metapipeline-ProteomicDDA'"
    ))
    edit_test_vault(path, "DELETE FROM protein_quant")
    edit_test_vault(path, "DELETE FROM precursor_quant")
    edit_test_vault(path, paste(
        "INSERT INTO source_files VALUES",
        "(110,'/sources/protein','protein','openms_protein_quant','quantitation',10,1,'imported',NULL),",
        "(111,'/sources/peptide','peptide','openms_peptide_quant','quantitation',10,1,'imported',NULL),",
        "(112,'/sources/channels','channels','tmt_channel_manifest','design',10,1,'imported',NULL),",
        "(113,'/sources/consensus','consensus','tmt_consensus_scan','quantitation',10,1,'imported',NULL)"
    ))
    edit_test_vault(path, "INSERT INTO quant_methods VALUES (3,'ProteinQuantifier')")
    edit_test_vault(path, "INSERT INTO quant_methods VALUES (4,'TMT.ScanBased')")
    edit_test_vault(path, "INSERT INTO grouping_methods VALUES (3,'ProteinQuantifier')")
    edit_test_vault(path, paste(
        "INSERT INTO protein_groups ",
        "(group_id,source_id,source_group_id,method_id,group_type,score,score_type,",
        "reported_n_proteins,reported_n_peptides) VALUES",
        "(95,110,'P1/P2',3,'quantification',0.8,'protein_score',2,1)"
    ))
    edit_test_vault(path, "INSERT INTO protein_group_members VALUES (95,61),(95,62)")
    edit_test_vault(path, "INSERT INTO samples VALUES (14,'plex_01|126'),(15,'plex_01|127N')")
    edit_test_vault(path, "INSERT INTO runs VALUES (24,NULL,'fraction_1','/input/fraction_1.mzML')")
    edit_test_vault(path, "INSERT INTO quant_channels VALUES (114,14,NULL,'TMT'),(115,15,NULL,'TMT')")
    edit_test_vault(path, "INSERT INTO tmt_plexes VALUES (201,'plex_01',114)")
    edit_test_vault(path, "INSERT INTO tmt_channels VALUES (114,201,'126','reference'),(115,201,'127N','treated')")
    edit_test_vault(path, "INSERT INTO tmt_runs VALUES (24,201)")
    edit_test_vault(path, "INSERT INTO spectra VALUES (72,24,'scan=8',130,2,510,2,104)")
    edit_test_vault(path, "INSERT INTO psms VALUES (83,72,51,'MSFragger',1,'target',104,2,0,1)")
    edit_test_vault(path, paste(
        "INSERT INTO precursor_quant ",
        "(quant_id,precursor_id,channel_id,source_id,method_id,quantity,",
        "source_element_id,source_map_id,rt_seconds,spectrum_id) VALUES",
        "(6,51,114,113,4,10,'e_1',0,120,72),",
        "(7,51,114,113,4,20,'e_2',0,180,72),",
        "(8,51,115,113,4,100,'e_1',1,120,72),",
        "(9,51,115,113,4,1,'e_2',1,180,72)"
    ))
    edit_test_vault(path, "INSERT INTO protein_quant VALUES (95,114,110,3,100),(95,115,110,3,200)")
    edit_test_vault(path, paste(
        "INSERT INTO peptide_quant_features VALUES",
        "(1001,111,1,31,'PEPTIDE(TMTpro)','P1/P2',2,0,'all')"
    ))
    edit_test_vault(path, "INSERT INTO peptide_quant VALUES (1001,114,111,3,10),(1001,115,111,3,20)")
    vault <- open_vault(path)
    on.exit(close_vault(vault), add = TRUE)

    expect_true(all(c("tmt_plexes", "tmt_channels", "tmt_runs",
        "peptide_quant_long") %in% list_vault_tables(vault)))
    expect_equal(dplyr::collect(tbl_protein_quant(vault))$quantity, c(100, 200))
    expect_equal(dplyr::collect(tbl_protein_groups(vault,
        method = "ProteinQuantifier"))$reported_n_peptides[[1]], 1)
    expect_equal(dplyr::collect(tbl_peptide_quant(vault))$source_peptide,
        rep("PEPTIDE(TMTpro)", 2))
    precursor_rows <- dplyr::collect(dplyr::arrange(tbl_precursor_quant(vault), quant_id))
    expect_equal(precursor_rows$quantity, c(10, 20, 100, 1))
    expect_equal(precursor_rows$rt_seconds, c(120, 180, 120, 180))
    expect_equal(precursor_rows$spectrum_id, bit64::as.integer64(rep(72, 4)))
    precursor_matrix <- matrix_precursor_abundance(vault)
    expect_equal(unname(precursor_matrix["51", c("plex_01|126", "plex_01|127N")]),
        c(30, 101))
    expect_equal(matrix_precursor_abundance(vault, duplicates = "sum"),
        precursor_matrix)
    expect_error(matrix_precursor_abundance(vault, duplicates = "best"),
        "TMT precursor duplicates")
    expect_error(matrix_precursor_abundance(vault, duplicates = "mean"),
        "TMT precursor duplicates")
    expect_error(matrix_precursor_abundance(vault, duplicates = "error"),
        "Duplicate feature/sample")
    precursor_se <- as_summarized_experiment(vault, level = "precursor")
    expect_equal(as.numeric(SummarizedExperiment::assay(precursor_se)["51", ]),
        c(30, 101))
    expect_error(as_summarized_experiment(vault, level = "precursor",
        duplicates = "best"), "TMT precursor duplicates")
    expect_equal(matrix_peptide_abundance(vault)["1001", "plex_01|127N"], 20)
    expect_equal(summarize_quantified_counts(vault, level = "peptide")$quantified_count,
        c(1, 1))
    tmt_psm <- dplyr::collect(dplyr::filter(tbl_psms(vault), native_id == "scan=8"))
    expect_equal(tmt_psm$sample_name, "plex_01")

    se <- as_summarized_experiment(vault, level = "peptide",
        samples = c("plex_01|126", "plex_01|127N"))
    expect_equal(SummarizedExperiment::assay(se)["1001", "plex_01|126"], 10)
    rows <- as.data.frame(SummarizedExperiment::rowData(se))
    columns <- as.data.frame(SummarizedExperiment::colData(se))
    expect_equal(rows$source_peptide, "PEPTIDE(TMTpro)")
    expect_equal(rows$leading_protein_accession, "P1")
    expect_equal(columns$tmt_plex, rep("plex_01", 2))
    expect_equal(columns$tmt_reporter, c("126", "127N"))
    expect_equal(columns$tmt_reference, c(TRUE, FALSE))
    expect_equal(columns$run_names[[2]], "fraction_1")
    expect_s3_class(plot_protein_quantity_distribution(vault), "ggplot")
    expect_s3_class(plot_protein_abundance(vault, "P1/P2", style = "ggplot2"), "ggplot")
})
