test_that("UniProt annotations remain lazy, one-to-many, and source-specific", {
    path <- make_test_vault()
    edit_test_vault(path, paste0(
        "INSERT INTO source_files VALUES ",
        "(109,'/sources/human.idmapping.gz','uniprot_idmapping/human.idmapping.gz',",
        "'uniprot_idmapping','annotation',10,1720000000000000001,'imported',NULL)"
    ))
    edit_test_vault(path, paste0(
        "INSERT INTO protein_identifiers VALUES ",
        "(61,109,'O75083','Gene_Name','WDR1'),",
        "(61,109,'O75083','Ensembl','ENSG1'),",
        "(61,109,'O75083','Ensembl','ENSG2')"
    ))
    vault <- open_vault(path)
    on.exit(close_vault(vault), add = TRUE)
    expect_true("uniprot-idmapping" %in% vault_info(vault)$layers)
    query <- tbl_protein_identifiers(vault, mapping_type = "Ensembl")
    expect_s3_class(query, "tbl_lazy")
    rows <- dplyr::collect(query)
    expect_equal(sort(rows$mapped_id), c("ENSG1", "ENSG2"))
    expect_equal(unique(rows$accession), "P1")
    expect_equal(unique(rows$source_accession), "O75083")
    expect_equal(nrow(dplyr::collect(tbl_protein_identifiers(vault, source_id = 109))), 3L)
    expect_equal(nrow(dplyr::collect(tbl_protein_identifiers(vault, source_id = 999))), 0L)
})

test_that("annotations explain how to add the missing layer", {
    vault <- open_vault(make_test_vault())
    on.exit(close_vault(vault), add = TRUE)
    expect_error(tbl_protein_identifiers(vault), "add-uniprot-idmapping")
    expect_error(tbl_vault(vault, "protein_identifiers"), "add-uniprot-idmapping")
})
