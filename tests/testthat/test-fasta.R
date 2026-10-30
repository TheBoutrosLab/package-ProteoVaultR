test_that("reference FASTA is reported as an imported layer", {
    path <- make_test_vault()
    edit_test_vault(
        path = path,
        statement = paste(
            "INSERT INTO source_files VALUES",
            "(109,'/reference/search.fasta','reference_fasta/search.fasta',",
            "'fasta','reference',100,1720000000000000001,'imported',NULL)"
        )
    )
    vault <- open_vault(path = path)
    on.exit(close_vault(vault = vault), add = TRUE)
    expect_true("fasta" %in% vault_info(vault = vault)$layers)
})
