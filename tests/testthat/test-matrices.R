test_that("matrices preserve values, identifiers, and requested sample order", {
    vault <- open_vault(make_test_vault())
    on.exit(close_vault(vault), add = TRUE)
    expect_message(matrix_protein_abundance(vault), "Collapsed 0 duplicated source_group_id")
    result <- matrix_protein_abundance(vault)
    expect_equal(dim(result), c(2L, 2L))
    expect_equal(rownames(result), c("P1;P2", "P3"))
    expect_equal(colnames(result), c("S1", "S2"))
    expect_equal(result["P1;P2", "S1"], 100)
    expect_true(is.na(result["P1;P2", "S2"]))
    expect_equal(result["P3", "S1"], 0)
    internal <- matrix_protein_abundance(vault, row_id = "group_id")
    expect_equal(rownames(internal), c("91", "92"))
    expect_equal(internal["91", "S1"], result["P1;P2", "S1"])
    expect_equal(colnames(matrix_protein_abundance(vault, samples = c("S2", "S1"))), c("S2", "S1"))
    expect_error(matrix_protein_abundance(vault, samples = "unknown"), "sample names")
    expect_error(matrix_protein_abundance(vault, samples = c("S1", "S1")), "unique")
    expect_error(matrix_protein_abundance(vault, row_id = "accession"), "arg")
    precursor <- matrix_precursor_abundance(vault)
    expect_equal(precursor["52", "S2"], 20)
    expect_true(is.na(precursor["51", "S2"]))
})

test_that("duplicate policies resolve each source group and sample", {
    path <- make_test_vault()
    edit_test_vault(
        path = path,
        statement = paste0(
            "INSERT INTO protein_groups ",
            "(group_id,source_id,source_group_id,method_id,group_type,score,score_type) VALUES ",
            "(95,101,'P1;P2',1,'quantification',NULL,NULL)"
        )
    )
    edit_test_vault(
        path = path,
        statement = paste0(
            "INSERT INTO protein_quant VALUES ",
            "(95,111,101,1,150),(95,112,101,1,NULL)"
        )
    )
    vault <- open_vault(path)
    on.exit(close_vault(vault), add = TRUE)
    expect_message(
        matrix_protein_abundance(vault),
        "Collapsed 1 duplicated source_group_id ID\\(s\\) using best"
    )
    expect_equal(matrix_protein_abundance(vault)["P1;P2", "S1"], 150)
    expect_true(is.na(matrix_protein_abundance(vault)["P1;P2", "S2"]))
    expect_equal(matrix_protein_abundance(vault, duplicates = "sum")["P1;P2", "S1"], 250)
    expect_equal(matrix_protein_abundance(vault, duplicates = "mean")["P1;P2", "S1"], 125)
    expect_error(matrix_protein_abundance(vault, duplicates = "error"), "Duplicate feature/sample")
    internal <- matrix_protein_abundance(vault, row_id = "group_id")
    expect_equal(rownames(internal), c("91", "92", "95"))
    expect_equal(internal["95", "S1"], 150)
})

test_that("duplicate policy also applies to precursor matrices", {
    path <- make_test_vault()
    edit_test_vault(path, "INSERT INTO quant_channels VALUES (114,11,21,'LFQ')")
    edit_test_vault(path, paste0(
        "INSERT INTO precursor_quant ",
        "(quant_id,precursor_id,channel_id,source_id,method_id,quantity,q_value,pep) ",
        "VALUES (6,51,114,102,1,12,NULL,NULL)"
    ))
    vault <- open_vault(path)
    on.exit(close_vault(vault), add = TRUE)
    expect_message(
        matrix_precursor_abundance(vault),
        "Collapsed 1 duplicated precursor_id ID\\(s\\) using best"
    )
    expect_equal(matrix_precursor_abundance(vault)["51", "S1"], 12)
    expect_equal(matrix_precursor_abundance(vault, duplicates = "sum")["51", "S1"], 22)
    expect_error(matrix_precursor_abundance(vault, duplicates = "error"), "Duplicate feature/sample")
})
