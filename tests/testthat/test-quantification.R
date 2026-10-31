test_that("quantification access stays lazy and separates report measurements", {
    vault <- open_vault(path = make_test_vault())
    on.exit(close_vault(vault), add = TRUE)
    expect_s3_class(tbl_protein_quant(vault), "tbl_lazy")
    expect_s3_class(tbl_precursor_quant(vault), "tbl_lazy")
    rows <- dplyr::collect(tbl_protein_quant(vault))
    expect_equal(nrow(rows), 4L)
    expect_s3_class(rows$group_id, "integer64")
    expect_equal(sum(is.na(rows$quantity)), 1L)
    expect_equal(sum(rows$quantity == 0, na.rm = TRUE), 1L)
    report <- dplyr::collect(tbl_precursor_quant(vault, method = "Precursor.Quantity"))
    expect_equal(report$quantity, 5)
    expect_equal(report$q_value, .001)
    expect_error(tbl_protein_quant(vault, method = "Precursor.Quantity"))
    expect_error(tbl_protein_quant(vault, source_id = 102), "compatible")
})

test_that("ambiguous sources require a selection even when one is empty", {
    path <- make_test_vault()
    edit_test_vault(
        path,
        paste0(
            "INSERT INTO source_files SELECT 108,'/other','other',",
            "format,stage,size_bytes,mtime_ns,status,details ",
            "FROM source_files WHERE source_id=101"
        )
    )
    vault <- open_vault(path)
    on.exit(close_vault(vault), add = TRUE)
    expect_error(tbl_protein_quant(vault), "Candidates:.*101.*101|Candidates:.*108")
    expect_equal(nrow(dplyr::collect(tbl_protein_quant(vault, source_id = 108))), 0L)
    expect_equal(nrow(dplyr::collect(tbl_protein_quant(vault, source_id = 101))), 4L)
})

test_that("per-sample counts distinguish positive, zero, and missing", {
    vault <- open_vault(make_test_vault())
    on.exit(close_vault(vault), add = TRUE)
    counts <- summarize_quantified_counts(vault)
    expect_equal(counts$sample_name, c("S1", "S2"))
    expect_equal(as.numeric(counts$quantified_count), c(1, 1))
    expect_equal(unique(counts$level), "protein_group")
    precursor_counts <- summarize_quantified_counts(vault, level = "precursor")
    expect_equal(as.numeric(precursor_counts$quantified_count), c(1, 1))
})


test_that("counts count features rather than repeated measurements", {
    path <- make_test_vault()
    edit_test_vault(path, "INSERT INTO quant_channels VALUES (114,11,21,'LFQ')")
    edit_test_vault(path, "INSERT INTO protein_quant VALUES (91,114,101,1,100)")
    vault <- open_vault(path)
    on.exit(close_vault(vault), add = TRUE)
    counts <- summarize_quantified_counts(vault)
    expect_equal(as.numeric(counts$quantified_count), c(1, 1))
    expect_equal(matrix_protein_abundance(vault)["P1;P2", "S1"], 100)
    expect_error(matrix_protein_abundance(vault, duplicates = "error"), "Duplicate feature/sample")
})


test_that("measurement resolution follows schema rows rather than a producer format name", {
    path <- make_test_vault()
    edit_test_vault(path, "UPDATE source_files SET format='another_quant_format' WHERE source_id=101")
    vault <- open_vault(path)
    on.exit(close_vault(vault), add = TRUE)
    expect_equal(nrow(dplyr::collect(tbl_protein_quant(vault))), 4L)
    expect_equal(dim(matrix_protein_abundance(vault)), c(2L, 2L))
})
