test_that("cohort QC metrics are lazy and keep unavailable values", {
    path <- make_test_vault()
    vault <- open_vault(path = path)
    expect_error(tbl_qc_metrics(vault = vault), "--layer qc")
    close_vault(vault = vault)
    edit_test_vault(
        path = path,
        statement = paste(
            "INSERT INTO source_files VALUES",
            "(109,'/sources/qc','qc','qc_metrics','qc',10,1720000000000000001,'imported',NULL)"
        )
    )
    edit_test_vault(
        path = path,
        statement = paste(
            "INSERT INTO qc_metrics (source_id,row_number,dataset_id,layer,source_process,sample_name,search_engine,fraction_id,metric_name,\"group\",value,unit,status,source_file) VALUES",
            "(109,1,'TEST','acquisition','','S1','','','ms2_count','all','1200','count','available','S1.mzML'),",
            paste0(
                "(109,2,'TEST','raw_identification','','S1','MSFragger',",
                "'','score_type','all','expect','','available','S1.idXML'),"
            ),
            paste0(
                "(109,3,'TEST','protein_inference','','','','','group_count',",
                "'all',NULL,'count','unavailable','groups.idXML')"
            )
        )
    )
    vault <- open_vault(path = path)
    on.exit(close_vault(vault = vault), add = TRUE)
    expect_true("qc" %in% vault_info(vault = vault)$layers)
    query <- tbl_qc_metrics(vault = vault)
    expect_s3_class(query, "tbl_lazy")
    rows <- dplyr::collect(query)
    expect_equal(rows$metric_name, c("ms2_count", "score_type", "group_count"))
    expect_equal(rows$value[[1L]], "1200")
    expect_equal(rows$value[[2L]], "expect")
    expect_true(is.na(rows$value[[3L]]))
    expect_equal(rows$status[[3L]], "unavailable")
})
