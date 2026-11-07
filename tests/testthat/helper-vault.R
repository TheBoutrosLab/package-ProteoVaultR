make_test_vault <- function(optional = FALSE, path = tempfile(fileext = ".duckdb")) {
    connection <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE), dbdir = path)
    on.exit(DBI::dbDisconnect(connection), add = TRUE)
    sql <- paste(readLines(testthat::test_path("fixtures", "schema-v1.sql")), collapse = "\n")
    for (statement in strsplit(sql, ";", fixed = TRUE)[[1L]]) {
        if (nzchar(trimws(statement))) DBI::dbExecute(connection, statement)
    }
    statements <- c(
        "INSERT INTO metadata VALUES ('format','ProteoVault'),('schema_version','1.1.0'),('producer_version','0.1.0')",
        "INSERT INTO analysis VALUES (1,'TEST','metapipeline-ProteomicDIA','1.3.0','/not/a/real/source/root')",
        "INSERT INTO samples VALUES (11,'S1'),(12,'S2'),(13,'S3')",
        paste0(
            "INSERT INTO runs VALUES (21,11,'S1','/input/S1.mzML'),",
            "(22,12,'S2','/input/S2.mzML'),(23,13,'S3','/input/S3.mzML')"
        ),
        "INSERT INTO peptides VALUES (31,'PEPTIDE'),(32,'OTHER')",
        "INSERT INTO peptidoforms VALUES (41,31,'PEPTIDE','DIA-NN'),(42,32,'OTHER','DIA-NN')",
        "INSERT INTO precursors VALUES (51,41,2),(52,42,3)",
        paste0(
            "INSERT INTO proteins VALUES (61,'P1','PEPTIDE','first',FALSE,FALSE),",
            "(62,'P2',NULL,'second',FALSE,FALSE),(63,'P3',NULL,'third',FALSE,FALSE)"
        )
    )
    for (statement in statements) DBI::dbExecute(connection, statement)
    sources <- list(
        s_pg = c("diann_pg_matrix", "quantitation"),
        s_pr = c("diann_pr_matrix", "quantitation"),
        s_report = c("diann_report", "quantitation"),
        s_valid = c("idxml", "validation"),
        s_infer = c("idxml", "inference")
    )
    if (optional) {
        sources$s_search <- c("idxml", "search")
        sources$s_pin <- c("pin", "MSBooster")
    }
    source_ids <- c(s_pg = 101, s_pr = 102, s_report = 103, s_valid = 104, s_infer = 105, s_search = 106, s_pin = 107)
    for (source in names(sources)) {
        DBI::dbExecute(
            connection,
            paste0(
                "INSERT INTO source_files VALUES (", source_ids[[source]], ",'/sources/", source,
                "','", source, "','", sources[[source]][1], "','", sources[[source]][2],
                "',10,1720000000000000001,'imported',NULL)"
            )
        )
    }
    statements <- c(
        "INSERT INTO quant_methods VALUES (1,'matrix'),(2,'Precursor.Quantity')",
        "INSERT INTO grouping_methods VALUES (1,'DIA-NN'),(2,'ProteinProphet')",
        "INSERT INTO processing_stages VALUES (1,'validation'),(2,'search')",
        "INSERT INTO score_types VALUES (1,'q-value'),(2,'Posterior Error Probability_score'),(3,'hyperscore')",
        "INSERT INTO spectra VALUES (71,21,'scan=7',120,2,500,2,104)",
        "INSERT INTO psms VALUES (81,71,51,'MSFragger',1,'target',104,1,0,1)",
        "INSERT INTO psm_scores VALUES (81,1,0.005,FALSE),(81,2,0.002,FALSE)",
        "INSERT INTO psm_proteins VALUES (81,61),(81,62)",
        "INSERT INTO peptide_proteins VALUES (31,61),(31,62)",
        paste0(
            "INSERT INTO protein_groups VALUES ",
            "(91,101,'P1;P2',1,'quantification',NULL,NULL),",
            "(92,101,'P3',1,'quantification',NULL,NULL),",
            "(93,105,'group_0',2,'protein_group',0.9,'probability'),",
            "(94,105,'indist_0',2,'indistinguishable_proteins',0.8,'probability')"
        ),
        "INSERT INTO protein_group_members VALUES (91,61),(91,62),(92,63),(93,61),(93,62),(94,61)",
        "INSERT INTO quant_channels VALUES (111,11,21,'LFQ'),(112,12,22,'LFQ'),(113,13,23,'LFQ')",
        "INSERT INTO protein_quant VALUES (91,111,101,1,100),(91,112,101,1,NULL),(92,111,101,1,0),(92,112,101,1,200)",
        paste0(
            "INSERT INTO precursor_quant VALUES ",
            "(51,111,102,1,10,NULL,NULL),(51,112,102,1,NULL,NULL,NULL),",
            "(52,111,102,1,0,NULL,NULL),(52,112,102,1,20,NULL,NULL),",
            "(51,111,103,2,5,0.001,0.002)"
        )
    )
    if (optional) {
        statements <- c(
            statements,
            "INSERT INTO psms VALUES (82,71,51,'MSFragger',1,'target',106,1,0,2)",
            "INSERT INTO psm_scores VALUES (82,3,25,TRUE)",
            "INSERT INTO psm_features VALUES (82,'feature',0.2)",
            "INSERT INTO pin_features VALUES (107,1,'S1.7.7.2','7',1,'deltaRTLOESS',0.2)"
        )
    }
    for (statement in statements) DBI::dbExecute(connection, statement)
    path
}

edit_test_vault <- function(path, statement) {
    connection <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE), dbdir = path)
    on.exit(DBI::dbDisconnect(connection), add = TRUE)
    DBI::dbExecute(connection, statement)
    invisible(path)
}
