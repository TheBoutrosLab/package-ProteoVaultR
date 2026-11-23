test_that("QC plots read selected metrics from the vault", {
    path <- make_test_vault()
    edit_test_vault(
        path = path,
        statement = paste(
            "INSERT INTO source_files VALUES",
            "(109,'/sources/qc','qc','qc_metrics','qc',10,1720000000000000001,'imported',NULL)"
        )
    )
    rows <- data.frame(
        metric_name = c(
            "scan_count",
            paste0("mass_error_", c("p05", "p25", "p50", "p75", "p95")),
            paste0("feature_charge_", 1:4, "_psm_count")
        ),
        value = c("100", "1", "2", "3", "4", "5", "10", "20", "30", "40"),
        source_process = c(rep("acquisition", 6), rep("MSBooster", 4))
    )
    statements <- vapply(
        X = seq_len(nrow(rows)),
        FUN = function(index) {
            sprintf(
                paste0(
                    "(109,%d,'TEST','identification','%s','S1','MSFragger','',",
                    "'%s','all','%s','count','available','S1.idXML')"
                ),
                index,
                rows$source_process[[index]],
                rows$metric_name[[index]],
                rows$value[[index]]
            )
        },
        FUN.VALUE = character(1)
    )
    edit_test_vault(
        path = path,
        statement = paste(
            "INSERT INTO qc_metrics ",
            "(source_id,row_number,dataset_id,layer,source_process,sample_name,",
            "search_engine,fraction_id,metric_name,\"group\",value,unit,status,source_file) VALUES",
            paste(statements, collapse = ",")
        )
    )
    edit_test_vault(path, paste(
        "INSERT INTO qc_metrics ",
        "(source_id,row_number,dataset_id,layer,source_process,sample_name,",
        "search_engine,fraction_id,metric_name,\"group\",value,unit,status,source_file) ",
        "SELECT source_id, row_number+100, dataset_id,",
        "layer, source_process, sample_name, search_engine, fraction_id,",
        "metric_name, 'target', value, unit, status, source_file",
        "FROM qc_metrics WHERE metric_name LIKE 'mass_error_%'"
    ))
    vault <- open_vault(path = path)
    on.exit(close_vault(vault = vault), add = TRUE)

    metric_bpg <- plot_qc_metric(vault = vault, metric_name = "scan_count")
    metric_ggplot <- plot_qc_metric(
        vault = vault,
        metric_name = "scan_count",
        style = "ggplot2"
    )
    box_ggplot <- plot_qc_boxplot(vault = vault, metric_name = "mass_error")
    box_bpg <- plot_qc_boxplot(
        vault = vault,
        metric_name = "mass_error",
        style = "BPG"
    )
    box_target <- plot_qc_boxplot(vault, "mass_error", levels = "target")
    box_ordered <- plot_qc_boxplot(vault, "mass_error",
        levels = c("target", "all"))
    box_target_bpg <- plot_qc_boxplot(vault, "mass_error",
        levels = "target", style = "BPG")
    charge_bpg <- plot_qc_charge_psm_counts(vault = vault)
    charge_ggplot <- plot_qc_charge_psm_counts(
        vault = vault,
        style = "ggplot2"
    )
    expect_s3_class(metric_bpg, "trellis")
    expect_s3_class(metric_ggplot, "ggplot")
    expect_s3_class(box_ggplot, "ggplot")
    expect_s3_class(box_bpg, "trellis")
    expect_equal(as.character(box_target$data$group), "target")
    expect_equal(levels(box_target$data$group), "target")
    expect_equal(levels(box_ordered$data$group), c("target", "all"))
    expect_s3_class(box_target_bpg, "trellis")
    expect_s3_class(charge_bpg, "trellis")
    expect_s3_class(charge_ggplot, "ggplot")
    grDevices::pdf(file = tempfile(fileext = ".pdf"))
    on.exit(grDevices::dev.off(), add = TRUE)
    expect_silent(print(box_bpg))
    expect_silent(print(box_target_bpg))
    expect_error(plot_qc_boxplot(vault, "mass_error", levels = "missing"),
        "levels must select")
    expect_error(
        plot_qc_metric(vault = vault, metric_name = "missing"),
        "No available QC rows"
    )
})

test_that("TMT QC channels are individual bars across plexes", {
    path <- make_test_vault()
    edit_test_vault(path, paste(
        "INSERT INTO source_files VALUES",
        "(109,'/sources/qc','qc','qc_metrics','qc',10,1720000000000000001,'imported',NULL)"
    ))
    edit_test_vault(path, paste(
        "INSERT INTO samples VALUES",
        "(14,'plex_01|126'),(15,'plex_01|131'),",
        "(16,'plex_02|126'),(17,'plex_02|127N')"
    ))
    edit_test_vault(path, paste(
        "INSERT INTO quant_channels VALUES",
        "(114,14,NULL,'TMT'),(115,15,NULL,'TMT'),",
        "(116,16,NULL,'TMT'),(117,17,NULL,'TMT')"
    ))
    edit_test_vault(path, paste(
        "INSERT INTO tmt_plexes VALUES",
        "(201,'plex_01',114),(202,'plex_02',116)"
    ))
    edit_test_vault(path, paste(
        "INSERT INTO tmt_channels VALUES",
        "(114,201,'126','reference'),(115,201,'131',NULL),",
        "(116,202,'126','reference'),(117,202,'127N','treated')"
    ))
    edit_test_vault(path, paste(
        "INSERT INTO qc_metrics",
        "(source_id,row_number,dataset_id,layer,source_process,sample_name,",
        "search_engine,fraction_id,metric_name,\"group\",value,unit,status,source_file) VALUES",
        "(109,1,'TEST','quantification','','plex_01','','','feature_count',",
        "'reference','100','count','available','plex_01.csv'),",
        "(109,2,'TEST','quantification','','plex_01','','','feature_count',",
        "'','150','count','available','plex_01.csv'),",
        "(109,3,'TEST','quantification','','plex_02','','','feature_count',",
        "'reference','200','count','available','plex_02.csv'),",
        "(109,4,'TEST','quantification','','plex_02','','','feature_count',",
        "'plex_02|127N','250','count','available','plex_02.csv')"
    ))
    vault <- open_vault(path)
    on.exit(close_vault(vault), add = TRUE)

    ggplot <- plot_qc_metric(vault, "feature_count", style = "ggplot2")
    bpg <- plot_qc_metric(vault, "feature_count", style = "BPG")
    expect_s3_class(ggplot, "ggplot")
    expect_s3_class(bpg, "trellis")
    expect_setequal(as.character(ggplot$data$sample), c(
        "plex_01|126", "plex_01|131", "plex_02|126", "plex_02|127N"
    ))
    expect_equal(rev(levels(ggplot$data$sample)), c(
        "plex_01|126", "plex_01|131", "plex_02|126", "plex_02|127N"
    ))
    expect_equal(unique(as.character(ggplot$data$group)), "all")
    expect_equal(ggplot$data$value, c(100, 150, 200, 250))
    one_plex <- plot_qc_metric(vault, "feature_count",
        source_file_pattern = "plex_01.csv", style = "ggplot2")
    one_plex_bpg <- plot_qc_metric(vault, "feature_count",
        source_file_pattern = "plex_01.csv", style = "BPG")
    expect_setequal(as.character(one_plex$data$sample),
        c("reference", "131"))
    expect_equal(rev(levels(one_plex$data$sample)), c("reference", "131"))
    expect_s3_class(one_plex_bpg, "trellis")
    grDevices::pdf(file = tempfile(fileext = ".pdf"))
    on.exit(grDevices::dev.off(), add = TRUE)
    expect_silent(print(bpg))
    expect_silent(print(one_plex_bpg))
})
