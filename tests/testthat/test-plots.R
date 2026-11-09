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
        statement = paste("INSERT INTO qc_metrics (source_id,row_number,dataset_id,layer,source_process,sample_name,search_engine,fraction_id,metric_name,\"group\",value,unit,status,source_file) VALUES", paste(statements, collapse = ","))
    )
    edit_test_vault(path, paste(
        "INSERT INTO qc_metrics (source_id,row_number,dataset_id,layer,source_process,sample_name,search_engine,fraction_id,metric_name,\"group\",value,unit,status,source_file) SELECT source_id, row_number+100, dataset_id,",
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
