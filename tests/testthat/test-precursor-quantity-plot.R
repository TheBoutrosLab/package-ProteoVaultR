test_that("precursor quantity distributions use all measured precursors", {
    path <- make_test_vault()
    edit_test_vault(path = path, statement = paste(
        "UPDATE precursor_quant SET quantity=15",
        "WHERE precursor_id=51 AND channel_id=112 AND source_id=102"
    ))
    edit_test_vault(
        path = path,
        statement = paste(
            "INSERT INTO precursor_quant ",
            "(quant_id,precursor_id,channel_id,source_id,method_id,quantity,q_value,pep) VALUES",
            "(6,51,112,103,2,8,NULL,NULL),",
            "(7,52,111,103,2,12,NULL,NULL)"
        )
    )
    vault <- open_vault(path = path)
    on.exit(close_vault(vault = vault), add = TRUE)

    boxplot <- plot_precursor_quantity_distribution(
        vault = vault,
        quantity_type = "both",
        plot_type = "boxplot"
    )
    violin <- plot_precursor_quantity_distribution(
        vault = vault,
        quantity_type = "normalized",
        plot_type = "violin"
    )

    expect_s3_class(boxplot, "ggplot")
    expect_s3_class(violin, "ggplot")
    expect_equal(nrow(boxplot$data), 6L)
    expect_setequal(as.character(unique(boxplot$data$precursor_id)), c("51", "52"))
    expect_equal(nrow(violin$data), 3L)
    expect_equal(boxplot$data$log2_quantity, log2(boxplot$data$quantity))
    raw_plot <- plot_precursor_quantity_distribution(vault, transform = "none")
    expect_equal(nrow(raw_plot$data), 4L)
    expect_true(0 %in% raw_plot$data$quantity)
    expect_equal(raw_plot$data$plot_quantity, raw_plot$data$quantity)
    expect_equal(raw_plot$labels$y, "precursor quantity")
    expect_silent(ggplot2::ggplot_build(boxplot))
    expect_s3_class(violin$layers[[1L]]$geom, "GeomViolin")

    all_precursors <- plot_precursor_quantity_distribution(
        vault = vault,
        quantity_type = "raw"
    )
    expect_equal(nrow(all_precursors$data), 3L)
    expect_error(
        plot_precursor_quantity_distribution(
            vault = vault,
            plot_type = "scatter"
        ),
        "boxplot.*violin"
    )
})

test_that("precursor distribution selects an exact report method", {
    path <- make_test_vault()
    edit_test_vault(path, "INSERT INTO quant_methods VALUES (3,'Precursor.Translated')")
    edit_test_vault(path, "INSERT INTO source_quant_methods VALUES (103,'precursor',3)")
    edit_test_vault(path,
        paste0(
            "INSERT INTO precursor_quant ",
            "(quant_id,precursor_id,channel_id,source_id,method_id,quantity,q_value,pep) ",
            "VALUES (6,51,111,103,3,25,0.005,0.002)"
        ))
    vault <- open_vault(path)
    on.exit(close_vault(vault), add = TRUE)

    plot <- plot_precursor_quantity_distribution(vault,
        method = "Precursor.Translated")
    expect_equal(plot$data$quantity, 25)
    expect_equal(as.character(plot$data$quantity_type), "Precursor.Translated")
    expect_error(plot_precursor_quantity_distribution(vault,
        quantity_type = "raw", method = "Precursor.Translated"), "either")
    expect_error(plot_precursor_quantity_distribution(vault,
        method = "matrix", max_q_value = 0.05), "report")
})

test_that("normalized distribution works without report quantities", {
    path <- make_test_vault()
    edit_test_vault(
        path = path,
        statement = "DELETE FROM precursor_quant WHERE source_id=103"
    )
    edit_test_vault(
        path = path,
        statement = "DELETE FROM source_files WHERE source_id=103"
    )
    edit_test_vault(
        path = path,
        statement = paste(
            "UPDATE precursor_quant SET quantity=15",
            "WHERE precursor_id=51 AND channel_id=112 AND source_id=102"
        )
    )
    vault <- open_vault(path = path)
    on.exit(close_vault(vault = vault), add = TRUE)

    plot <- plot_precursor_quantity_distribution(vault = vault)
    expect_s3_class(plot, "ggplot")
    expect_equal(nrow(plot$data), 3L)
    expect_equal(plot$data$quantity, c(10, 15, 20))
    expect_error(
        plot_precursor_quantity_distribution(
            vault = vault,
            quantity_type = "raw"
        ),
        "compatible precursor source"
    )
    expect_error(
        plot_precursor_quantity_distribution(
            vault = vault,
            max_q_value = 0.05
        ),
        "compatible precursor source"
    )
})

test_that("q-value cutoff filters precursor quantities from the report", {
    path <- make_test_vault()
    edit_test_vault(
        path = path,
        statement = "INSERT INTO quant_methods VALUES (3,'Precursor.Normalised')"
    )
    edit_test_vault(path, "INSERT INTO source_quant_methods VALUES (103,'precursor',3)")
    edit_test_vault(
        path = path,
        statement = paste(
            "INSERT INTO precursor_quant ",
            "(quant_id,precursor_id,channel_id,source_id,method_id,quantity,q_value,pep) VALUES",
            "(6,51,111,103,3,6,0.005,NULL),",
            "(7,51,112,103,3,9,0.02,NULL),",
            "(8,52,111,103,3,14,0.5,NULL)"
        )
    )
    vault <- open_vault(path = path)
    on.exit(close_vault(vault = vault), add = TRUE)

    filtered <- plot_precursor_quantity_distribution(
        vault = vault,
        max_q_value = 0.05
    )
    all_reported <- plot_precursor_quantity_distribution(
        vault = vault,
        max_q_value = 1
    )
    expect_equal(filtered$data$quantity, c(6, 9))
    expect_equal(all_reported$data$quantity, c(6, 9, 14))
    expect_equal(unique(as.character(filtered$data$quantity_type)), "normalized")
    expect_error(
        plot_precursor_quantity_distribution(
            vault = vault,
            max_q_value = 1.1
        ),
        "between 0 and 1"
    )
})
