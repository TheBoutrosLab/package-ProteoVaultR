# ProteoVaultR

ProteoVaultR reads linked proteomics results from ProteoVault DuckDB files. It provides read-only connections, lazy queries, quantitative matrices, summaries, and QC plots.

## Install

Requires R 4.2 or newer, DuckDB 2.0 or newer, and dbplyr 2.6 or newer. Install the latest version from GitHub:

```r
install.packages("remotes")
remotes::install_github("TheBoutrosLab/package-ProteoVaultR")
```

## Quick start

```r
library(ProteoVaultR)

vault <- open_vault(path = "results.proteovault.duckdb")
vault_info(vault)
counts <- summarize_quantified_counts(vault = vault)
abundance <- matrix_protein_abundance(vault = vault)
genes <- dplyr::collect(tbl_protein_identifiers(vault, mapping_type = "Gene_Name"))
close_vault(vault)
```

`tbl_protein_identifiers()` reads UniProt annotations when the optional
`uniprot-idmapping` layer is present. It retains multiple mappings per protein
and the mapping archive's `source_id`.

Attach sample metadata to the R handle to show shorter IDs on quantity plots.
The `sample_name` column must match the names stored in the vault. Assign the
returned handle; the DuckDB file is unchanged.

```r
vault <- open_vault(path = "results.proteovault.duckdb")
vault <- set_sample_metadata(vault, data.frame(
    sample_name = c("long_raw_sample_name_1", "long_raw_sample_name_2"),
    sample_id = c("Control", "Treatment")
))
plot_protein_quantity_distribution(vault)
close_vault(vault)
```

See the vignettes for [getting started](vignettes/getting-started.Rmd), [querying results](vignettes/querying-results.Rmd), [optional layers](vignettes/optional-layers.Rmd), [protein coverage](vignettes/protein-coverage.Rmd), and [quality control](vignettes/qc.Rmd). These can be incorporated into a bookdown site later. Function-specific help is available with `?open_vault` and similar R help topics.

## License

GNU GPL version 2 or later. Copyright 2026 Boutros Lab, Sanford Burnham Prebys Medical Discovery Institute.
