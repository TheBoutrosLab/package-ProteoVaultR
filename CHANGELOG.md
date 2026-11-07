# Changelog

All notable changes to ProteoVaultR are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased] - 2026-11-04

### Added

- Fit limma differential-abundance contrasts on exported assays with validated designs and annotated results.
- Plot limma contrasts as BPG or ggplot volcano plots with transformed probability axes and selected protein annotations.
- Export protein-group, precursor, and aggregated peptide quantities with sample and feature annotations as a `SummarizedExperiment`.
- List the tables and views present in a vault with `list_vault_tables()`.
- Plot a selected protein group's abundance across samples with `plot_protein_abundance()`, optionally coloring bars by attached sample metadata.

### Changed

- Select exact protein-group or precursor quantification methods in quantity distribution plots.
- Plot protein-group and precursor quantity distributions without log transformation when requested.
- Draw visible black borders around bars in BPG and ggplot QC plots.
- Select and order QC boxplot group levels with `plot_qc_boxplot(levels = ...)`.
- Read DIA-NN report precursor quantities from core vaults without a `quant-details` layer.
- Attach sample metadata to vault objects in memory and use sample IDs in quantity plot labels.
- Read ProteoVault schema 1.1.0 and expose UniProt protein annotations through `tbl_protein_identifiers()`.

- Use source group labels as default protein matrix row names and allow selecting internal group IDs and duplicate handling.

## [0.1.0] - 2026-10-30

### Added

- Read-only access to ProteoVault schema 1.0.0 with lazy queries for PSMs, scores, protein groups, and quantification results.
- Protein-group and precursor abundance matrices, quantified-feature summaries, and source selection.
- QC metric access and plots for charge counts, metric distributions, precursor quantities, and protein-group quantities.
- Protein sequence coverage plots with validated peptide and PTM annotations.
- R help pages, vignettes, a synthetic DuckDB fixture, tests, and R CMD check CI.


[0.1.0]: https://github.com/TheBoutrosLab/package-ProteoVaultR/releases/tag/v0.1.0
[unreleased]: https://github.com/TheBoutrosLab/package-ProteoVaultR/compare/v0.1.0...HEAD
