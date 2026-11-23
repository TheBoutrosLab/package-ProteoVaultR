# Contributing

Document public functions with roxygen2 and regenerate NAMESPACE/man files with `roxygen2::roxygenise()`. Run `testthat::test_local()` and R CMD check before proposing changes.

Use four-space indentation, `%>%`, and named arguments in multiline calls. Preserve lazy queries until an explicit collection, matrix, or summary operation. Keep source/method/stage distinctions and protein-group memberships explicit. Analysis functions must not save figures or rewrite input databases.

## Preview the documentation

Install the package dependencies and pkgdown (`install.packages("pkgdown")`),
plus Pandoc (available with RStudio).
From the package directory, build the website with:

```r
pkgdown::build_site(preview = FALSE)
pkgdown::preview_site()
```

Generated pages are written to `docs/site/`, which is ignored by Git. README is
the homepage; navigation and the grouped function reference are configured in
`_pkgdown.yml`. The demo-data article shows download code with `eval = FALSE`;
rendering it does not download the databases. Add documentation `.Rmd` sources
directly under `vignettes/`, then add
their navigation entries to `_pkgdown.yml` and the README table of contents as
each page is ready. List website-only
sources individually in `.Rbuildignore` so package builds exclude them.
Number article section headings as `1.`, `2.`, and so on, with subsections such as
`1.1.`. Keep page titles unnumbered and leave code blocks unchanged.

The pkgdown GitHub Actions workflow caches the HepG2 DIA database with a key
based on the dataset manifest and downloads from public S3 when cached data is
missing or invalid. The verified input is cached before dependency installation
and rendering.
Publish only the generated site, not the databases or local planning files.
The local documentation plan is kept in ignored `docs/_build/`.

For a terminal build on this Mac, RStudio's bundled Pandoc can be exposed with:

```bash
export RSTUDIO_PANDOC=/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64
Rscript -e 'pkgdown::build_site(preview = FALSE)'
```

Use your own Pandoc installation path on other systems. The `_pkgdown.yml` `url`
field and the published website link in `DESCRIPTION` should be added together
when the hosting address is confirmed. Until then, pkgdown reports missing URL
metadata during its site readiness checks; this does not prevent a local build.
The README currently uses site-relative article links for the pkgdown homepage.
Once hosting is confirmed, replace them with published article URLs so the links
also work when the README is viewed on GitHub.

Use Bootstrap alert boxes sparingly for tips, notes, and consequential choices.
Use `alert-info` for tips or context and `alert-warning` for pitfalls. Include a
short bold label so the meaning does not depend on color. For example:

```markdown
::: {.alert .alert-info role="note"}
**Tip:** Filter results before bringing them into R.
:::
```

Callout labels are not headings and do not receive section numbers. Keep ordinary
explanations in the body text rather than putting every paragraph in a box.

Execute demo-analysis examples during article builds. Keep installation and
large-file download commands unevaluated. Give each visible example chunk one
result: one printed object, table, or plot. Keep library loading and other setup
in a hidden chunk. Display result tables with `knitr::kable()` and descriptive
column names rather than raw console tibbles. Omit table captions. Print lists, vectors,
and other non-tabular outputs in their native form; do not convert them into
data frames just for display. Pipe tables into `knitr::kable()` with `%>%`,
placing the formatting call after the data-selection steps.

For local builds, place `hepg2_ifna-dia.proteovault.duckdb` in
`vignettes/demo-data/`, or set `PROTEOVAULTR_DEMO_DIR` to an absolute path
containing that file before building. The current analysis articles use only
HepG2 DIA; the download guide also offers TMT and DDA LFQ data for readers.
Analysis articles must fail when their required input is missing; do not substitute
synthetic data or silently skip execution. The Actions workflow restores these files to the configured directory before
rendering. The local data
directory is excluded from both Git and package archives.

## Documentation workflow

`.github/workflows/pkgdown.yaml` builds on pull requests, pushes to `main`, and
manual dispatch. It uses the maintained r-lib setup actions and caches R
dependencies separately from the DIA demo database. Analysis code executes on
every build. Pull requests build without publishing or saving demo-data caches;
every successful build uploads `docs/site/` as a Pages artifact. A separate
job publishes that artifact only from `main`, using GitHub's official Pages
actions. The build job has read-only repository permissions.

`.github/demo-data.json` records the inputs required by executable articles,
with their public URLs, sizes, and SHA-256 checksums.
When replacing a dataset, update its size and checksum together; changing the
manifest gives the required inputs a new cache key. If future articles need
another database, add it to the manifest and the workflow cache paths. To prepare inputs locally from the repository root:

```bash
python3 .github/scripts/prepare-demo-data.py
```

The script verifies existing files before reusing them, downloads invalid or
missing files into temporary paths, and replaces them only after verification.
A cache hit never bypasses checksum verification. GitHub caches are immutable;
if a cache entry is repeatedly invalid, delete that entry in Actions caches and
rerun so a verified copy can be saved under the same key.

GitHub Pages is not yet enabled for this repository. Before merging, configure
Settings → Pages → Build and deployment → Source to use **GitHub Actions**
and confirm the intended site visibility. This setting is needed for publishing;
pull-request builds do not require it. Then set the actual site URL in
`_pkgdown.yml` and `DESCRIPTION`. Verify both a cold build and a cache-hit build
on GitHub before treating deployment as validated.

Use `theme_boutroslab()` for ggplot2 figures in documentation. It follows the
existing BPG-style plots: white panels without grids, black borders and ticks,
bold black 12-point titles, and 11-point axis/legend labels. Add plot-specific
settings such as label rotation after the shared theme. Data annotations should
also use 11-point bold text where practical.

For package plotting helpers, lead with a single-line call containing only the
required arguments. Let the functions supply their defaults and styling. Add
optional arguments only when the example specifically needs to demonstrate them.
