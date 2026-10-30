# Contributing

Document public functions with roxygen2 and regenerate NAMESPACE/man files with `roxygen2::roxygenise()`. Run `testthat::test_local()` and R CMD check before proposing changes.

Use four-space indentation, `%>%`, and named arguments in multiline calls. Preserve lazy queries until an explicit collection, matrix, or summary operation. Keep source/method/stage distinctions and protein-group memberships explicit. Analysis functions must not save figures or rewrite input databases.
