
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cyCompare

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/cyCompare)](https://CRAN.R-project.org/package=cyCompare)
[![R-CMD-check](https://github.com/ggrlab/cyCompare/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ggrlab/cyCompare/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**cyCompare** is a computational toolkit for comparing flow cytometry
data across runs, instruments, and sites. Core capabilities include:

  - Consistent pre-processing (e.g., configurable `asinh`/cofactor
    transforms and channel/marker standardisation).
  - Cross-instrument/site alignment and comparability checks.
  - Clustering integration (e.g., FlowSOM) and downstream agreement
    analyses (e.g., Bland-Altman).
  - Variability quantification using optimal transport distance (OTD)
    for both device and sample variability.
  - Reporting utilities for cohort-level comparisons and QC.

> **Scope.** cyCompare is about computational normalisation and
> comparability for cytometry-like data. Many cells, low to moderate
> number of features per cell.

-----

## Installation

You can install the development version of **cyCompare** from GitHub:

``` r
# using pak (recommended)
install.packages("pak")
pak::pak("ggrlab/cyCompare")

# or using remotes
install.packages("remotes")
remotes::install_github("ggrlab/cyCompare")
```

If you plan to build vignettes locally:

``` r
remotes::install_github(
  "ggrlab/cyCompare",
  build_vignettes = TRUE,
  dependencies = TRUE
)
```

## Quick start

Please check the
[vignettes](https://ggrlab.github.io/cyCompare/articles/Tutorial.html)
for a more detailed introduction and examples.—

# How was this package created?

``` r
# for VSCode
install.packages("languageserver")
install.packages("devtools")
usethis::create_tidy_package("/home/gugl/clonedgit/ggrlab/cyCompare")
usethis::proj_activate("/home/gugl/clonedgit/ggrlab/cyCompare")
usethis::use_tidy_style(strict = TRUE)
usethis::use_git()
```

`usethis` tells you to envoke further github-related commands. There is
two ways to continue: 1. Create a personal access token (PAT) and use it
to authenticate with github 2. Manually push the package to github

Pushing manually works fine, but some advanced `usethis` commands won’t
work properly, therefore I will continue with the PAT.

``` r
#
usethis::create_github_token() # if done already, use "github token, ggrlab, PAT" password
gitcreds::gitcreds_set() # Then enter the freshly generated token
usethis::use_github(
    organisation = "ggrlab",
    private = TRUE,
    visibility = "private"
)
```

``` r
usethis::use_tidy_github()
# # 2023-10-23: Do NOT use github action checks as it cannot cope with git lfs,
# # therefore always fails. - Unsolved 2024-05-22
usethis::use_github_action("check-standard")
# usethis::use_tidy_github_labels()
usethis::use_pkgdown_github_pages() # Disallowed  by github for private: "GitHub API error (422):Your current plan does not support GitHub Pages for this repository."
```

Additional information:

``` r
usethis::use_author(
    given = "Gunther",
    family = "Glehr",
    email = "gunthergl@gmx.net",
    role = c("aut", "cre"),
    comment = c("ORCID" = "0000-0002-1495-9162")
)
usethis::use_news_md()
lintr::use_lintr(type = "tidyverse")
# Change file .lintr manually to:
# linters: linters_with_defaults(line_length_linter = line_length_linter(120),indentation_linter = indentation_linter(4)) # see vignette("lintr")
# encoding: "UTF-8"
```

precommit is a wonderful tool to check your code before committing it.

``` r
# https://lorenzwalthert.github.io/precommit/articles/precommit.html
# install.packages("precommit")
# install.packages("reticulate")
# bash::$ conda deactivate
# bash::$ pip3 install pre-commit
precommit::install_precommit()
precommit::use_precommit()
## Use pre-commit-config.yaml from restrictedROC
```

Before committing: `pre-commit install --hook-type pre-push`, then
commit.

Used packages:

``` r
usethis::use_data_table()
# usethis::use_package("devtools")
# usethis::use_package("lifecycle")
precommit::snippet_generate("additional-deps-roxygenize")
```
