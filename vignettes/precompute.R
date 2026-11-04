#!/usr/bin/env Rscript
## Precompute the nlmixr2bayes vignettes (and re-render the README).
##
## The vignette sources are `vignettes/*.Rmd.orig`.  This script knits each one
## to the `vignettes/*.Rmd` that ships in the package, with every chunk's output
## (and any figure, as `vignettes/fig-*.png`) baked in as static markdown.  CRAN
## and `R CMD build` then only convert markdown to HTML -- no Stan program is
## compiled and nothing is sampled at package-build time.  That is deliberate:
## a Stan compile inside a vignette build fails on Windows and runs out of
## memory on macOS CI, and would cost CRAN minutes per platform besides.
##
## The expensive fits inside the documents are run for real and cached in
## `inst/cache` with nlmixr2save's `:=` operator (see cacheSetup.R), so a
## re-knit reloads them instead of re-sampling.  `inst/cache` is
## `.Rbuildignore`d: it only feeds this script and the README.
##
##   Rscript vignettes/precompute.R            # fit only what is missing
##   Rscript vignettes/precompute.R --clean    # refit everything
##
## Re-run this, and commit the regenerated `vignettes/*.Rmd` and figures,
## whenever a `.Rmd.orig` changes or a change alters the numbers it shows.
##
## Each document is knit in its OWN fresh R subprocess: building many rxode2
## models in one long-lived session accumulates loaded model DLLs and
## eventually fails with "error building model", and a Stan compile leaves
## PKG_CPPFLAGS pointing at Stan's headers, so isolating each knit keeps this
## reliable.

if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("precompute.R needs the 'rmarkdown' package")
}

## Run from the package root whichever way it was invoked.
if (basename(getwd()) == "vignettes") {
  setwd("..")
}
if (!file.exists("DESCRIPTION")) {
  stop("run precompute.R from the nlmixr2bayes package root (or vignettes/)")
}

rscript <- file.path(R.home("bin"), "Rscript")

## Vignettes to precompute, by base name (vignettes/<name>.Rmd.orig).
vignettes <- c(
  "nlmixr2bayes",
  "rxstan-handcoded",
  "rxstan-dde",
  "torsten-comparison"
)

args <- commandArgs(trailingOnly = TRUE)
if ("--clean" %in% args) {
  unlink(list.files("inst/cache", pattern = "\\.(zip|rds)$", full.names = TRUE))
  message("precompute.R: cleared inst/cache/")
}
dir.create("inst/cache", showWarnings = FALSE, recursive = TRUE)

runSub <- function(label, cmd, wd = ".") {
  .old <- setwd(wd)
  on.exit(setwd(.old))
  status <- system2(rscript, c("-e", shQuote(cmd)))
  if (!identical(status, 0L)) {
    message("precompute.R: FAILED on ", label, " (exit ", status, ")")
    return(label)
  }
  character(0)
}

failed <- character(0)

## README.md is a committed artifact of README.Rmd, so it renders in place.
message("precompute.R: rendering README.Rmd ...")
failed <- c(
  failed,
  runSub(
    "README.Rmd",
    'rmarkdown::render("README.Rmd", quiet = TRUE, envir = new.env(parent = globalenv()))'
  )
)

for (v in vignettes) {
  .orig <- file.path("vignettes", paste0(v, ".Rmd.orig"))
  if (!file.exists(.orig)) {
    warning("precompute.R: skipping missing vignette source ", .orig)
    next
  }
  message("precompute.R: knitting ", .orig, " ...")
  ## Knit from inside vignettes/, so `source("cacheSetup.R")` resolves and the
  ## figures land next to the .Rmd that references them.  Stale figures for
  ## this vignette are removed first so a dropped plot does not linger.
  unlink(list.files("vignettes", pattern = paste0("^fig-", v, "-.*\\.png$"), full.names = TRUE))
  failed <- c(
    failed,
    runSub(
      .orig,
      sprintf(
        paste0(
          ## knit() defaults to error = TRUE, which would bake a failing
          ## chunk's "Error:" into the shipped vignette and still exit 0;
          ## abort instead, like rmarkdown::render().  Chunks that demonstrate
          ## an error set error = TRUE themselves.  fig.path keeps every
          ## figure at vignettes/fig-<name>-*.png, which the cleanup above and
          ## CLAUDE.md assume.
          'knitr::opts_chunk$set(error = FALSE, fig.path = "fig-%s-"); ',
          'knitr::knit("%s.Rmd.orig", output = "%s.Rmd", quiet = TRUE, envir = new.env(parent = globalenv()))'
        ),
        v,
        v,
        v
      ),
      wd = "vignettes"
    )
  )
}

## A bare `x := rstan::sampling(...)` is saved by nlmixr2save as a plain
## saveRDS() of the whole stanfit -- compiled shared object included (~23 MB),
## because saveFitItem() only runs for items inside an nlmixr2 fit.  Trim those
## the same way saveFitItem.stanfit() does, so the committed cache stays small.
trimCmd <- paste0(
  'pkgload::load_all(quiet = TRUE); ',
  'for (f in list.files("inst/cache", pattern = "\\\\.rds$", full.names = TRUE)) { ',
  '  x <- readRDS(f); ',
  '  if (is.list(x) && methods::is(x$ret, "stanfit")) { ',
  '    x$ret <- .stanTrimStanfit(x$ret); saveRDS(x, f) ',
  '  } ',
  '}'
)
failed <- c(failed, runSub("(stanfit trim)", trimCmd))

message("precompute.R: done. Cached fits:")
print(list.files("inst/cache", pattern = "\\.(zip|rds)$"))
if (length(failed)) {
  message("precompute.R: documents that did NOT build cleanly: ", paste(failed, collapse = ", "))
  quit(status = 1L)
}
