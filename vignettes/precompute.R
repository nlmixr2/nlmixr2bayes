#!/usr/bin/env Rscript
## Populate inst/cache/ with the fits the nlmixr2bayes documentation shows.
##
## Every fit in the README and the vignettes is run for real and cached with
## nlmixr2save's `:=` operator (see cacheSetup.R).  This script renders each
## document to a throwaway location; the `:=` calls inside the documents do the
## actual fitting and write the results into inst/cache.  Because the render
## runs the real document code there is a single source of truth -- the .Rmd --
## and no way for the cache to drift from what the page claims produced it.
##
##   Rscript vignettes/precompute.R            # fit only what is missing
##   Rscript vignettes/precompute.R --clean    # refit everything
##
## Each document is rendered in its OWN fresh R subprocess: building many
## rxode2 models in one long-lived session accumulates loaded model DLLs and
## eventually fails with "error building model", and a Stan compile leaves
## PKG_CPPFLAGS pointing at Stan's headers, so isolating each render keeps this
## reliable.

if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("precompute.R needs the 'rmarkdown' package")
}

## Run from the package root whichever way it was invoked.
if (basename(getwd()) == "vignettes") setwd("..")
if (!file.exists("DESCRIPTION")) {
  stop("run precompute.R from the nlmixr2bayes package root (or vignettes/)")
}

rscript <- file.path(R.home("bin"), "Rscript")

## Documents whose fits are cached.  Add one here once its expensive calls use
## the `:=` operator.
docs <- c("README.Rmd",
          "vignettes/nlmixr2bayes.Rmd",
          "vignettes/rxstan-handcoded.Rmd",
          "vignettes/rxstan-dde.Rmd",
          "vignettes/torsten-comparison.Rmd")

args <- commandArgs(trailingOnly = TRUE)
if ("--clean" %in% args) {
  unlink(list.files("inst/cache", pattern = "\\.(zip|rds)$", full.names = TRUE))
  message("precompute.R: cleared inst/cache/")
}
dir.create("inst/cache", showWarnings = FALSE, recursive = TRUE)

outDir <- tempfile("nlmixr2bayes-precompute-")
dir.create(outDir)

failed <- character(0)
for (d in docs) {
  if (!file.exists(d)) {
    warning("precompute.R: skipping missing document ", d)
    next
  }
  message("precompute.R: rendering ", d, " to populate the cache ...")
  ## README.md is a committed artifact of README.Rmd, so that one renders in
  ## place; the vignettes only need to run, so they go to a throwaway dir.
  cmd <- if (identical(d, "README.Rmd")) {
    sprintf('rmarkdown::render("%s", quiet=TRUE, envir=new.env(parent=globalenv()))', d)
  } else {
    sprintf('rmarkdown::render("%s", output_dir="%s", quiet=TRUE, envir=new.env(parent=globalenv()))',
            d, outDir)
  }
  status <- system2(rscript, c("-e", shQuote(cmd)))
  if (!identical(status, 0L)) {
    message("precompute.R: FAILED on ", d, " (exit ", status, ")")
    failed <- c(failed, d)
  }
}

message("precompute.R: done. Cached fits:")
print(list.files("inst/cache", pattern = "\\.(zip|rds)$"))
if (length(failed)) {
  message("precompute.R: documents that did NOT render cleanly: ",
          paste(failed, collapse = ", "))
  quit(status = 1L)
}
