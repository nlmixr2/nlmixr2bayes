## Shared fit cache for the nlmixr2bayes documentation.
##
## Every expensive computation in the README and the vignettes is run FOR REAL
## and then cached with nlmixr2save's `:=` operator, so each document renders
## with live output instead of hand-transcribed numbers, and without refitting
## on every build.  `:=` in trusted-cache mode (`nlmixr2save.check = FALSE`)
## loads the committed result when the cache file is there and runs + saves it
## when it is not, so the .Rmd stays the single source of truth for what
## produced the numbers.
##
## The cache lives in `inst/cache`, so it is *installed with the package*: a
## reader who re-renders these documents from an installed nlmixr2bayes gets
## the committed fits rather than an hour of sampling.  It stays small because
## nlmixr2bayes registers a `saveFitItem()` method that drops rstan's compiled
## shared object from the saved `stanfit` (R/stanSave.R) -- roughly 15 MB per
## fit down to 1.5 MB, with every read-only rstan accessor still working.
##
## Usage, from a document's setup chunk:
##
##   source("cacheSetup.R")     # "vignettes/cacheSetup.R" from the package root
##   bayesCache("nlmixr2bayes-")
##
## Repopulate with `Rscript vignettes/precompute.R` (add `--clean` to refit
## everything).

#' Point nlmixr2save's `:=` at the nlmixr2bayes documentation cache
#'
#' @param prefix file-name prefix for this document's cache entries, so two
#'   documents can both cache a variable called `fit` without colliding.
#' @return the cache directory, invisibly
bayesCache <- function(prefix = "") {
  if (!requireNamespace("nlmixr2save", quietly = TRUE)) {
    stop("the nlmixr2bayes documentation caches its fits with nlmixr2save;\n",
         "  install it with pak::pak(\"nlmixr2/nlmixr2save\")", call. = FALSE)
  }
  ## `:=` is an operator, so it has to be on the search path, not merely
  ## installed, for `fit := nlmixr2(...)` to parse into nlmixr2save's method.
  suppressPackageStartupMessages(library("nlmixr2save", character.only = TRUE))
  ## Prefer the source tree, so `precompute.R` and a pkgdown build write into
  ## the copy that is under version control; fall back to the installed cache.
  .src <- Filter(dir.exists, c("inst/cache", "../inst/cache", "../../inst/cache"))
  .src <- if (length(.src)) {
    normalizePath(.src[[1]])
  } else {
    system.file("cache", package = "nlmixr2bayes")
  }
  if (!nzchar(.src)) {
    stop("cannot find the nlmixr2bayes fit cache (inst/cache)", call. = FALSE)
  }
  ## Loading a cached fit unzips it in place, so the directory has to be
  ## writable.  An installed package under a system-wide library is not, so
  ## work from a scratch copy there.
  .dir <- if (file.access(.src, 2L) == 0L) {
    .src
  } else {
    .tmp <- file.path(tempdir(), "nlmixr2bayesCache")
    dir.create(.tmp, showWarnings = FALSE, recursive = TRUE)
    file.copy(list.files(.src, full.names = TRUE), .tmp, overwrite = TRUE)
    .tmp
  }
  options(nlmixr2save.dir = .dir,
          nlmixr2save.prefix = prefix,
          ## trusted cache: load what is there, run and save what is not
          nlmixr2save.check = FALSE,
          ## a committed cache necessarily outlives the nlmixr2est/rxode2 build
          ## that produced it; that is the point, so do not warn about it
          nlmixr2save.checkVersion = FALSE)
  invisible(.dir)
}
