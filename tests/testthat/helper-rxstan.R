## Shared fixtures and a compiled-model cache.  testthat sources helper*.R
## before any test file, and runs every file in one R session, so a plain
## memo here deduplicates compilation across files as well as within them.

## `nlmixr2 = FALSE` for the files that fit hand-written .stan programs and so
## do not need the codegen path at all.
skipUnlessStan <- function(nlmixr2 = TRUE) {
  skip_if_not_installed("rstan")
  if (nlmixr2) skip_if_not_installed("nlmixr2")
  if (!nzchar(Sys.getenv("RXSTAN_STAN_TESTS"))) {
    skip("set RXSTAN_STAN_TESTS=1 to run the Stan compilation tests")
  }
}

.stanCache <- new.env(parent = emptyenv())

## Compiling dominates this suite (644s of a 668s run), and the same program
## is compiled repeatedly -- inst/stan/pk_1cmt_oral.stan in three files, and
## one generated program three times in test-loo.R.
##
## Keyed on the program text AND the bridge header, because the compiled model
## links that header: a key that ignored it would let a header change go
## unnoticed, which is the one failure only this suite catches.  `tag` and any
## other rxsStanModel() argument are in the key too -- the same .stan file
## compiled against finite_diff_tag is a DIFFERENT model, and handing back the
## analytic one would quietly turn that test into a tautology.
##
## The cached model may have been compiled under a different `modelName` than
## the caller asked for; that only changes the generated C++ namespace.
stanModelFor <- function(src, modelName, ...) {
  isPath <- length(src) == 1L && !grepl("\n", src, fixed = TRUE) &&
    file.exists(src)
  code <- if (isPath) readLines(src) else src
  hdr <- system.file("include", "rxstan", "rxstan.hpp", package = "nlmixr2bayes")
  key <- paste(.textHash(code), .textHash(readLines(hdr)),
               as.character(utils::packageVersion("rstan")),
               paste(vapply(list(...), function(a) paste(deparse(a), collapse = ""),
                            character(1)), collapse = "|"),
               sep = "-")

  hit <- .stanCache[[key]]
  if (!is.null(hit)) return(hit)

  ## In-session only, deliberately.  Two attempts at persisting compiled
  ## models across runs both failed silently: rstan::is_sm_valid is not
  ## exported, so a tryCatch around it turned "no such function" into "cache
  ## miss" and recompiled every run; and rstan_options(auto_write) writes
  ## nothing here because rxsStanModel() hands stan_model() the code rather
  ## than a file path.  Both looked like working caches -- the timings were
  ## the only thing that gave them away.  Worth revisiting only with a
  ## fresh-process timing check as the acceptance test.
  f <- tempfile(fileext = ".stan")
  writeLines(code, f)
  sm <- rxsStanModel(f, modelName = modelName, ...)
  assign(key, sm, envir = .stanCache)
  sm
}


## tools::md5sum only works on files, and adding a hashing dependency for the
## test suite is not worth it.
.textHash <- function(lines) {
  f <- tempfile()
  on.exit(unlink(f))
  writeLines(lines, f)
  unname(tools::md5sum(f))
}

## --- model fixtures ----------------------------------------------------

pkModel <- "
ka  <- exp(lka)
cl  <- exp(lcl)
v   <- exp(lv)
d/dt(depot)  <- -ka * depot
d/dt(center) <-  ka * depot - (cl / v) * center
cp <- center / v
"

ddeModel <- "
kin  <- exp(lkin)
kout <- exp(lkout)
tau  <- exp(ltau)
R(0) <- 1
d/dt(R) <- kin - kout * delay(R, tau)
"
