# The NIMBLE shim: a tiny standalone C entry point (inst/nimble/*.c) that
# resolves nlmixr2bayes_cond_batch_theta via R_GetCCallable at runtime,
# exactly like inst/include/nlmixr2bayes_lp.hpp does for Stan.  Compiled on
# demand -- nimbleExternalCall(oFile=) is a strict make-target dependency
# that wants a plain .o, never a linked .so, and linking directly against a
# package build's own .o would reach a disconnected copy of its static state
# (see the header comment in nlmixr2bayes_nimble_shim.c) -- and cached by
# content hash, mirroring stanCompile()'s cache (R/stanCompile.R).

#' Path to the shipped shim source
#' @noRd
.nimbleShimSource <- function() {
  .c <- system.file("nimble", "nlmixr2bayes_nimble_shim.c", package = "nlmixr2bayes")
  if (!nzchar(.c)) stop("nlmixr2bayes_nimble_shim.c not found", call. = FALSE) # nocov
  .c
}

#' Path to the shipped shim header (what `nimbleExternalCall(headerFile=)` needs)
#' @noRd
.nimbleShimHeader <- function() {
  .h <- system.file("nimble", "nlmixr2bayes_nimble_shim.h", package = "nlmixr2bayes")
  if (!nzchar(.h)) stop("nlmixr2bayes_nimble_shim.h not found", call. = FALSE) # nocov
  normalizePath(.h, winslash = "/")
}

#' Compile the NIMBLE shim to an object file, cached by content hash
#'
#' @param cache use (and populate) the on-disk cache
#' @param cacheDir cache directory (default
#'   `tools::R_user_dir("nlmixr2bayes", "cache")`)
#' @return path to the compiled `.o`
#' @noRd
.nimbleShimCompile <- function(cache = TRUE, cacheDir = NULL) {
  rxode2::rxReq("nimble")
  .src <- .nimbleShimSource()
  .key <- digest::digest(list(readLines(.src),
                              as.character(utils::packageVersion("nimble")),
                              as.character(utils::packageVersion("nlmixr2bayes")),
                              R.version.string))
  .dir <- if (is.null(cacheDir)) tools::R_user_dir("nlmixr2bayes", "cache") else cacheDir
  .o <- file.path(.dir, paste0("nimbleshim-", .key, ".o"))
  if (cache && file.exists(.o)) return(.o)
  dir.create(.dir, recursive = TRUE, showWarnings = FALSE)
  .tmp <- tempfile("nlmixr2bayes_nimble_shim_", fileext = ".c")
  file.copy(.src, .tmp, overwrite = TRUE)
  .wd <- dirname(.tmp)
  .base <- sub("\\.c$", "", basename(.tmp))
  .old <- setwd(.wd)
  on.exit(setwd(.old), add = TRUE)
  .rbin <- file.path(R.home("bin"),
                     if (identical(.Platform$OS.type, "windows")) "R.exe" else "R")
  .out <- system2(.rbin, c("CMD", "SHLIB", paste0(.base, ".c")),
                  stdout = TRUE, stderr = TRUE)
  .oGen <- file.path(.wd, paste0(.base, ".o"))
  if (!file.exists(.oGen)) {
    stop("failed to compile the nimble shim:\n", paste(.out, collapse = "\n"),
         call. = FALSE)
  }
  if (cache) {
    file.copy(.oGen, .o, overwrite = TRUE)
    return(.o)
  }
  .oGen
}
