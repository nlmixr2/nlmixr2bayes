# Build-environment probes (issue #1, Spec 3).  Each probe costs milliseconds
# (rstan::stanc only, no compilation) and turns a class of silent future
# breakage -- a stanc3 that changes the external-function call shape, an rstan
# that moves the includes= injection point -- into a load-time refusal with an
# actionable message.

.stanBuildEnv <- new.env(parent = emptyenv())

#' The fixture used by the stanc probes: minimal program declaring the
#' external function exactly as the generator emits it.
#' @noRd
.stanProbeCode <- function() {
  paste(
    "functions {",
    "  vector nlmixr2_cond_all(matrix etaMat);",
    "}",
    "data {",
    "  int<lower=1> N;",
    "  int<lower=1> P;",
    "}",
    "parameters {",
    "  matrix[N, P] eta;",
    "}",
    "model {",
    "  target += sum(nlmixr2_cond_all(eta));",
    "}",
    sep = "\n")
}

#' Probe the Stan build environment
#'
#' Runs `rstan::stanc(allow_undefined=TRUE)` on a fixture and records three
#' facts about the generated C++ (cached per session):
#'
#' 1. `callShape`: stanc3 still emits the external call as
#'    `nlmixr2_cond_all(<expr>, pstream__)` with no template-argument list --
#'    the shape `inst/include/nlmixr2stan_lp.hpp`'s concrete overloads bind to.
#' 2. `injectionPoint`: the first `class <ident>` in the generated code sits
#'    inside `namespace <model>_namespace`, which is where
#'    `rstan::stan_model(includes=)` splices the header.
#' 3. `stancVersion`/`rstanVersion` for the record.
#'
#' @param force re-run the probe even if cached
#' @return a list with `ok`, `callShape`, `injectionPoint`, `stancVersion`,
#'   `rstanVersion`; `ok=NA` when rstan is not installed
#' @export
#' @author Matthew L. Fidler
nlmixr2stanBuildInfo <- function(force = FALSE) {
  if (!force && !is.null(.stanBuildEnv$info)) return(.stanBuildEnv$info)
  if (!requireNamespace("rstan", quietly = TRUE)) {
    .ret <- list(ok = NA, callShape = NA, injectionPoint = NA,
                 stancVersion = NA_character_, rstanVersion = NA_character_)
    .stanBuildEnv$info <- .ret
    return(.ret)
  }
  .sc <- rstan::stanc(model_code = .stanProbeCode(), model_name = "nlmixr2stan_probe",
                      allow_undefined = TRUE)
  .cpp <- .sc$cppcode
  # (1) call shape: the plain two-argument form, no template arguments
  .callShape <- grepl("nlmixr2_cond_all(eta, pstream__)", .cpp, fixed = TRUE) &&
    !grepl("nlmixr2_cond_all<", .cpp, fixed = TRUE)
  # (2) injection point: rstan splices `includes` before the first
  # `class <ident>`; that match must be inside the model namespace so the
  # header's unqualified definitions resolve
  .clsPos <- regexpr("class[[:space:]]+[A-Za-z_][A-Za-z0-9_]*", .cpp)
  .nsPos <- regexpr("namespace[[:space:]]+[A-Za-z_][A-Za-z0-9_]*_namespace", .cpp)
  .injection <- .clsPos > 0 && .nsPos > 0 && .nsPos < .clsPos
  .ret <- list(ok = isTRUE(.callShape) && isTRUE(.injection),
               callShape = .callShape,
               injectionPoint = .injection,
               stancVersion = tryCatch(as.character(rstan::stan_version()),
                                       error = function(e) NA_character_),
               rstanVersion = as.character(utils::packageVersion("rstan")))
  .stanBuildEnv$info <- .ret
  .ret
}

#' Error unless the Stan build environment probes pass
#' @noRd
.stanAssertBuildOk <- function() {
  .i <- nlmixr2stanBuildInfo()
  if (isTRUE(is.na(.i$ok))) {
    stop("rstan is not installed; install rstan to compile the linked Stan model",
         call. = FALSE)
  }
  if (!isTRUE(.i$callShape)) {
    stop("this stanc3 no longer emits the external-function call shape ",
         "nlmixr2stan's header binds to (nlmixr2_cond_all(eta, pstream__)); ",
         "nlmixr2stan needs updating for rstan ", .i$rstanVersion,
         call. = FALSE)
  }
  if (!isTRUE(.i$injectionPoint)) {
    stop("this rstan no longer injects includes= inside the model namespace; ",
         "nlmixr2stan needs updating for rstan ", .i$rstanVersion,
         call. = FALSE)
  }
  invisible(TRUE)
}
