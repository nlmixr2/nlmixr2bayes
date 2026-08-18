# These shims are added for completion
nmObjGet.stanfit <- function(x, ...) {
  .obj <- x[[1]]
  .env <- .obj$env
  if (exists("stanfit", envir = .env, inherits = FALSE)) {
    get("stanfit", envir = .env, inherits = FALSE)
  } else {
    NULL
  }
}
attr(nmObjGet.stanfit, "desc") <- "Stan fit object"
attr(nmObjGet.stanfit, "rstudio") <- list(0)

nmObjGet.posteriorSummary <- function(x, ...) {
  .obj <- x[[1]]
  .env <- .obj$env
  if (exists("posteriorSummary", envir = .env, inherits = FALSE)) {
    get("posteriorSummary", envir = .env, inherits = FALSE)
  } else {
    NULL
  }
}
attr(nmObjGet.stanfit, "desc") <- "Stan posterior summary"
