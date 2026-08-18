#' Path to the rxstan C++ headers
#' @author Lukas A. Widmer
#' @export
rxsIncludePath <- function() {
  system.file("include", package = "nlmixr2bayes", mustWork = TRUE)
}

#' Build the `includes=` string for [rstan::stan_model()]
#'
#' rstan splices `includes` immediately before the model class, which is
#' *inside* `model_<name>_namespace`.  The bridge needs global scope (it defines
#' a namespace, template specializations, and static function-pointer caches),
#' so the returned text closes that namespace, emits the bridge at global scope,
#' and reopens it with a `using` declaration that makes `rx_solve` visible to
#' the generated unqualified call.
#'
#' Ordering matters and is why `pre` and `post` are separate: a custom tag type
#' must be declared *before* the header (the primary template is named in
#' `rx_solve`'s body, which is bound at definition time), while a
#' `solve_policy` specialization for it can only be written *after* the primary
#' template exists.
#'
#' @param modelName model name passed to [rstan::stanc()]; must have been
#'   compiled with `obfuscate_model_name = FALSE` so the namespace is
#'   predictable.
#' @param tag optional C++ tag type selecting a `rxstan::solve_policy`
#'   specialization.  Defaults to rxstan's analytic policy; pass
#'   `"::rxstan::finite_diff_tag"` to compile the same program against central
#'   differences instead.
#' @param pre C++ emitted at global scope *before* the bridge header, e.g. the
#'   declaration of a custom tag type.
#' @param post C++ emitted at global scope *after* the bridge header, e.g. a
#'   generated `rxstan::solve_policy` specialization.
#' @return a single string
#' @author Lukas A. Widmer
#' @export
rxsStanIncludes <- function(modelName, tag = NULL, pre = NULL, post = NULL) {
  stopifnot(is.character(modelName), length(modelName) == 1L)
  header <- file.path(rxsIncludePath(), "rxstan", "rxstan.hpp")

  paste0(
    "\n} // close the model namespace so the bridge lands at global scope\n",
    if (!is.null(pre)) paste0(pre, "\n") else "",
    if (!is.null(tag)) paste0("#define RXSTAN_MODEL_TAG ", tag, "\n") else "",
    "#include \"", header, "\"\n",
    if (!is.null(post)) paste0(post, "\n") else "",
    "namespace model_", modelName, "_namespace {\n",
    "using ::rxstan::rx_solve;\n")
}
