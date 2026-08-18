#' Build a Stan `init` function from a model's `ini()` values
#'
#' Population PK models are routinely multimodal -- a one-compartment oral
#' model can swap absorption and elimination and fit almost as well, the
#' "flip-flop" -- so Stan's default inits, drawn uniformly on
#' `(-2, 2)` unconstrained, regularly send two chains to two different modes
#' and report `Rhat` in the twenties.  Starting from the `ini()` estimates
#' avoids that.
#'
#' Jitter is applied to the fixed effects only.  The between-subject SDs and
#' the residual error are left at their `ini()` values because they are
#' positively constrained and a jitter large enough to be useful is also large
#' enough to push them past zero; the `z` matrix starts at 0, which is the
#' prior mean of the non-centred parameterisation.
#'
#' @param gen the list returned by [rxsStanFromUi()]
#' @param jitter SD of the normal noise added to each fixed effect per chain.
#'   `0` gives every chain the same starting point, which makes `Rhat`
#'   optimistic -- prefer a small positive value.
#' @return a function of `chain_id`, suitable for `rstan::sampling(init = )`
#' @examples
#' \donttest{
#' if (requireNamespace("nlmixr2", quietly = TRUE)) {
#'   # gen <- rxsStanFromUi(model, data)
#'   # fit <- rstan::sampling(sm, data = gen$standata, init = rxsInit(gen))
#' }
#' }
#' @author Lukas A. Widmer
#' @export
rxsInit <- function(gen, jitter = 0.2) {
  if (is.null(gen$inits)) {
    stop("rxsInit() needs the list returned by rxsStanFromUi()", call. = FALSE)
  }
  base <- gen$inits
  thetas <- gen$stanNames$theta
  if (jitter < 0) stop("rxsInit(): jitter must be >= 0", call. = FALSE)

  function(chain_id = 1L) {
    out <- base
    if (jitter > 0 && length(thetas)) {
      for (nm in intersect(thetas, names(out))) {
        out[[nm]] <- out[[nm]] + stats::rnorm(1L, 0, jitter)
      }
    }
    out
  }
}
