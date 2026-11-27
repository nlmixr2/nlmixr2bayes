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
#' prior mean of the non-centered parameterization.
#'
#' @param gen the list returned by [rxsStanFromUi()]
#' @param jitter SD of the normal noise added to each fixed effect per chain.
#'   `0` gives every chain the same starting point, which makes `Rhat`
#'   optimistic -- prefer a small positive value.
#' @return a function of `chain_id`, suitable for `rstan::sampling(init = )`
#' @examples
#' \donttest{
#' mod <- function() {
#'   ini({
#'     tka <- 0.5
#'     tcl <- 1
#'     tv <- 3.4
#'     eta.cl ~ 0.1
#'     add.sd <- 0.2
#'   })
#'   model({
#'     ka <- exp(tka)
#'     cl <- exp(tcl + eta.cl)
#'     v <- exp(tv)
#'     d/dt(depot) <- -ka * depot
#'     d/dt(center) <- ka * depot - cl / v * center
#'     cp <- center / v
#'     cp ~ add(add.sd)
#'   })
#' }
#' ev <- rxode2::et(amt = 100, cmt = "depot")
#' ev <- as.data.frame(rxode2::et(ev, c(0.5, 1, 2, 4, 8, 12, 24)))
#' d <- rbind(cbind(ev, id = 1L), cbind(ev, id = 2L))
#' d$dv <- 1
#' gen <- rxsStanFromUi(mod, d)
#' # one init list per chain; the fixed effects are jittered
#' init <- rxsInit(gen)
#' str(init(chain_id = 1))
#' # jitter = 0 starts every chain at the ini() values
#' str(rxsInit(gen, jitter = 0)(chain_id = 2))
#' }
#' @author Lukas A. Widmer
#' @export
rxsInit <- function(gen, jitter = 0.2) {
  if (is.null(gen$inits)) {
    stop("rxsInit() needs the list returned by rxsStanFromUi()", call. = FALSE)
  }
  base <- gen$inits
  thetas <- gen$stanNames$theta
  if (jitter < 0) {
    stop("rxsInit(): jitter must be >= 0", call. = FALSE)
  }

  # rstan calls the init function with this exact argument name.
  # nolint start: object_name_linter.
  function(chain_id = 1L) {
    # nolint end
    out <- base
    if (jitter > 0 && length(thetas)) {
      for (nm in intersect(thetas, names(out))) {
        out[[nm]] <- out[[nm]] + stats::rnorm(1L, 0, jitter)
      }
    }
    out
  }
}
