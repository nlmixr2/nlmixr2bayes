# nlmixr2bayes 0.0.1

* First release.

* Adds a Stan-based Bayesian estimation method to 'nlmixr2' through two
  independent subsystems.

* **Likelihood-level linking** (`est = "stan"`, and the algorithm-named
  `"nuts"`, `"advi"` and `"pathfinder"` variants): Stan is linked directly
  to the compiled 'rxode2'/'nlmixr2est' likelihood instead of the model
  being re-expressed in the Stan language.  ODE solving, the
  residual-error models and censoring stay in 'rxode2'/'nlmixr2est'; Stan
  receives the log-likelihood and analytic gradients through
  `precomputed_gradients()` and samples theta, Omega and eta jointly.
  Prior distributions are read from the model's `ini({})` block.

* **Hand-coded Stan models** (the "rxstan" bridge, `rxsRegister()` /
  `rxsStanModel()` / `rxsStanFromUi()`): write your own `.stan` program
  and declare 'rxode2' as the ODE solver backend via `rx_solve()`, for
  models needing Stan code the automatic path does not generate.

* Both paths support dosing events (bolus, infusion, `addl`, steady
  state), derivatives with respect to modeled lag time, bioavailability,
  duration and rate, delay differential equations via `delay()` /
  `past()`, and analytic parameter sensitivities through forward
  sensitivity equations rather than autodiff through the solver.

* Fits carry `$stanfit`, `$stanCode` and `$posteriorSummary`, and WAIC /
  LOO are available through the 'loo' package.

* Any prior Stan samples can be written in `ini({})` -- `dnorm()`,
  `dcauchy()`, `dbeta()`, `dgamma()`, `dlnorm()`, `dexp()`, `dunif()`,
  `dweibull()`, `dlogis()`, `studentT()`, ... -- and all of them are part of
  the posterior.  With `ofv = "focei"` (the default) the fit also reports a
  FOCEi objective row at the posterior point estimate; that row is a plug-in
  likelihood criterion, comparable across nlmixr2 methods, and is evaluated
  *without* the `ini({})` priors.  A fit with priors says so in `$runInfo`.
