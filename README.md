# nlmixr2stan

A Stan interface for nlmixr2: `nlmixr2(model, data, est = "stan")`.

The idea is to **link** Stan to the likelihoods that already exist rather than
re-express the model in Stan's language: rxode2 compiles the ODE/solving code
and nlmixr2est holds the likelihood machinery, so Stan calls into those
directly.  That keeps one source of truth for the likelihood and avoids
reimplementing solving, censoring and the residual-error models.

## What linking buys

Everything rxode2/nlmixr2est can express flows through `est="stan"` without a
translation step, including things Stan's own ODE interface cannot express:

- **dosing events** -- bolus, infusion, `addl`, steady state, multiple
  compartments, exactly as the event table describes them (value and eta
  gradient FD-verified on dosed models);
- **derivatives with respect to modelled dose handling** — an *estimated*
  lag time, bioavailability, duration or rate, whose dependence runs
  through the event time rather than the ODE right-hand side and so needs
  a **jump condition** at the event.  rxode2's analytic event
  sensitivities (`eventSens = "jump"`) supply it, and nlmixr2est's
  theta-sensitivity model carries it through the linked gradient
  (nlmixr2est#946; the `alag` theta column FD-agrees at ~1e-5 where it was
  once silently zero -- this package's finite-difference gate caught it and
  now locks it in);
- **delay differential equations**, via rxode2's `delay()` / `past()`;
- **the residual-error and censoring machinery** (M2/M3/M4) already
  validated in nlmixr2est.

The model is compiled and linked **without engaging Stan's threading** —
a second threading runtime layered on top of rxode2/nlmixr2est's would
contend with it.  Parallelism has two axes instead: by default (unix)
the **chains run in forked processes** (`stanControl(chainCores=)`,
inheriting the `rxode2::getRxThreads()` budget capped at the chain
count; each fork gets a copy-on-write duplicate of the linked state,
and draws are bit-identical to a sequential run), and with
`chainCores=1` the chains run sequentially with **subject-parallel
OpenMP inside each likelihood evaluation** (`stanControl(cores=)`).
The rate-limiting step in most pharmacometric models is the ODE
solving; these are the two places it parallelizes.

Two gradient tiers are provided, both supplying analytic derivatives to Stan
via `precomputed_gradients` (nothing is re-derived by autodiff):

1. per-individual conditional gradients `d/d(eta) log p(y_i|eta_i)` through
   nlmixr2est;
2. plus `d/d(theta)` of the conditional at fixed eta (forward sensitivities
   for non-mu structural and residual thetas; the mu-reference identity
   `d/dtheta_p = d/deta_k` for the rest, extended to mu-referenced
   covariate coefficients as `d/dtheta_p = cov_i * d/deta_k` when the
   sensitivity model does not already carry them), so Stan samples theta +
   Omega + etas jointly.  Stan owns the full, normalized eta prior
   (non-centred, `eta = z L'`), so its autodiff composes the supplied
   gradients through the Omega parameterization for free.

The bridge is **first-order only**: `precomputed_gradients()` builds a
reverse-mode `var`, and Stan Math has no forward-mode (`fvar`) counterpart
for externally supplied partials, so anything needing higher-order autodiff
of the target through the external term — Riemannian HMC's metric, the
embedded-Laplace (`laplace_marginal`) machinery — is out of scope.  That
costs nothing in practice: NUTS, ADVI, optimization and Pathfinder are all
first-order, and a Laplace-*marginalized* target is exactly what
nlmixr2est's own (NONMEM-validated) FOCEi/Laplace machinery computes
natively — linking that marginal would be the sensible route, not pushing
second-order AD through the bridge.

A planned companion (merging an existing standalone rxode2--Stan
bridge) links at the **ODE level** instead: hand-written Stan programs
that call rxode2's compiled solver and analytic sensitivities directly
-- your own Stan code, rxode2's events/DDEs/jump sensitivities.  This
package links at the **likelihood level**: no Stan code is written at
all.  Same architecture, two entry points.

Beyond NUTS, `stanControl(algorithm=)` runs Stan's ADVI variational
approximations (`"meanfield"`, `"fullrank"`) and multi-path Pathfinder
(implemented in-package against the model's exact log density and
analytic gradient, since rstan does not expose the Pathfinder service)
on the same generated program and linked likelihood, returning the same
complete nlmixr2 fit 3-10x faster; the Pareto-k̂ diagnostic replaces
Rhat/ESS and warns when the approximation is not trustworthy (use NUTS
for the final answer).  `est="advi"` and `est="pathfinder"` are sugar
for these, the way `"foce"`/`"focei"` are members of the focei family.

Prior distributions come from the `ini({})` block (lotri/rxode2), whose
distribution catalogue is deliberately one-to-one with Stan's — see
`lotri::lotriPriorDists()` and `rxode2::rxUiPriors()`.  A model without
priors is refused, printing the exact `prior()` lines to add.

## Is it fast?

Benchmarked against the *same* 1-compartment oral ODE population model
(theo_sd, 12 subjects) written natively in Stan with `ode_rk45_tol` --
identical priors, identical non-centred eta structure, matched inits,
matched tolerances (1e-8), matched sampler settings (Stan defaults),
and the same 4-core budget spent the way each method naturally can
(native: parallel chains; nlmixr2stan: subject-parallel threads inside
each gradient, chains sequential).  The linked solver is rxode2's dense
Dormand-Prince (`dop853`, `dense=TRUE`), the twin of Stan's `ode_rk45`:

|  | native Stan | nlmixr2stan (linked) |
|---|---|---|
| wall, 2 chains x 1000 iter (chains sequential, 4 subject threads) | — | 408.8 s |
| wall, 2 chains x 1000 iter (chains FORKED, serial inner) | 142.7 s | 302.6 s |
| worst bulk ESS | 236 | **286-378** |
| worst ESS / s | 1.66 | 0.70 sequential / **1.25 forked** |
| gradient evaluation | 1.98 ms | 2.80 ms via `grad_log_prob` (1.36 ms at the C level) |
| posterior means | agree to < 0.01 on every parameter | |

The posteriors are the same; the linked run extracts *more* effective
samples per draw; forked chains (now the default) close most of
native's wall-clock edge, and disabling the failure cascade
(`maxOdeRecalc=0, fallbackFD=FALSE`) changes nothing measurable -- the
cascade never fires on a healthy trajectory (the forked retry-on run
reproduced the sequential retry-off draws bit-for-bit).  The former
two-integrations-per-gradient gap is closed: with nlmixr2est's combined
eta+theta sensitivity build (nlmixr2est#958) the whole tier-2 gradient
comes from ONE solve per subject through a fused batch entry --
0.62 ms per full-population gradient at the C level (vs 1.26 ms
two-model, and native's 1.98 ms coupled solve).  nlmixr2stan negotiates
it automatically when the loaded nlmixr2est provides it.  Per gradient the linked C path is
cheaper than the native solve -- and for models with closed-form
solutions the `linCmt()` tier evaluates ~5-10x cheaper still, an
option a native ODE implementation does not have.  ADVI completes the
same fit in ~30-40 s and Pathfinder in ~90 s (with the Pareto-k-hat
diagnostic guarding both).  Reproduce with
`Rscript inst/bench/native-vs-linked.R`.

## Is it right?

Because the bridge operates at the **likelihood level** (it exposes
`log p(y_i|eta_i)` and its gradients, not raw ODE solutions), the oracles are
likelihood-level.  Independent oracles first — checks against something the
bridge does not share code with:

| check | agreement |
|---|---|
| Stan's assembled gradient (`grad_log_prob`) vs. Richardson finite differences of its own log density | ~3e-8 |
| `d/d(eta)` and `d/d(theta)` vs. central differences of the linked conditional value | ~3e-8 |
| quadrature marginal over the linked conditional vs. nlmixr2's own `est="agq"` (nAGQ=101) at matched theta | 3e-7 (independent *marginalizations* — trapezoid vs adaptive Gauss–Hermite; the integrand engine is shared, and is anchored absolutely by the hand-density row below) |
| censored (M2/M3/M4) and `ll()` conditionals vs. hand-computed textbook densities | exact up to a pinned, parameter-free constant; all gradients FD-verified, including the residual-SD dependence of the censored CDF terms |
| prior-only sampling (likelihood stripped) vs. every declared prior | one-sample KS per parameter, including the default LKJ + half-Cauchy omega path |
| full posterior vs. a hand-written **native-Stan** implementation of the same model and priors (no external function) | within 3× combined MCSE on every parameter |
| the linked likelihood engine itself | NONMEM-validated upstream in nlmixr2est (Wang 2007 objective, per-subject ETA/CWRES) |

And internal-consistency checks — exact identities the implementation must
satisfy:

| check | agreement |
|---|---|
| linked value vs. `nlmixr2est::foceiLikRun(type="cond")` | exact (same engine, one code path) |
| mu-referenced theta-gradient column vs. the eta gradient | exact |
| assembled gradient under Omega^-1 swaps spanning 4 orders of magnitude | exact (the conditional is Omega-free by construction) |
| `log_prob` decomposes into conditional + normalized eta prior | ~1e-6 (difference form, constants cancel) |
| bitwise determinism: 500+ interleaved evaluations of both batch entries at alternating points | exact (bitwise identical — the target is a pure function of state, which NUTS's reversibility assumes) |

Plus simulate-and-recover fits, including a correlated omega block whose
generative correlation lands inside the 90% credible interval, and
FOCEi-agreement runs on `nlmixr2data::theo_sd` under weak priors.

Wrong-by-construction variants are asserted to **fail**: the historical
sign-flipped gradient assembly (off by `2 Omega^-1 eta`) fails the FD test,
and `likelihood="focep"` — whose value and gradient are gradients of
different functions — is refused by capability flags, with the underlying FD
mismatch itself locked in as a test upstream.

## Current scope

Supported: mixed AND population-only models (a no-eta model -- e.g. a
single-subject Bayesian fit -- runs as "tier 0" through nlmixr2est's nlm
path: one external scalar `nlmixr2_pop_ll(theta)` carrying the complete
data log-likelihood and its analytic gradient, value tied to a
hand-written density and FD-gate-verified; needs nlmixr2est with
nlmixr2est#953), ODE and `linCmt()` models, normal residual models
(add/prop/combined and transforms with fixed lambda), censored data
(M2/M3/M4 via CENS/LIMIT) and user-written `ll()` endpoints (both
gate-verified: values tie to textbook densities up to a parameter-free
constant, gradients FD-verified), mu-referenced and non-mu etas,
mu-referenced covariates (subject-constant and time-varying), all 39
real-valued prior distributions in the lotri catalogue (univariate,
multivariate normal families, LKJ/Wishart on omega blocks).

Mu-referenced covariates work for both subject-constant covariates (the
`cov_i * d/deta` scatter identity) and time-varying covariates (the
coefficient rides the forward-sensitivity model like any other structural
theta, the same way the other nlmixr2est methods treat a time-varying
regressor); both routes are FD-verified.

Estimated transform-both-sides lambda (Box-Cox / Yeo-Johnson) is
supported: the linked conditional supplies the transformed-scale density
with an exact d/dlambda column (nlmixr2est#949; FD-verified at ~1e-10),
and the generator adds the DV-transform Jacobian Stan-side as
`target += (lambda - 1) * sumLogJac` — a pure data statistic, so lambda's
full gradient is exact end to end (the assembled target is gate-verified
against the untransformed-scale density in difference form).

Finite mixtures (`mix()`, 2 components) are supported: the linked batch
evaluates the component-conditional likelihoods in nlmixr2est's blessed
component-major layout (nlmixr2est#955), the generator marginalizes with
the Stan Users Guide `log_sum_exp` pattern (component-specific etas whose
priors factor out; the mixing probability's gradient is pure Stan
autodiff -- the conditional is provably p-free, FD-verified), and the
per-subject membership posteriors land in `fit$env$mixProb`.  The whole
assembled mixture target is FD-verified through `grad_log_prob`.

Refused with an explanatory error (not silently wrong): mixtures with
more than 2 components (for now), IOV (until nlmixr2est#952), and the 8
discrete distributions (every `ini({})` parameter is real-valued).

This package requires nlmixr2est with the FOCEi conditional-likelihood C API
(nlmixr2est#937, #939, #941).  See the issue tracker for the roadmap
(WAIC/LOO, covariate gradients, SBC).

Supersedes the prior-specification blocker in nlmixr2/nlmixr2est#799.
