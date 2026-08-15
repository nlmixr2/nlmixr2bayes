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
rxode2 and nlmixr2est parallelize over subjects inside each likelihood
evaluation, and a second threading runtime layered on top would contend with
it.  Stan runs single-threaded; chains run sequentially; the cores go to the
subject loop (`stanControl(likCores=)`). The rate limiting step in most
pharamcometric models is the ODE solving; therefore this is the approach that
will likely give the best speed.

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

Prior distributions come from the `ini({})` block (lotri/rxode2), whose
distribution catalogue is deliberately one-to-one with Stan's — see
`lotri::lotriPriorDists()` and `rxode2::rxUiPriors()`.  A model without
priors is refused, printing the exact `prior()` lines to add.

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

Supported: ODE and `linCmt()` models, normal residual models
(add/prop/combined and transforms with fixed lambda), censored data
(M2/M3/M4 via CENS/LIMIT) and user-written `ll()` endpoints (both
gate-verified: values tie to textbook densities up to a parameter-free
constant, gradients FD-verified), mu-referenced and non-mu etas,
mu-referenced covariates (subject-constant), all 39 real-valued prior
distributions in the lotri catalogue (univariate, multivariate normal
families, LKJ/Wishart on omega blocks).

Refused with an explanatory error (not silently wrong): estimated
transform-both-sides lambda (the linked conditional omits the DV-transform
Jacobian), a time-varying covariate on a mu-referenced coefficient (the
`cov_i * d/deta` factorization needs the covariate constant within
subject), mixture models, IOV, and the 8 discrete distributions (every
`ini({})` parameter is real-valued).

This package requires nlmixr2est with the FOCEi conditional-likelihood C API
(nlmixr2est#937, #939, #941).  See the issue tracker for the roadmap
(WAIC/LOO, covariate gradients, SBC).

Supersedes the prior-specification blocker in nlmixr2/nlmixr2est#799.
