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

- **dosing events** — bolus, infusion, `addl`, steady state, multiple
  compartments, exactly as the event table describes them (value and eta
  gradient FD-verified on dosed models);
- **delay differential equations**, via rxode2's `delay()` / `past()`;
- **the residual-error and censoring machinery** (M2/M3/M4) already
  validated in nlmixr2est.

One capability is machinery-ready but deliberately **refused** for now:
*estimating* a dose-handling parameter (lag time, bioavailability, duration,
rate).  Its derivative needs a jump condition at the event — rxode2's
event-sensitivity machinery (`eventSens = "jump"`) provides exactly that,
but the linked theta-sensitivity model does not yet request it, so the
column comes back silently zero while the true derivative is not
(nlmixr2est#946, caught by this package's finite-difference gate).
`est="stan"` refuses such models with an explanation rather than sampling a
value/gradient mismatch; `fix()` the parameter to proceed.

The model is compiled and linked **without engaging Stan's threading** —
rxode2 and nlmixr2est parallelize over subjects inside each likelihood
evaluation, and a second threading runtime layered on top would contend with
it.  Stan runs single-threaded; chains run sequentially; the cores go to the
subject loop (`stanControl(likCores=)`).

Two gradient tiers are provided, both supplying analytic derivatives to Stan
via `precomputed_gradients` (nothing is re-derived by autodiff):

1. per-individual conditional gradients `d/d(eta) log p(y_i|eta_i)` through
   nlmixr2est;
2. plus `d/d(theta)` of the conditional at fixed eta (forward sensitivities
   for non-mu structural and residual thetas; the mu-reference identity
   `d/dtheta_p = d/deta_k` for the rest), so Stan samples theta + Omega +
   etas jointly.  Stan owns the full, normalized eta prior (non-centred,
   `eta = z L'`), so its autodiff composes the supplied gradients through the
   Omega parameterization for free.

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
| bitwise determinism and thread-count invariance of every batch evaluation | exact |

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
(add/prop/combined and transforms with fixed lambda), mu-referenced and
non-mu etas, all 39 real-valued prior distributions in the lotri catalogue
(univariate, multivariate normal families, LKJ/Wishart on omega blocks).

Refused with an explanatory error (not silently wrong): estimated
transform-both-sides lambda (the linked conditional omits the DV-transform
Jacobian), mu-referenced covariates (the covariate-coefficient gradient is
not wired yet), mixture models, IOV, and the 8 discrete distributions
(every `ini({})` parameter is real-valued).

This package requires nlmixr2est with the FOCEi conditional-likelihood C API
(nlmixr2est#937, #939, #941).  See the issue tracker for the roadmap
(WAIC/LOO, covariate gradients, SBC).

Supersedes the prior-specification blocker in nlmixr2/nlmixr2est#799.
