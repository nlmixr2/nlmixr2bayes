# nlmixr2stan

A Stan interface for nlmixr2.

The idea is to **link** Stan to the likelihoods that already exist rather than
re-express the model in Stan's language: rxode2 compiles the ODE/solving code
and nlmixr2est holds the likelihood machinery, so Stan should call into those
directly. That keeps one source of truth for the likelihood and avoids
reimplementing solving, censoring and the residual-error models.

The model is compiled and linked **without Stan's threading (TBB) attached** —
rxode2 and nlmixr2est do their own parallelism over subjects, and a second
threading runtime layered on top would contend with it.

Three tiers of gradient provision are planned:

1. models providing only the outer gradient through nlmixr2est
2. models providing per-individual inner gradients through nlmixr2est
3. models providing inner gradients through nlmixr2 plus the outer likelihood
   gradient through nlmixr2est

Prior distributions come from the `ini({})` block (lotri/rxode2), whose
distribution catalogue is deliberately one-to-one with Stan's — see
`lotri::lotriPriorDists()` and `rxode2::rxUiPriors()`.

This package is in the design stage; see the issue tracker.

Supersedes the prior-specification blocker in nlmixr2/nlmixr2est#799.
