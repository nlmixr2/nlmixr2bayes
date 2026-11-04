# Register an rxode2 model as a Stan-callable solver

Compiles \`model\` with analytic first-order sensitivities with respect
to \`sens\`, and returns an integer handle that a Stan program can pass
to \`rx_solve()\`.

## Usage

``` r
rxsRegister(
  model,
  events,
  sens,
  output,
  params = numeric(0),
  atol = 1e-08,
  rtol = 1e-08,
  method = "liblsoda",
  eventSens = NULL,
  perSubject = FALSE,
  quiet = TRUE,
  fast = TRUE,
  ...
)
```

## Arguments

- model:

  rxode2 model source: a string, a \`rxode2\` expression, or an
  already-compiled model. It is recompiled with \`calcSens = sens\`.

- events:

  event table (see \[rxode2::et()\]) or a data frame of events.

- sens:

  character vector of parameter names Stan will supply and differentiate
  with respect to, in the exact order Stan will pass them.

- output:

  character vector of state names to return, in order.

- params:

  named numeric vector of any remaining fixed parameters.

- atol, rtol:

  solver tolerances. Defaults are deliberately tight: HMC needs the
  solver error well below the gradient signal or the trajectories stop
  being reversible.

- method:

  rxode2 integration method.

- eventSens:

  \`"jump"\` for analytic dosing-event sensitivities (modeled lag time,
  bioavailability, duration, rate), \`"fd"\` for the legacy
  finite-difference behavior, or \`NULL\` to leave rxode2's default
  alone.

- perSubject:

  give every subject its own block of \`sens\` parameters, as a
  population model needs. Stan then supplies \`nsub \* length(sens)\`
  values, block-major (subject 1's parameters first). Because subjects
  are independent the returned Jacobian stays \`ny x length(sens)\`:
  each output only depends on its own subject's block.

- quiet:

  suppress rxode2's per-solve diagnostics. HMC proposes bad parameters
  routinely and each rejected proposal otherwise prints a block of
  \`lsoda\`/\`intdy\` warnings, which buries everything else over a real
  fit. Only that repeated chatter is silenced: a failed solve is still
  turned into an error that Stan sees and rejects on, and
  \[rxsSolveStats()\] keeps counting failures either way, so nothing is
  lost.

- fast:

  use Path A, which re-drives the solve rxode2 has already built instead
  of calling \[rxode2::rxSolve()\] again. That skips roughly 2.3 ms of
  per-call R setup and teardown, which dominates small-to-medium
  problems. Any other \`rxSolve()\` in the session invalidates it; the
  bridge detects that and falls back to the slow path for one solve,
  which re-arms it. It is switched off automatically when rxode2 reports
  a non-zero \`needSort\` for the sensitivity-augmented model \*and\*
  the event table contains dose rows – a modeled \`alag()\`, \`f()\`,
  \`dur()\` or \`rate()\`, whose event ordering Path A cannot refresh.
  Both conditions are required because all four quantities act on doses,
  so an observation-only event table cannot be re-ordered; rxode2 sets
  the flag for a delay differential equation, which has no doses at all.
  \[print()\] on the handle reports when the fast path was refused.

- ...:

  passed through to \[rxode2::rxSolve()\] on every solve.

## Value

an integer handle, with attributes describing the layout

## Details

The returned outputs are ODE \*states\*. rxode2 emits sensitivities only
for states (\`rx\_\_sens\_\<state\>\_BY\_\<param\>\_\_\`), never for
derived \`lhs\` quantities, so any derived prediction (a concentration
\`center / v\`, say) must be formed in the Stan program, where Stan
differentiates it itself.

## Author

Lukas A. Widmer
