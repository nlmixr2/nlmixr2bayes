# Package index

## All functions

- [`adviControl()`](https://nlmixr2.github.io/nlmixr2bayes/reference/adviControl.md)
  : Control options for \`est="advi"\` (Stan's ADVI variational Bayes)
- [`getValidNlmixrCtl(`*`<stan>`*`)`](https://nlmixr2.github.io/nlmixr2bayes/reference/getValidNlmixrCtl.stan.md)
  : Validate a stanControl for the nlmixr2 dispatch (mandatory: the
  default errors before \`nlmixr2Est\` dispatch)
- [`nlmixr2bayesBuildInfo()`](https://nlmixr2.github.io/nlmixr2bayes/reference/nlmixr2bayesBuildInfo.md)
  : Probe the Stan build environment
- [`nutsControl()`](https://nlmixr2.github.io/nlmixr2bayes/reference/nutsControl.md)
  : Control options for \`est="nuts"\` (Stan's NUTS sampler)
- [`pathfinderControl()`](https://nlmixr2.github.io/nlmixr2bayes/reference/pathfinderControl.md)
  : Control options for \`est="pathfinder"\` (Stan's multi-path
  Pathfinder)
- [`rxsCheckGradient()`](https://nlmixr2.github.io/nlmixr2bayes/reference/rxsCheckGradient.md)
  : Compare a Stan model's gradient against finite differences of its
  log density
- [`rxsFastAvailable()`](https://nlmixr2.github.io/nlmixr2bayes/reference/rxsFastAvailable.md)
  : Is Path A currently armed for this handle?
- [`rxsHandles()`](https://nlmixr2.github.io/nlmixr2bayes/reference/rxsHandles.md)
  : List registered handles
- [`rxsIncludePath()`](https://nlmixr2.github.io/nlmixr2bayes/reference/rxsIncludePath.md)
  : Path to the rxstan C++ headers
- [`rxsInit()`](https://nlmixr2.github.io/nlmixr2bayes/reference/rxsInit.md)
  : Build a Stan \`init\` function from a model's \`ini()\` values
- [`rxsInvalidateFast()`](https://nlmixr2.github.io/nlmixr2bayes/reference/rxsInvalidateFast.md)
  : Force the next solve back through the slow path
- [`rxsProbeRxode2()`](https://nlmixr2.github.io/nlmixr2bayes/reference/rxsProbeRxode2.md)
  : Report whether rxode2's C function-pointer table was installed
- [`rxsRegister()`](https://nlmixr2.github.io/nlmixr2bayes/reference/rxsRegister.md)
  : Register an rxode2 model as a Stan-callable solver
- [`rxsRelease()`](https://nlmixr2.github.io/nlmixr2bayes/reference/rxsRelease.md)
  : Release a registered handle
- [`rxsReleaseAll()`](https://nlmixr2.github.io/nlmixr2bayes/reference/rxsReleaseAll.md)
  : Release every registered handle
- [`rxsResetStats()`](https://nlmixr2.github.io/nlmixr2bayes/reference/rxsResetStats.md)
  : Reset a handle's solve counters
- [`rxsSolve()`](https://nlmixr2.github.io/nlmixr2bayes/reference/rxsSolve.md)
  : Solve a registered handle from R
- [`rxsSolveStats()`](https://nlmixr2.github.io/nlmixr2bayes/reference/rxsSolveStats.md)
  : How often a handle has solved, and how often that failed
- [`rxsStanFromUi()`](https://nlmixr2.github.io/nlmixr2bayes/reference/rxsStanFromUi.md)
  : Generate a Stan program and a bridge handle from an nlmixr2 model
- [`rxsStanIncludes()`](https://nlmixr2.github.io/nlmixr2bayes/reference/rxsStanIncludes.md)
  : Build the \`includes=\` string for \[rstan::stan_model()\]
- [`rxsStanModel()`](https://nlmixr2.github.io/nlmixr2bayes/reference/rxsStanModel.md)
  : Compile a Stan program against the rxstan bridge
- [`stanCompile()`](https://nlmixr2.github.io/nlmixr2bayes/reference/stanCompile.md)
  : Compile a Stan program against the linked likelihood
- [`stanControl()`](https://nlmixr2.github.io/nlmixr2bayes/reference/stanControl.md)
  : Control options for the \`est="stan"\` estimation method
- [`stanLinkFree()`](https://nlmixr2.github.io/nlmixr2bayes/reference/stanLinkFree.md)
  : Free the linked likelihood
- [`stanLinkSetup()`](https://nlmixr2.github.io/nlmixr2bayes/reference/stanLinkSetup.md)
  : Load the linked likelihood for a model/data pair
- [`stanPopLinkSetup()`](https://nlmixr2.github.io/nlmixr2bayes/reference/stanPopLinkSetup.md)
  : Load the tier-0 (population-only) likelihood for a model/data pair
- [`stanPriors()`](https://nlmixr2.github.io/nlmixr2bayes/reference/stanPriors.md)
  : Translate a model's \`ini()\` priors into Stan
- [`stanSetCores()`](https://nlmixr2.github.io/nlmixr2bayes/reference/stanSetCores.md)
  : Set the subject-parallel thread count for the linked likelihood
