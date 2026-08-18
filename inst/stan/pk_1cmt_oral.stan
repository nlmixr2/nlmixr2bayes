// One-compartment oral PK, population-free, fitted through rxstan.
//
// rx_solve() is declared but not defined here: --allow-undefined leaves it to
// the C++ bridge, which returns the `center` state at every observation time
// together with rxode2's analytic d center / d (lka, lcl, lv).
//
// Note that the concentration is formed HERE, not in rxode2.  rxode2 emits
// sensitivities for states only, so dividing by v is left to Stan, which
// differentiates that step itself.
functions {
  vector rx_solve(int handle, vector p);
}
data {
  int<lower=1> nObs;
  int<lower=1> handle;
  vector<lower=0>[nObs] cpObs;
}
parameters {
  real lka;
  real lcl;
  real lv;
  real<lower=0> sigma;
}
model {
  vector[3] p = [lka, lcl, lv]';
  vector[nObs] center = rx_solve(handle, p);
  vector[nObs] cp = center / exp(lv);

  lka ~ normal(0, 1);
  lcl ~ normal(1.4, 1);
  lv ~ normal(3.4, 1);
  sigma ~ normal(0, 1);

  log(cpObs) ~ normal(log(cp), sigma);
  target += -sum(log(cpObs));  // Jacobian of the log transform
}
