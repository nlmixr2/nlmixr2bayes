// Reference implementation of the same PK model using Stan's own ODE solver
// and its autodiff sensitivities, for benchmarking against the rxstan bridge.
//
// Matched to inst/stan/pk_1cmt_oral.stan: same parameters, same likelihood,
// and the single bolus expressed as an initial condition, which is the only
// way Stan can express a dose.
functions {
  vector pk_rhs(real t, vector y, real ka, real ke) {
    vector[2] dydt;
    dydt[1] = -ka * y[1];
    dydt[2] = ka * y[1] - ke * y[2];
    return dydt;
  }
}
data {
  int<lower=1> nObs;
  array[nObs] real ts;
  vector<lower=0>[nObs] cpObs;
  real<lower=0> amt;
  real<lower=0> relTol;
  real<lower=0> absTol;
  int<lower=1> maxSteps;
  int<lower=0, upper=1> stiff;
}
parameters {
  real lka;
  real lcl;
  real lv;
  real<lower=0> sigma;
}
model {
  vector[2] y0 = [amt, 0.0]';
  real ka = exp(lka);
  real ke = exp(lcl - lv);

  array[nObs] vector[2] ys;
  if (stiff == 1) {
    ys = ode_bdf_tol(pk_rhs, y0, 0.0, ts, relTol, absTol, maxSteps, ka, ke);
  } else {
    ys = ode_rk45_tol(pk_rhs, y0, 0.0, ts, relTol, absTol, maxSteps, ka, ke);
  }

  vector[nObs] cp;
  for (i in 1 : nObs) {
    cp[i] = ys[i][2] / exp(lv);
  }

  lka ~ normal(0, 1);
  lcl ~ normal(1.4, 1);
  lv ~ normal(3.4, 1);
  sigma ~ normal(0, 1);

  log(cpObs) ~ normal(log(cp), sigma);
  target += -sum(log(cpObs));  // Jacobian of the log transform
}
