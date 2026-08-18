// Population one-compartment oral PK with between-subject random effects.
//
// This is the shape nlmixr2 models actually take, and the reason the bridge
// hands back a block-diagonal Jacobian: every subject gets its own parameter
// block, and rxode2 solves all of them in a single call.
//
// As in the single-subject program, the concentration is formed HERE rather
// than in rxode2, because rxode2 emits sensitivities for states only.
functions {
  vector rx_solve(int handle, vector p);
}
data {
  int<lower=1> nSub;
  int<lower=1> nObs;
  int<lower=1> handle;
  vector<lower=0>[nObs] cpObs;
  array[nObs] int<lower=1, upper=nSub> subj;
}
parameters {
  vector[3] theta;                 // population log ka, log cl, log v
  vector<lower=0>[3] omega;        // between-subject SD
  matrix[3, nSub] z;               // non-centered random effects
  real<lower=0> sigma;
}
transformed parameters {
  matrix[3, nSub] phi;
  for (s in 1 : nSub) {
    phi[ : , s] = theta + omega .* z[ : , s];
  }
}
model {
  vector[3 * nSub] p;
  for (s in 1 : nSub) {
    p[((s - 1) * 3 + 1) : (s * 3)] = phi[ : , s];
  }

  vector[nObs] center = rx_solve(handle, p);
  vector[nObs] cp;
  for (i in 1 : nObs) {
    cp[i] = center[i] / exp(phi[3, subj[i]]);
  }

  theta[1] ~ normal(0, 1);
  theta[2] ~ normal(1.4, 1);
  theta[3] ~ normal(3.4, 1);
  omega ~ normal(0, 0.5);
  to_vector(z) ~ std_normal();
  sigma ~ normal(0, 1);

  log(cpObs) ~ normal(log(cp), sigma);
  target += -sum(log(cpObs));  // Jacobian of the log transform
}
