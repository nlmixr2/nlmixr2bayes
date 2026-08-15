#ifndef NLMIXR2STAN_LP_HPP
#define NLMIXR2STAN_LP_HPP
// Injected by rstan::stan_model(includes=) IMMEDIATELY BEFORE `class model...`,
// i.e. INSIDE namespace model<hash>_namespace.  Consequences, all load-bearing:
//   * NO #include of any kind here: an #include inside a namespace is a
//     disaster.  Everything needed is already at global scope -- rstan
//     prepends <Rcpp.h> (which pulls <R_ext/Rdynload.h>, declaring
//     R_GetCCallable) and <stan/model/model_header.hpp> (Eigen,
//     stan::math::var, precomputed_gradients, throw_domain_error).
//   * everything here must be legal inside a namespace.
//   * each compiled model gets its own copy in its own namespace -- no ODR
//     question, and no Stan symbol is ever exported from nlmixr2stan's own
//     shared object (which contains no Stan code at all).
//
// The C ABI it resolves lives in nlmixr2stan's DLL
// (src/nlmixr2stanPtr.c::nlmixr2stan_cond_batch), which forwards to
// nlmixr2est's FOCEi conditional-likelihood entry-point table: value is
// log p(y_i | eta_i) (identical to foceiLikRun(type="cond")) and grad is
// d/d(eta_i) of that value, assembled inside nlmixr2est so the sign
// convention is never reconstructed here.  Row-major nid x neta.

typedef int (*nlmixr2stan_cond_batch_t)(const double *eta, int nid, int neta,
                                        double *value, double *grad);

inline nlmixr2stan_cond_batch_t nlmixr2stan__fn() {
  static nlmixr2stan_cond_batch_t fn = 0;
  if (fn == 0) {
    // R_GetCCallable Rf_error()s (longjmp) if nlmixr2stan is not loaded.
    // That can only happen on the FIRST log_prob call, before any sampling
    // work; the R-side wrapper forces one evaluation up front so the failure
    // surfaces there rather than mid-chain.
    fn = reinterpret_cast<nlmixr2stan_cond_batch_t>(
      R_GetCCallable("nlmixr2stan", "nlmixr2stan_cond_batch"));
  }
  return fn;
}

// Shared driver: copy the (column-limited) Eigen matrix out row-major, call
// the C ABI, surface hard failures as a Stan domain error (which Stan treats
// as a rejection, not a crash).  rc > 0 means some subjects' values are -Inf
// with zeroed gradient rows -- also a rejection, handled by Stan when the
// -Inf reaches the target.
inline int nlmixr2stan__batch(const double *eta, int nid, int neta,
                              double *value, double *grad,
                              std::ostream *pstream__) {
  (void)pstream__;
  const int rc = nlmixr2stan__fn()(eta, nid, neta, value, grad);
  if (rc < 0) {
    stan::math::throw_domain_error("nlmixr2_cond_all", "nlmixr2est status",
                                   rc, "was ", " (<0 is a hard failure)");
  }
  return rc;
}

// --- double overload: log_prob<false,false,double>, write_array, and
// rstan::log_prob ------------------------------------------------------------
inline Eigen::Matrix<double, Eigen::Dynamic, 1>
nlmixr2_cond_all(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& eta,
                 std::ostream* pstream__) {
  const int nid = static_cast<int>(eta.rows());
  const int neta = static_cast<int>(eta.cols());
  std::vector<double> e(static_cast<size_t>(nid) * neta);
  std::vector<double> v(static_cast<size_t>(nid));
  std::vector<double> g(static_cast<size_t>(nid) * neta);
  for (int i = 0; i < nid; ++i) {
    for (int j = 0; j < neta; ++j) {
      e[static_cast<size_t>(i) * neta + j] = eta(i, j);
    }
  }
  nlmixr2stan__batch(e.data(), nid, neta, v.data(), g.data(), pstream__);
  Eigen::Matrix<double, Eigen::Dynamic, 1> out(nid);
  for (int i = 0; i < nid; ++i) out(i) = v[i];
  return out;
}

// --- var overload: the NUTS path --------------------------------------------
// One vari per subject via precomputed_gradients, carrying that subject's
// neta partials.  Nothing in nlmixr2est is templated on var; this is the
// entire autodiff surface, and Stan's own autodiff composes it onward (e.g.
// through eta = L * z in the transformed parameters block).
inline Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>
nlmixr2_cond_all(const Eigen::Matrix<stan::math::var, Eigen::Dynamic,
                                     Eigen::Dynamic>& eta,
                 std::ostream* pstream__) {
  using stan::math::var;
  const int nid = static_cast<int>(eta.rows());
  const int neta = static_cast<int>(eta.cols());
  std::vector<double> e(static_cast<size_t>(nid) * neta);
  std::vector<double> v(static_cast<size_t>(nid));
  std::vector<double> g(static_cast<size_t>(nid) * neta);
  for (int i = 0; i < nid; ++i) {
    for (int j = 0; j < neta; ++j) {
      e[static_cast<size_t>(i) * neta + j] = eta(i, j).val();
    }
  }
  nlmixr2stan__batch(e.data(), nid, neta, v.data(), g.data(), pstream__);
  Eigen::Matrix<var, Eigen::Dynamic, 1> out(nid);
  std::vector<var> ops(static_cast<size_t>(neta));
  std::vector<double> gi(static_cast<size_t>(neta));
  for (int i = 0; i < nid; ++i) {
    for (int j = 0; j < neta; ++j) {
      ops[j] = eta(i, j);
      gi[j] = g[static_cast<size_t>(i) * neta + j];
    }
    out(i) = stan::math::precomputed_gradients(v[i], ops, gi);
  }
  return out;
}

// ============================================================================
// Tier 2: theta sampled too.  nlmixr2stan_cond_batch_theta writes the
// natural-scale theta into the loaded problem (the omega tail of the
// parameter vector stays at its link-time values; the conditional does not
// depend on it) and returns value + d/d(eta) + d/d(theta), the theta columns
// combining nlmixr2est's forward sensitivities (non-mu structural and
// residual thetas) with the mu-reference identity d/dtheta_p = d/deta_k.

typedef int (*nlmixr2stan_cond_batch_theta_t)(const double *theta, int ntheta,
                                              const double *eta, int nid,
                                              int neta, double *value,
                                              double *gradEta,
                                              double *gradTheta);

inline nlmixr2stan_cond_batch_theta_t nlmixr2stan__fnTheta() {
  static nlmixr2stan_cond_batch_theta_t fn = 0;
  if (fn == 0) {
    fn = reinterpret_cast<nlmixr2stan_cond_batch_theta_t>(
      R_GetCCallable("nlmixr2stan", "nlmixr2stan_cond_batch_theta"));
  }
  return fn;
}

inline int nlmixr2stan__batchTheta(const double *theta, int ntheta,
                                   const double *eta, int nid, int neta,
                                   double *value, double *gradEta,
                                   double *gradTheta,
                                   std::ostream *pstream__) {
  (void)pstream__;
  const int rc = nlmixr2stan__fnTheta()(theta, ntheta, eta, nid, neta,
                                        value, gradEta, gradTheta);
  if (rc < 0) {
    stan::math::throw_domain_error("nlmixr2_cond_all2", "nlmixr2est status",
                                   rc, "was ", " (<0 is a hard failure)");
  }
  return rc;
}

// --- double overload --------------------------------------------------------
inline Eigen::Matrix<double, Eigen::Dynamic, 1>
nlmixr2_cond_all2(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& eta,
                  const Eigen::Matrix<double, Eigen::Dynamic, 1>& theta,
                  std::ostream* pstream__) {
  const int nid = static_cast<int>(eta.rows());
  const int neta = static_cast<int>(eta.cols());
  const int nth = static_cast<int>(theta.size());
  std::vector<double> e(static_cast<size_t>(nid) * neta);
  std::vector<double> th(static_cast<size_t>(nth));
  std::vector<double> v(static_cast<size_t>(nid));
  std::vector<double> ge(static_cast<size_t>(nid) * neta);
  std::vector<double> gt(static_cast<size_t>(nid) * nth);
  for (int i = 0; i < nid; ++i) {
    for (int j = 0; j < neta; ++j) {
      e[static_cast<size_t>(i) * neta + j] = eta(i, j);
    }
  }
  for (int j = 0; j < nth; ++j) th[j] = theta(j);
  nlmixr2stan__batchTheta(th.data(), nth, e.data(), nid, neta,
                          v.data(), ge.data(), gt.data(), pstream__);
  Eigen::Matrix<double, Eigen::Dynamic, 1> out(nid);
  for (int i = 0; i < nid; ++i) out(i) = v[i];
  return out;
}

// --- var overload: one vari per subject carrying neta + ntheta partials -----
inline Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>
nlmixr2_cond_all2(const Eigen::Matrix<stan::math::var, Eigen::Dynamic,
                                      Eigen::Dynamic>& eta,
                  const Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>& theta,
                  std::ostream* pstream__) {
  using stan::math::var;
  const int nid = static_cast<int>(eta.rows());
  const int neta = static_cast<int>(eta.cols());
  const int nth = static_cast<int>(theta.size());
  std::vector<double> e(static_cast<size_t>(nid) * neta);
  std::vector<double> th(static_cast<size_t>(nth));
  std::vector<double> v(static_cast<size_t>(nid));
  std::vector<double> ge(static_cast<size_t>(nid) * neta);
  std::vector<double> gt(static_cast<size_t>(nid) * nth);
  for (int i = 0; i < nid; ++i) {
    for (int j = 0; j < neta; ++j) {
      e[static_cast<size_t>(i) * neta + j] = eta(i, j).val();
    }
  }
  for (int j = 0; j < nth; ++j) th[j] = theta(j).val();
  nlmixr2stan__batchTheta(th.data(), nth, e.data(), nid, neta,
                          v.data(), ge.data(), gt.data(), pstream__);
  Eigen::Matrix<var, Eigen::Dynamic, 1> out(nid);
  std::vector<var> ops(static_cast<size_t>(neta + nth));
  std::vector<double> gi(static_cast<size_t>(neta + nth));
  for (int j = 0; j < nth; ++j) ops[neta + j] = theta(j);
  for (int i = 0; i < nid; ++i) {
    for (int j = 0; j < neta; ++j) {
      ops[j] = eta(i, j);
      gi[j] = ge[static_cast<size_t>(i) * neta + j];
    }
    for (int j = 0; j < nth; ++j) {
      gi[neta + j] = gt[static_cast<size_t>(i) * nth + j];
    }
    out(i) = stan::math::precomputed_gradients(v[i], ops, gi);
  }
  return out;
}
#endif
