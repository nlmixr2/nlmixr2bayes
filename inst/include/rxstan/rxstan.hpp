#ifndef RXSTAN_RXSTAN_HPP
#define RXSTAN_RXSTAN_HPP

// The Stan-facing half of the bridge.  This header is spliced into the
// translation unit rstan generates for a model, and supplies the body of the
// `rx_solve` function that the .stan program declares but does not define.
//
// Stan never differentiates through the solver.  rxode2 returns the states and
// the analytic dY/dtheta it already computes, and those partials are pushed
// straight onto the reverse tape with precomputed_gradients().

// rxode2's solve structure is a PROCESS SINGLETON: getRxSolve_(), the dydt /
// Jacobian function pointers and the eventSens shape are all file-scope
// globals in rxode2, and the Layout cache below is a plain static.  Two
// threads inside the bridge at once corrupt all of them, and the result is
// silently wrong gradients that HMC will happily sample from.
//
// This CANNOT be checked with `#ifdef STAN_THREADS`: rstan compiles every
// model with -DSTAN_THREADS regardless of rstan_options("threads_per_chain"),
// so the macro says only that the AD tape is thread_local, not that anything
// runs concurrently.  The hazard is concurrent *execution*, so detect that
// directly -- see solve_guard below.

#include <stan/math/rev.hpp>

#include <R_ext/Rdynload.h>

#include <cstddef>
#include <atomic>
#include <map>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace rxstan {

// Refuses concurrent entry to the bridge.  rx_solve() is never recursive, so a
// non-zero depth on arrival means another thread is already inside -- i.e. the
// bridge was used under map_rect()/reduce_sum() or a threaded sampler.
// Throwing turns silent corruption into a rejected proposal with a message.
class solve_guard {
 public:
  solve_guard() {
    if (depth().fetch_add(1, std::memory_order_acq_rel) != 0) {
      depth().fetch_sub(1, std::memory_order_acq_rel);
      throw std::domain_error(
          "rxstan: rx_solve() was entered from two threads at once. rxode2's "
          "solve structure is a process singleton, so it cannot be called "
          "inside map_rect() or reduce_sum(), or from a threaded sampler. Run "
          "chains as separate processes and leave subject-level parallelism to "
          "rxode2.");
    }
  }
  ~solve_guard() { depth().fetch_sub(1, std::memory_order_acq_rel); }

  solve_guard(const solve_guard &) = delete;
  solve_guard &operator=(const solve_guard &) = delete;

 private:
  static std::atomic<int> &depth() {
    static std::atomic<int> d(0);
    return d;
  }
};

// --- C ABI exposed by the nlmixr2bayes package -----------------------------
//
// The bridge used to ship as a package called `rxstan`; it now lives inside
// nlmixr2bayes, and src/init.c registers these entry points under that name
// (R_RegisterCCallable keys on the string, so the two have to agree).

typedef int (*rxs_layout_t)(int, int *, int *, int *, int *);
typedef int (*rxs_out_block_t)(int, int *, int);
typedef int (*rxs_solve_sens_t)(int, const double *, int, double *, double *,
                                int, char *, int);

inline rxs_layout_t layout_fn() {
  static rxs_layout_t f =
      reinterpret_cast<rxs_layout_t>(R_GetCCallable("nlmixr2bayes", "rxs_layout"));
  return f;
}

inline rxs_out_block_t out_block_fn() {
  static rxs_out_block_t f = reinterpret_cast<rxs_out_block_t>(
      R_GetCCallable("nlmixr2bayes", "rxs_out_block"));
  return f;
}

inline rxs_solve_sens_t solve_fn() {
  static rxs_solve_sens_t f = reinterpret_cast<rxs_solve_sens_t>(
      R_GetCCallable("nlmixr2bayes", "rxs_solve_sens"));
  return f;
}

// A population model gives every subject its own block of parameters.  Since
// subjects are independent the Jacobian is block diagonal, so output i depends
// only on block outBlock[i] -- and each autodiff node needs just that block's
// operands rather than the whole parameter vector.
struct Layout {
  int ny = 0;
  int np = 0;
  int nBlock = 0;
  int nBlocks = 0;
  std::vector<int> outBlock;
};

inline const Layout &layout(int handle) {
  // Deliberately not thread_local: solve_guard makes concurrent entry an
  // error, and every caller reaches this through rx_solve().
  static std::map<int, Layout> cache;
  std::map<int, Layout>::const_iterator it = cache.find(handle);
  if (it != cache.end()) return it->second;

  Layout L;
  if (layout_fn()(handle, &L.ny, &L.np, &L.nBlock, &L.nBlocks) != 0) {
    std::ostringstream m;
    m << "rxstan: unknown handle " << handle;
    throw std::domain_error(m.str());
  }
  L.outBlock.resize(L.ny);
  if (out_block_fn()(handle, L.outBlock.data(), L.ny) != 0) {
    throw std::domain_error("rxstan: could not read the output block map");
  }
  return cache.insert(std::make_pair(handle, L)).first->second;
}

// Values plus analytic sensitivities, straight from rxode2.  A failed solve
// becomes a domain_error so HMC rejects the proposal instead of aborting.
inline void solve_analytic(int handle, const Eigen::VectorXd &p,
                           Eigen::VectorXd &y, Eigen::MatrixXd &dydp) {
  const Layout &L = layout(handle);
  if (p.size() != L.np) {
    std::ostringstream m;
    m << "rxstan: handle " << handle << " expects " << L.np
      << " parameters, got " << p.size();
    throw std::invalid_argument(m.str());
  }

  y.resize(L.ny);
  dydp.resize(L.ny, L.nBlock);

  char err[512];
  err[0] = '\0';
  const int rc = solve_fn()(handle, p.data(), L.np, y.data(), dydp.data(), L.ny,
                            err, static_cast<int>(sizeof(err)));
  if (rc != 0) {
    std::ostringstream m;
    m << "rxstan: solve failed (" << rc << "): " << err;
    throw std::domain_error(m.str());
  }
}

inline Eigen::VectorXd solve_value(int handle, const Eigen::VectorXd &p) {
  Eigen::VectorXd y;
  Eigen::MatrixXd dydp;
  solve_analytic(handle, p, y, dydp);
  return y;
}

// Central differences over the solver.  This is the fallback the primary
// solve_policy template uses, and doubles as the independent oracle the
// gradient tests check the analytic path against.  Only the block a given
// output belongs to is perturbed, matching the analytic layout.
inline void solve_finite_diff(int handle, const Eigen::VectorXd &p,
                              Eigen::VectorXd &y, Eigen::MatrixXd &dydp) {
  const Layout &L = layout(handle);

  y = solve_value(handle, p);
  dydp.resize(L.ny, L.nBlock);

  for (int j = 0; j < L.nBlock; ++j) {
    Eigen::VectorXd pp = p;
    Eigen::VectorXd pm = p;
    std::vector<double> h(L.nBlocks);
    for (int b = 0; b < L.nBlocks; ++b) {
      const int k = b * L.nBlock + j;
      h[b] = 1e-6 * std::max(1.0, std::abs(p(k)));
      pp(k) += h[b];
      pm(k) -= h[b];
    }
    const Eigen::VectorXd up = solve_value(handle, pp);
    const Eigen::VectorXd um = solve_value(handle, pm);
    for (int i = 0; i < L.ny; ++i) {
      dydp(i, j) = (up(i) - um(i)) / (2 * h[L.outBlock[i]]);
    }
  }
}

// --- The specialization hook -----------------------------------------------
//
// Primary template: the default, no-analytic-information-supplied path.  A
// generated model header overrides it for its own tag,
//
//   struct my_model_tag {};
//   namespace rxstan {
//   template <> struct solve_policy<my_model_tag> {
//     static void jacobian(int h, const Eigen::VectorXd& p,
//                          Eigen::VectorXd& y, Eigen::MatrixXd& J) {
//       solve_analytic(h, p, y, J);
//     }
//   };
//   }
//   #define RXSTAN_MODEL_TAG my_model_tag
//
// which is what makes the analytic-vs-numeric choice a compile-time property
// of the model rather than a runtime branch.
template <typename Tag>
struct solve_policy {
  static void jacobian(int handle, const Eigen::VectorXd &p,
                       Eigen::VectorXd &y, Eigen::MatrixXd &dydp) {
    solve_finite_diff(handle, p, y, dydp);
  }
};

struct analytic_tag {};
struct finite_diff_tag {};

template <>
struct solve_policy<analytic_tag> {
  static void jacobian(int handle, const Eigen::VectorXd &p,
                       Eigen::VectorXd &y, Eigen::MatrixXd &dydp) {
    solve_analytic(handle, p, y, dydp);
  }
};

#ifndef RXSTAN_MODEL_TAG
#define RXSTAN_MODEL_TAG ::rxstan::analytic_tag
#endif

// --- Stan entry points ------------------------------------------------------
//
// Satisfies `vector rx_solve(int handle, vector p);` declared in the .stan
// functions block under --allow-undefined.

// Data-only instantiation: no tape, no sensitivities needed.
template <typename T, stan::require_col_vector_t<T> * = nullptr,
          stan::require_st_arithmetic<T> * = nullptr>
inline Eigen::VectorXd rx_solve(const int handle, const T &p,
                                std::ostream *pstream__) {
  (void)pstream__;
  const solve_guard guard;
  const Eigen::VectorXd pd = p;
  return solve_value(handle, pd);
}

// Autodiff instantiation.  Each output only depends on its own parameter
// block, so its node carries just that block's operands -- for a population
// model that is nBlock partials per output rather than nBlock * nSubjects.
template <typename T, stan::require_col_vector_t<T> * = nullptr,
          stan::require_st_var<T> * = nullptr>
inline Eigen::Matrix<stan::math::var, -1, 1> rx_solve(const int handle,
                                                      const T &p,
                                                      std::ostream *pstream__) {
  (void)pstream__;
  using stan::math::var;

  const solve_guard guard;
  const Eigen::Matrix<var, -1, 1> p_ref = p;
  const Eigen::VectorXd pd = stan::math::value_of(p_ref);

  Eigen::VectorXd y;
  Eigen::MatrixXd dydp;
  solve_policy<RXSTAN_MODEL_TAG>::jacobian(handle, pd, y, dydp);

  const Layout &L = layout(handle);
  std::vector<var> operands(L.nBlock);
  std::vector<double> grad(L.nBlock);

  Eigen::Matrix<var, -1, 1> out(y.size());
  for (Eigen::Index i = 0; i < y.size(); ++i) {
    const int base = L.outBlock[i] * L.nBlock;
    for (int j = 0; j < L.nBlock; ++j) {
      operands[j] = p_ref(base + j);
      grad[j] = dydp(i, j);
    }
    out(i) = stan::math::precomputed_gradients(y(i), operands, grad);
  }
  return out;
}

}  // namespace rxstan

#endif
