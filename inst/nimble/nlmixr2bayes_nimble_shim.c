/* Standalone shim consumed by NIMBLE's nimbleExternalCall(): no
 * rxode2/nlmixr2est headers and no dependency on any nlmixr2bayes build
 * artifact -- just R.h/Rdynload.h, resolving the registered C-callable at
 * first call, exactly like inst/include/nlmixr2bayes_lp.hpp does for Stan.
 *
 * Compiled on demand by .nimbleShimCompile() (R/nimbleShim.R) into a
 * throwaway .o that nimbleExternalCall(oFile=) links against.  Linking
 * directly against a package build's own .o/.so instead would be WRONG: it
 * creates a second, disconnected copy of nlmixr2bayesPtr.c's translation
 * unit (with its own static FOCEi pointer-table globals) that the real
 * R-side link setup (stanLinkSetup() + the theta-base/mu-ref .Call()s) can
 * never reach, since that setup only ever touches the actual loaded
 * nlmixr2bayes DLL.  Going through R_GetCCallable is the only way to reach
 * that real, singleton instance. */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

typedef int (*nlmixr2bayes_cond_batch_theta_t)(
    const double *theta, int ntheta, const double *eta, int nid, int neta,
    double *value, double *gradEta, double *gradTheta);

#ifdef __cplusplus
extern "C" {
#endif

int nlmixr2bayes_nimble_cond_batch_theta(const double *theta, int ntheta,
                                          const double *eta, int nid,
                                          int neta, double *value,
                                          double *gradEta,
                                          double *gradTheta) {
  /* Cached for this DSO's lifetime, same tradeoff inst/include/
   * nlmixr2bayes_lp.hpp already accepts for Stan: unloading and reloading
   * nlmixr2bayes (e.g. repeated pkgload::load_all() in one R session) while
   * a compiled model built against the OLD load is still alive would leave
   * this pointing at unmapped memory. Not addressed here, matching the
   * existing accepted precedent -- a fix would need to cover both entry
   * points together. */
  static nlmixr2bayes_cond_batch_theta_t fn = NULL;
  if (fn == NULL) {
    fn = (nlmixr2bayes_cond_batch_theta_t)
      R_GetCCallable("nlmixr2bayes", "nlmixr2bayes_cond_batch_theta");
  }
  if (fn == NULL) return -200; /* nlmixr2bayes not loaded / symbol missing */
  return fn(theta, ntheta, eta, nid, neta, value, gradEta, gradTheta);
}

#ifdef __cplusplus
}
#endif
