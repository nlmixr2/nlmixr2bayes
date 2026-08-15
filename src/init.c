#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

SEXP iniNlmixr2estFocei(SEXP p);
SEXP _nlmixr2stan_setCores(SEXP n);
SEXP _nlmixr2stan_apiVersion(void);
SEXP _nlmixr2stan_dims(void);
SEXP _nlmixr2stan_setTheta(SEXP thetaS);
SEXP _nlmixr2stan_setOmegaInv(SEXP m);
SEXP _nlmixr2stan_condBatch(SEXP etaS);
SEXP _nlmixr2stan_setThetaBase(SEXP thetaS);
SEXP _nlmixr2stan_setMuRef(SEXP idxS);
SEXP _nlmixr2stan_setMuRefCov(SEXP thetaIdxS, SEXP etaIdxS, SEXP covValS);
SEXP _nlmixr2stan_clearThetaBase(void);
SEXP _nlmixr2stan_condBatchTheta(SEXP thetaS, SEXP etaS);
int nlmixr2stan_cond_batch(const double *eta, int nid, int neta,
                           double *value, double *grad);
int nlmixr2stan_cond_batch_theta(const double *theta, int ntheta,
                                 const double *eta, int nid, int neta,
                                 double *value, double *gradEta,
                                 double *gradTheta);

static const R_CallMethodDef CallEntries[] = {
  {"_nlmixr2stan_iniFoceiPtrs", (DL_FUNC) &iniNlmixr2estFocei, 1},
  {"_nlmixr2stan_setCores", (DL_FUNC) &_nlmixr2stan_setCores, 1},
  {"_nlmixr2stan_apiVersion", (DL_FUNC) &_nlmixr2stan_apiVersion, 0},
  {"_nlmixr2stan_dims", (DL_FUNC) &_nlmixr2stan_dims, 0},
  {"_nlmixr2stan_setTheta", (DL_FUNC) &_nlmixr2stan_setTheta, 1},
  {"_nlmixr2stan_setOmegaInv", (DL_FUNC) &_nlmixr2stan_setOmegaInv, 1},
  {"_nlmixr2stan_condBatch", (DL_FUNC) &_nlmixr2stan_condBatch, 1},
  {"_nlmixr2stan_setThetaBase", (DL_FUNC) &_nlmixr2stan_setThetaBase, 1},
  {"_nlmixr2stan_setMuRef", (DL_FUNC) &_nlmixr2stan_setMuRef, 1},
  {"_nlmixr2stan_setMuRefCov", (DL_FUNC) &_nlmixr2stan_setMuRefCov, 3},
  {"_nlmixr2stan_clearThetaBase", (DL_FUNC) &_nlmixr2stan_clearThetaBase, 0},
  {"_nlmixr2stan_condBatchTheta", (DL_FUNC) &_nlmixr2stan_condBatchTheta, 2},
  {NULL, NULL, 0}
};

void R_init_nlmixr2stan(DllInfo *dll) {
  /* the symbols the runtime-compiled Stan model resolves */
  R_RegisterCCallable("nlmixr2stan", "nlmixr2stan_cond_batch",
                      (DL_FUNC) &nlmixr2stan_cond_batch);
  R_RegisterCCallable("nlmixr2stan", "nlmixr2stan_cond_batch_theta",
                      (DL_FUNC) &nlmixr2stan_cond_batch_theta);
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
