#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

SEXP iniNlmixr2estFocei(SEXP p);
SEXP iniNlmixr2estNlm(SEXP p);
SEXP _nlmixr2bayes_setCores(SEXP n);
SEXP _nlmixr2bayes_apiVersion(void);
SEXP _nlmixr2bayes_nMix(void);
SEXP _nlmixr2bayes_nlmApiVersion(void);
SEXP _nlmixr2bayes_nlmDims(void);
SEXP _nlmixr2bayes_popEval(SEXP thetaS);
SEXP _nlmixr2bayes_dims(void);
SEXP _nlmixr2bayes_setTheta(SEXP thetaS);
SEXP _nlmixr2bayes_setOmegaInv(SEXP m);
SEXP _nlmixr2bayes_condBatch(SEXP etaS);
SEXP _nlmixr2bayes_setThetaBase(SEXP thetaS);
SEXP _nlmixr2bayes_setMuRef(SEXP idxS);
SEXP _nlmixr2bayes_setMuRefCov(SEXP thetaIdxS, SEXP etaIdxS, SEXP covValS);
SEXP _nlmixr2bayes_clearThetaBase(void);
SEXP _nlmixr2bayes_condBatchTheta(SEXP thetaS, SEXP etaS);
int nlmixr2bayes_cond_batch(const double *eta, int nid, int neta,
                           double *value, double *grad);
int nlmixr2bayes_pop_eval(const double *theta, int ntheta,
                         double *value, double *dTheta);
int nlmixr2bayes_cond_batch_theta(const double *theta, int ntheta,
                                 const double *eta, int nid, int neta,
                                 double *value, double *gradEta,
                                 double *gradTheta);
int nlmixr2bayes_iter_tick(const double *par, int n, double objf);

static const R_CallMethodDef CallEntries[] = {
  {"_nlmixr2bayes_iniFoceiPtrs", (DL_FUNC) &iniNlmixr2estFocei, 1},
  {"_nlmixr2bayes_iniNlmPtrs", (DL_FUNC) &iniNlmixr2estNlm, 1},
  {"_nlmixr2bayes_nlmApiVersion", (DL_FUNC) &_nlmixr2bayes_nlmApiVersion, 0},
  {"_nlmixr2bayes_nlmDims", (DL_FUNC) &_nlmixr2bayes_nlmDims, 0},
  {"_nlmixr2bayes_popEval", (DL_FUNC) &_nlmixr2bayes_popEval, 1},
  {"_nlmixr2bayes_setCores", (DL_FUNC) &_nlmixr2bayes_setCores, 1},
  {"_nlmixr2bayes_apiVersion", (DL_FUNC) &_nlmixr2bayes_apiVersion, 0},
  {"_nlmixr2bayes_nMix", (DL_FUNC) &_nlmixr2bayes_nMix, 0},
  {"_nlmixr2bayes_dims", (DL_FUNC) &_nlmixr2bayes_dims, 0},
  {"_nlmixr2bayes_setTheta", (DL_FUNC) &_nlmixr2bayes_setTheta, 1},
  {"_nlmixr2bayes_setOmegaInv", (DL_FUNC) &_nlmixr2bayes_setOmegaInv, 1},
  {"_nlmixr2bayes_condBatch", (DL_FUNC) &_nlmixr2bayes_condBatch, 1},
  {"_nlmixr2bayes_setThetaBase", (DL_FUNC) &_nlmixr2bayes_setThetaBase, 1},
  {"_nlmixr2bayes_setMuRef", (DL_FUNC) &_nlmixr2bayes_setMuRef, 1},
  {"_nlmixr2bayes_setMuRefCov", (DL_FUNC) &_nlmixr2bayes_setMuRefCov, 3},
  {"_nlmixr2bayes_clearThetaBase", (DL_FUNC) &_nlmixr2bayes_clearThetaBase, 0},
  {"_nlmixr2bayes_condBatchTheta", (DL_FUNC) &_nlmixr2bayes_condBatchTheta, 2},
  {NULL, NULL, 0}
};

void R_init_nlmixr2bayes(DllInfo *dll) {
  /* the symbols the runtime-compiled Stan model resolves */
  R_RegisterCCallable("nlmixr2bayes", "nlmixr2bayes_cond_batch",
                      (DL_FUNC) &nlmixr2bayes_cond_batch);
  R_RegisterCCallable("nlmixr2bayes", "nlmixr2bayes_cond_batch_theta",
                      (DL_FUNC) &nlmixr2bayes_cond_batch_theta);
  R_RegisterCCallable("nlmixr2bayes", "nlmixr2bayes_pop_eval",
                      (DL_FUNC) &nlmixr2bayes_pop_eval);
  R_RegisterCCallable("nlmixr2bayes", "nlmixr2bayes_iter_tick",
                      (DL_FUNC) &nlmixr2bayes_iter_tick);
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
