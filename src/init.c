#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

SEXP iniNlmixr2estFocei(SEXP p);
SEXP iniNlmixr2estNlm(SEXP p);
SEXP iniRxodePtrs(SEXP p);

/* rxstan bridge entry points (bridge.cpp, fast.cpp) */
SEXP rxstanProbeRxode2(void);
SEXP rxstanSetDims(SEXP handleSXP, SEXP nySXP, SEXP nBlockSXP, SEXP nBlocksSXP,
                   SEXP outBlockSXP);
SEXP rxstanClearDims(SEXP handleSXP);
SEXP rxstanSolve(SEXP handleSXP, SEXP pSXP);
SEXP rxstanSolveSlow(SEXP handleSXP, SEXP pSXP);
SEXP rxstanStats(SEXP handleSXP);
SEXP rxstanResetStats(SEXP handleSXP);
SEXP rxstanSetSilent(SEXP silentSXP);
SEXP rxstanFastSetup(SEXP handleSXP, SEXP sensIdxSXP, SEXP outIdxSXP,
                     SEXP sensStateSXP, SEXP nobsSXP, SEXP blockOfSXP,
                     SEXP quietSXP);
SEXP rxstanFastInvalidate(void);
SEXP rxstanFastAvailable(SEXP handleSXP);
SEXP rxstanProbeSolve(void);
SEXP rxstanProbeParams(void);
SEXP rxstanProbeStates(void);

/* C ABI symbols the Stan translation unit resolves via R_GetCCallable */
int rxs_dims(int handle, int *ny, int *np);
int rxs_layout(int handle, int *ny, int *np, int *nBlock, int *nBlocks);
int rxs_out_block(int handle, int *out, int n);
int rxs_stats(int handle, double *nSolve, double *nFail);
int rxs_solve_sens(int handle, const double *p, int np, double *y, double *dydp,
                   int ny, char *errbuf, int errlen);
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
SEXP _nlmixr2bayes_resetEvalCount(void);
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
  {"_nlmixr2bayes_iniRxodePtrs", (DL_FUNC) &iniRxodePtrs, 1},
  /* rxstan bridge */
  {"iniRxodePtrs",              (DL_FUNC) &iniRxodePtrs,           1},
  {"C_rxstanProbeRxode2",       (DL_FUNC) &rxstanProbeRxode2,      0},
  {"C_rxstanSetDims",           (DL_FUNC) &rxstanSetDims,          5},
  {"C_rxstanClearDims",         (DL_FUNC) &rxstanClearDims,        1},
  {"C_rxstanSolve",             (DL_FUNC) &rxstanSolve,            2},
  {"C_rxstanSolveSlow",         (DL_FUNC) &rxstanSolveSlow,        2},
  {"C_rxstanStats",             (DL_FUNC) &rxstanStats,            1},
  {"C_rxstanResetStats",        (DL_FUNC) &rxstanResetStats,       1},
  {"C_rxstanSetSilent",         (DL_FUNC) &rxstanSetSilent,        1},
  {"C_rxstanFastSetup",         (DL_FUNC) &rxstanFastSetup,        7},
  {"C_rxstanFastInvalidate",    (DL_FUNC) &rxstanFastInvalidate,   0},
  {"C_rxstanFastAvailable",     (DL_FUNC) &rxstanFastAvailable,    1},
  {"C_rxstanProbeSolve",        (DL_FUNC) &rxstanProbeSolve,       0},
  {"C_rxstanProbeParams",       (DL_FUNC) &rxstanProbeParams,      0},
  {"C_rxstanProbeStates",       (DL_FUNC) &rxstanProbeStates,      0},
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
  {"_nlmixr2bayes_resetEvalCount", (DL_FUNC) &_nlmixr2bayes_resetEvalCount, 0},
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
  /* rxstan bridge C ABI: the Stan translation unit resolves these */
  R_RegisterCCallable("nlmixr2bayes", "rxs_dims",       (DL_FUNC) &rxs_dims);
  R_RegisterCCallable("nlmixr2bayes", "rxs_layout",     (DL_FUNC) &rxs_layout);
  R_RegisterCCallable("nlmixr2bayes", "rxs_out_block",  (DL_FUNC) &rxs_out_block);
  R_RegisterCCallable("nlmixr2bayes", "rxs_stats",      (DL_FUNC) &rxs_stats);
  R_RegisterCCallable("nlmixr2bayes", "rxs_solve_sens", (DL_FUNC) &rxs_solve_sens);
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
