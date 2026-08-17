/* Consume nlmixr2est's FOCEi conditional-likelihood C API through its
 * external-pointer table (inst/include/nlmixr2estFoceiPtr.h, installed at
 * .onLoad from nlmixr2est::.nlmixr2estFoceiPtrs()), and re-export ONE plain-C
 * symbol via R_RegisterCCallable for the injected Stan header
 * (inst/include/nlmixr2bayes_lp.hpp) to find with R_GetCCallable.
 *
 * Deliberately NO Rcpp, NO Armadillo, NO Eigen, NO StanHeaders here: the
 * external function is compiled into rstan's own translation unit, and
 * everything crossing a shared-object boundary is a double pointer or an int,
 * so no C++ type ever appears in more than one DSO (see the design notes in
 * nlmixr2/nlmixr2bayes#1). */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <time.h>
#include <rxode2ptr.h>
iniRxode2ptr
#include <nlmixr2estFoceiPtr.h>
#include <nlmixr2estNlmPtr.h>

/* instantiates the 7 function-pointer globals + SEXP iniNlmixr2estFocei(SEXP) */
iniNlmixr2estFoceiGlobals

/* instantiates the 3 nlm (tier-0) globals + SEXP iniNlmixr2estNlm(SEXP) */
iniNlmixr2estNlmGlobals

/* Subject-parallel thread count for the batch entries, set from R before
 * sampling (rxode2/nlmixr2est parallelize over subjects inside each call;
 * Stan itself stays single-threaded). */
static int nlmixr2bayesCores = 1;

SEXP _nlmixr2bayes_setCores(SEXP n) {
  int c = Rf_asInteger(n);
  if (c == NA_INTEGER || c < 1) c = 1;
  nlmixr2bayesCores = c;
  return Rf_ScalarInteger(nlmixr2bayesCores);
}

/* Iteration tick from inside the compiled Stan model's model block: the
 * natural-scale display vector (thetas + actual omega entries) and the
 * current -2*sum(conditional log-lik).  Forwards to nlmixr2est's scale.h
 * iteration print (entry 8 of the FOCEi table); no-op (-1) when that
 * entry is absent (older nlmixr2est) or printing is not armed. */
int nlmixr2bayes_iter_tick(const double *par, int n, double objf) {
  if (nlmixr2FoceiIterPrintRowP == NULL) return -1;
  return nlmixr2FoceiIterPrintRowP(par, n, objf);
}

SEXP _nlmixr2bayes_apiVersion(void) {
  if (nlmixr2FoceiApiVersionP == NULL) return Rf_ScalarInteger(-1);
  return Rf_ScalarInteger(nlmixr2FoceiApiVersionP());
}

/* mixture component count of the loaded problem (#955): 1 non-mixture, K
 * mixture, -1 not loaded, -2 when the loaded nlmixr2est predates the entry
 * (treat as refuse-mixtures) */
SEXP _nlmixr2bayes_nMix(void) {
  if (nlmixr2FoceiNMixP == NULL) return Rf_ScalarInteger(-2);
  return Rf_ScalarInteger(nlmixr2FoceiNMixP());
}

/* ---- tier 0: population-only models via the nlm C API (#953) -------------
 * The nlm entry returns MINUS log-likelihood + its analytic theta gradient
 * (natural scale); the injected header negates.  -100 = table not
 * installed (older nlmixr2est without the nlm API). */
int nlmixr2bayes_pop_eval(const double *theta, int ntheta,
                         double *value, double *dTheta) {
  if (nlmixr2NlmEvalP == NULL) return -100;
  return nlmixr2NlmEvalP(theta, ntheta, value, dTheta);
}

SEXP _nlmixr2bayes_nlmApiVersion(void) {
  if (nlmixr2NlmApiVersionP == NULL) return Rf_ScalarInteger(-1);
  return Rf_ScalarInteger(nlmixr2NlmApiVersionP());
}

SEXP _nlmixr2bayes_nlmDims(void) {
  SEXP ret, nm;
  int ntheta = 0, nobs = 0, flags = 0, rc;
  if (nlmixr2NlmDimsP == NULL) Rf_error("nlmixr2est nlm pointer table not installed");
  rc = nlmixr2NlmDimsP(&ntheta, &nobs, &flags);
  ret = PROTECT(Rf_allocVector(INTSXP, 4));
  INTEGER(ret)[0] = rc;
  INTEGER(ret)[1] = ntheta;
  INTEGER(ret)[2] = nobs;
  INTEGER(ret)[3] = flags;
  nm = PROTECT(Rf_allocVector(STRSXP, 4));
  SET_STRING_ELT(nm, 0, Rf_mkChar("status"));
  SET_STRING_ELT(nm, 1, Rf_mkChar("ntheta"));
  SET_STRING_ELT(nm, 2, Rf_mkChar("nobs"));
  SET_STRING_ELT(nm, 3, Rf_mkChar("flags"));
  Rf_setAttrib(ret, R_NamesSymbol, nm);
  UNPROTECT(2);
  return ret;
}

/* R-callable mirror over the exact C path the compiled Stan model uses;
 * returns the RAW nlm convention (value = -logLik) -- the header negates */
SEXP _nlmixr2bayes_popEval(SEXP thetaS) {
  int ntheta, rc, i;
  double value = NA_REAL;
  SEXP val, grd, ret, nm;
  if (TYPEOF(thetaS) != REALSXP) Rf_error("'theta' must be a numeric vector");
  ntheta = (int) Rf_xlength(thetaS);
  grd = PROTECT(Rf_allocVector(REALSXP, ntheta));
  for (i = 0; i < ntheta; i++) REAL(grd)[i] = 0.0;
  rc = nlmixr2bayes_pop_eval(REAL(thetaS), ntheta, &value, REAL(grd));
  if (rc < 0) {
    UNPROTECT(1);
    Rf_error("nlmixr2bayes_pop_eval failed with status %d", rc);
  }
  val = PROTECT(Rf_ScalarReal(value));
  ret = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(ret, 0, val);
  SET_VECTOR_ELT(ret, 1, grd);
  SET_VECTOR_ELT(ret, 2, Rf_ScalarInteger(rc));
  nm = PROTECT(Rf_allocVector(STRSXP, 3));
  SET_STRING_ELT(nm, 0, Rf_mkChar("value"));
  SET_STRING_ELT(nm, 1, Rf_mkChar("grad"));
  SET_STRING_ELT(nm, 2, Rf_mkChar("nBad"));
  Rf_setAttrib(ret, R_NamesSymbol, nm);
  UNPROTECT(4);
  return ret;
}

/* THE symbol the compiled Stan model resolves via
 * R_GetCCallable("nlmixr2bayes", "nlmixr2bayes_cond_batch").  Signature frozen
 * by the injected header; never change it without bumping the header's guard.
 * eta/grad are nid x neta ROW-MAJOR (the nlmixr2est ABI). */
int nlmixr2bayes_cond_batch(const double *eta, int nid, int neta,
                           double *value, double *grad) {
  if (nlmixr2FoceiCondBatchP == NULL) return -100; /* table not installed */
  return nlmixr2FoceiCondBatchP(eta, nid, neta, nlmixr2bayesCores, value, grad);
}

/* ---- tier-2 state: theta base vector + mu-reference map -------------------
 * Stan samples only the model thetas; the loaded problem's parameter vector
 * is length npars = ntheta + <omega block>.  The omega tail stays at the
 * link-time values (the conditional likelihood does not depend on it; the
 * gradient assembly subtracts the same Omega^-1 term nlmixr2est added, so
 * the pair stays exact).  muRef[k] holds the 1-based theta index that
 * mu-references eta k (0 = none): for such a theta,
 * d/dtheta_p log p(y_i|eta_i) equals the eta_k gradient. */
static double *nlmixr2bayesThetaBase = NULL;
static int nlmixr2bayesNpars = 0;
static int *nlmixr2bayesMuRef = NULL;
static int nlmixr2bayesNeta = 0;
/* mu-referenced covariate coefficients: for mu_k = ... + theta_p * cov + eta_k
 * the chain rule gives d/dtheta_p = cov_i * d/deta_k.  covTheta[c] is the
 * 1-based theta index, covEta[c] the 0-based eta index, covVal the nid-per-
 * entry per-subject covariate values (entry-major: covVal[c*nid + i]). */
static int nlmixr2bayesNCov = 0;
static int nlmixr2bayesCovNid = 0;
static int *nlmixr2bayesCovTheta = NULL;
static int *nlmixr2bayesCovEta = NULL;
static double *nlmixr2bayesCovVal = NULL;

SEXP _nlmixr2bayes_setThetaBase(SEXP thetaS) {
  int n = (int) Rf_xlength(thetaS), i;
  if (TYPEOF(thetaS) != REALSXP || n < 1) Rf_error("'theta' must be a numeric vector");
  if (nlmixr2bayesThetaBase != NULL) { R_Free(nlmixr2bayesThetaBase); nlmixr2bayesThetaBase = NULL; }
  nlmixr2bayesThetaBase = R_Calloc((size_t) n, double);
  for (i = 0; i < n; i++) nlmixr2bayesThetaBase[i] = REAL(thetaS)[i];
  nlmixr2bayesNpars = n;
  return Rf_ScalarInteger(n);
}

SEXP _nlmixr2bayes_setMuRef(SEXP idxS) {
  int n = (int) Rf_xlength(idxS), i;
  if (TYPEOF(idxS) != INTSXP) Rf_error("'muRef' must be an integer vector");
  if (nlmixr2bayesMuRef != NULL) { R_Free(nlmixr2bayesMuRef); nlmixr2bayesMuRef = NULL; }
  nlmixr2bayesMuRef = R_Calloc((size_t) (n > 0 ? n : 1), int);
  for (i = 0; i < n; i++) nlmixr2bayesMuRef[i] = INTEGER(idxS)[i];
  nlmixr2bayesNeta = n;
  return Rf_ScalarInteger(n);
}

/* thetaIdxS: 1-based theta index per entry; etaIdxS: 0-based eta index per
 * entry; covValS: nid x nCov (col-major, which IS entry-major here) */
SEXP _nlmixr2bayes_setMuRefCov(SEXP thetaIdxS, SEXP etaIdxS, SEXP covValS) {
  int nc = (int) Rf_xlength(thetaIdxS), nid, i;
  SEXP dim;
  if (TYPEOF(thetaIdxS) != INTSXP || TYPEOF(etaIdxS) != INTSXP ||
      (int) Rf_xlength(etaIdxS) != nc) {
    Rf_error("'thetaIdx'/'etaIdx' must be integer vectors of equal length");
  }
  dim = Rf_getAttrib(covValS, R_DimSymbol);
  if (TYPEOF(covValS) != REALSXP || Rf_length(dim) != 2 ||
      INTEGER(dim)[1] != nc) {
    Rf_error("'covVal' must be a numeric nid x length(thetaIdx) matrix");
  }
  nid = INTEGER(dim)[0];
  if (nlmixr2bayesCovTheta != NULL) { R_Free(nlmixr2bayesCovTheta); nlmixr2bayesCovTheta = NULL; }
  if (nlmixr2bayesCovEta != NULL) { R_Free(nlmixr2bayesCovEta); nlmixr2bayesCovEta = NULL; }
  if (nlmixr2bayesCovVal != NULL) { R_Free(nlmixr2bayesCovVal); nlmixr2bayesCovVal = NULL; }
  nlmixr2bayesNCov = 0;
  nlmixr2bayesCovNid = 0;
  if (nc > 0) {
    nlmixr2bayesCovTheta = R_Calloc((size_t) nc, int);
    nlmixr2bayesCovEta = R_Calloc((size_t) nc, int);
    nlmixr2bayesCovVal = R_Calloc((size_t) nc * nid, double);
    for (i = 0; i < nc; i++) {
      nlmixr2bayesCovTheta[i] = INTEGER(thetaIdxS)[i];
      nlmixr2bayesCovEta[i] = INTEGER(etaIdxS)[i];
    }
    for (i = 0; i < nc * nid; i++) nlmixr2bayesCovVal[i] = REAL(covValS)[i];
    nlmixr2bayesNCov = nc;
    nlmixr2bayesCovNid = nid;
  }
  return Rf_ScalarInteger(nc);
}

SEXP _nlmixr2bayes_clearThetaBase(void) {
  if (nlmixr2bayesThetaBase != NULL) { R_Free(nlmixr2bayesThetaBase); nlmixr2bayesThetaBase = NULL; }
  if (nlmixr2bayesMuRef != NULL) { R_Free(nlmixr2bayesMuRef); nlmixr2bayesMuRef = NULL; }
  if (nlmixr2bayesCovTheta != NULL) { R_Free(nlmixr2bayesCovTheta); nlmixr2bayesCovTheta = NULL; }
  if (nlmixr2bayesCovEta != NULL) { R_Free(nlmixr2bayesCovEta); nlmixr2bayesCovEta = NULL; }
  if (nlmixr2bayesCovVal != NULL) { R_Free(nlmixr2bayesCovVal); nlmixr2bayesCovVal = NULL; }
  nlmixr2bayesNpars = 0;
  nlmixr2bayesNeta = 0;
  nlmixr2bayesNCov = 0;
  nlmixr2bayesCovNid = 0;
  return R_NilValue;
}

/* Tier 2: write theta (natural scale, first ntheta slots of the base
 * vector), then value + d/d(eta) + d/d(theta) of the conditional at the
 * supplied etas.  gradTheta is nid x ntheta row-major: the sensitivity
 * columns come from nlmixr2est's forward sensitivities, the mu-referenced
 * columns from the eta gradient via the identity above.  Returns >=0
 * (number of non-finite subjects) or <0 on hard failure; -101 when the
 * tier-2 state was never installed, -102 on a theta-set failure, -103 when
 * the theta-sensitivity model is not wired but non-mu thetas exist. */
int nlmixr2bayes_cond_batch_theta(const double *theta, int ntheta,
                                 const double *eta, int nid, int neta,
                                 double *value, double *gradEta,
                                 double *gradTheta) {
  int rc, rcT, i, k, p;
  size_t j;
  if (nlmixr2FoceiCondBatchP == NULL) return -100;
  if (nlmixr2bayesThetaBase == NULL || nlmixr2bayesMuRef == NULL ||
      ntheta > nlmixr2bayesNpars || neta != nlmixr2bayesNeta) return -101;
  for (j = 0; j < (size_t) ntheta; j++) nlmixr2bayesThetaBase[j] = theta[j];
  rc = nlmixr2FoceiSetThetaP(nlmixr2bayesThetaBase, nlmixr2bayesNpars);
  if (rc != 0) return -102;
  /* #958 fused entry: value + d/d(eta) + d/d(theta) from ONE combined-model
   * solve per subject.  Self-negotiating: -5 (not a combined build) and -4
   * (sensitivities not wired) fall back to the two-call path; other
   * negative codes are hard failures either way. */
  rc = -5;
  if (nlmixr2FoceiCondBatchThetaGradP != NULL) {
    rc = nlmixr2FoceiCondBatchThetaGradP(eta, nid, neta, nlmixr2bayesCores,
                                         value, gradEta, gradTheta);
    if (rc < 0 && rc != -5 && rc != -4) return rc;
  }
  if (rc < 0) {
    rc = nlmixr2FoceiCondBatchP(eta, nid, neta, nlmixr2bayesCores, value,
                                gradEta);
    if (rc < 0) return rc;
    /* forward-sensitivity theta columns (zero-filled by the callee); -4 =
     * model not wired, tolerable only when no sensitivity thetas exist */
    rcT = nlmixr2FoceiCondThetaGradP(eta, nid, neta, nlmixr2bayesCores,
                                     gradTheta);
    if (rcT == -4) {
      if (nlmixr2FoceiThetaSensIdxP(NULL, 0) != 0) return -103;
      for (j = 0; j < (size_t) nid * ntheta; j++) gradTheta[j] = 0.0;
    } else if (rcT < 0) {
      return rcT;
    }
  }
  /* mu-referenced columns: d/dtheta_p = d/deta_k of the conditional */
  for (i = 0; i < nid; i++) {
    for (k = 0; k < neta; k++) {
      p = nlmixr2bayesMuRef[k];
      if (p >= 1 && p <= ntheta) {
        gradTheta[(size_t) i * ntheta + (p - 1)] +=
          gradEta[(size_t) i * neta + k];
      }
    }
  }
  /* mu-referenced covariate coefficients: d/dtheta_p = cov_i * d/deta_k */
  if (nlmixr2bayesNCov > 0) {
    int c;
    if (nlmixr2bayesCovNid != nid) return -104;
    for (c = 0; c < nlmixr2bayesNCov; c++) {
      p = nlmixr2bayesCovTheta[c];
      k = nlmixr2bayesCovEta[c];
      if (p < 1 || p > ntheta || k < 0 || k >= neta) return -105;
      for (i = 0; i < nid; i++) {
        gradTheta[(size_t) i * ntheta + (p - 1)] +=
          nlmixr2bayesCovVal[(size_t) c * nid + i] *
          gradEta[(size_t) i * neta + k];
      }
    }
  }
  return rc;
}

/* R-callable mirror (col-major in/out) for tests and the FD oracle */
SEXP _nlmixr2bayes_condBatchTheta(SEXP thetaS, SEXP etaS) {
  int nid, neta, ntheta, i, j, rc;
  double *eta, *value, *gradEta, *gradTheta, *in;
  SEXP dim, val, ge, gt, ret, nm;
  dim = Rf_getAttrib(etaS, R_DimSymbol);
  if (TYPEOF(etaS) != REALSXP || Rf_length(dim) != 2) Rf_error("'eta' must be a numeric matrix");
  if (TYPEOF(thetaS) != REALSXP) Rf_error("'theta' must be a numeric vector");
  nid = INTEGER(dim)[0];
  neta = INTEGER(dim)[1];
  ntheta = (int) Rf_xlength(thetaS);
  in = REAL(etaS);
  eta = (double *) R_alloc((size_t) nid * neta, sizeof(double));
  value = (double *) R_alloc((size_t) nid, sizeof(double));
  gradEta = (double *) R_alloc((size_t) nid * neta, sizeof(double));
  gradTheta = (double *) R_alloc((size_t) nid * ntheta, sizeof(double));
  for (i = 0; i < nid; i++) {
    for (j = 0; j < neta; j++) eta[(size_t) i * neta + j] = in[i + (size_t) j * nid];
  }
  rc = nlmixr2bayes_cond_batch_theta(REAL(thetaS), ntheta, eta, nid, neta,
                                    value, gradEta, gradTheta);
  if (rc < 0) Rf_error("nlmixr2bayes_cond_batch_theta failed with status %d", rc);
  val = PROTECT(Rf_allocVector(REALSXP, nid));
  ge = PROTECT(Rf_allocMatrix(REALSXP, nid, neta));
  gt = PROTECT(Rf_allocMatrix(REALSXP, nid, ntheta));
  for (i = 0; i < nid; i++) {
    REAL(val)[i] = value[i];
    for (j = 0; j < neta; j++) REAL(ge)[i + (size_t) j * nid] = gradEta[(size_t) i * neta + j];
    for (j = 0; j < ntheta; j++) REAL(gt)[i + (size_t) j * nid] = gradTheta[(size_t) i * ntheta + j];
  }
  ret = PROTECT(Rf_allocVector(VECSXP, 4));
  SET_VECTOR_ELT(ret, 0, val);
  SET_VECTOR_ELT(ret, 1, ge);
  SET_VECTOR_ELT(ret, 2, gt);
  SET_VECTOR_ELT(ret, 3, Rf_ScalarInteger(rc));
  nm = PROTECT(Rf_allocVector(STRSXP, 4));
  SET_STRING_ELT(nm, 0, Rf_mkChar("value"));
  SET_STRING_ELT(nm, 1, Rf_mkChar("gradEta"));
  SET_STRING_ELT(nm, 2, Rf_mkChar("gradTheta"));
  SET_STRING_ELT(nm, 3, Rf_mkChar("nBad"));
  Rf_setAttrib(ret, R_NamesSymbol, nm);
  UNPROTECT(5);
  return ret;
}

/* ---- R-callable mirrors over the same entry points ------------------------
 * These are what the R layer and the tests use, so the exact C path the Stan
 * model will take is validated (finite differences, foceiLikRun agreement)
 * without Stan in the picture.  R matrices are COLUMN-major; the ABI is
 * row-major, so the wrappers transpose on the way in and out. */

SEXP _nlmixr2bayes_dims(void) {
  SEXP ret, nm;
  int nid = 0, neta = 0, ntheta = 0, npars = 0, flags = 0, rc;
  if (nlmixr2FoceiDimsP == NULL) Rf_error("nlmixr2est FOCEi pointer table not installed");
  rc = nlmixr2FoceiDimsP(&nid, &neta, &ntheta, &npars, &flags);
  ret = PROTECT(Rf_allocVector(INTSXP, 6));
  INTEGER(ret)[0] = rc;
  INTEGER(ret)[1] = nid;
  INTEGER(ret)[2] = neta;
  INTEGER(ret)[3] = ntheta;
  INTEGER(ret)[4] = npars;
  INTEGER(ret)[5] = flags;
  nm = PROTECT(Rf_allocVector(STRSXP, 6));
  SET_STRING_ELT(nm, 0, Rf_mkChar("status"));
  SET_STRING_ELT(nm, 1, Rf_mkChar("nid"));
  SET_STRING_ELT(nm, 2, Rf_mkChar("neta"));
  SET_STRING_ELT(nm, 3, Rf_mkChar("ntheta"));
  SET_STRING_ELT(nm, 4, Rf_mkChar("npars"));
  SET_STRING_ELT(nm, 5, Rf_mkChar("flags"));
  Rf_setAttrib(ret, R_NamesSymbol, nm);
  UNPROTECT(2);
  return ret;
}

SEXP _nlmixr2bayes_setTheta(SEXP thetaS) {
  if (nlmixr2FoceiSetThetaP == NULL) Rf_error("nlmixr2est FOCEi pointer table not installed");
  if (TYPEOF(thetaS) != REALSXP) Rf_error("'theta' must be a numeric vector");
  return Rf_ScalarInteger(nlmixr2FoceiSetThetaP(REAL(thetaS), (int)Rf_xlength(thetaS)));
}

SEXP _nlmixr2bayes_setOmegaInv(SEXP m) {
  int n, i, j, rc;
  double *oi, *in;
  SEXP dim;
  if (nlmixr2FoceiSetOmegaInvP == NULL) Rf_error("nlmixr2est FOCEi pointer table not installed");
  dim = Rf_getAttrib(m, R_DimSymbol);
  if (TYPEOF(m) != REALSXP || Rf_length(dim) != 2 ||
      INTEGER(dim)[0] != INTEGER(dim)[1]) {
    Rf_error("'omegaInv' must be a square numeric matrix");
  }
  n = INTEGER(dim)[0];
  in = REAL(m);
  oi = (double *) R_alloc((size_t) n * n, sizeof(double));
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      oi[(size_t) i * n + j] = in[i + (size_t) j * n]; /* col- -> row-major */
    }
  }
  rc = nlmixr2FoceiSetOmegaInvP(oi, n);
  return Rf_ScalarInteger(rc);
}

/* value + d/d(eta) of the conditional log-likelihood; etaS is an R (col-major)
 * nid x neta matrix.  Returns list(value=, grad=, nBad=); errors on a hard
 * failure code so R callers cannot mistake it for a rejection. */
SEXP _nlmixr2bayes_condBatch(SEXP etaS) {
  int nid, neta, i, j, rc;
  double *eta, *value, *grad, *in, *gv, *gg;
  SEXP dim, val, grd, nbad, ret, nm;
  dim = Rf_getAttrib(etaS, R_DimSymbol);
  if (TYPEOF(etaS) != REALSXP || Rf_length(dim) != 2) {
    Rf_error("'eta' must be a numeric matrix");
  }
  nid = INTEGER(dim)[0];
  neta = INTEGER(dim)[1];
  in = REAL(etaS);
  eta = (double *) R_alloc((size_t) nid * neta, sizeof(double));
  value = (double *) R_alloc((size_t) nid, sizeof(double));
  grad = (double *) R_alloc((size_t) nid * neta, sizeof(double));
  for (i = 0; i < nid; i++) {
    for (j = 0; j < neta; j++) {
      eta[(size_t) i * neta + j] = in[i + (size_t) j * nid];
    }
  }
  rc = nlmixr2bayes_cond_batch(eta, nid, neta, value, grad);
  if (rc < 0) Rf_error("nlmixr2bayes_cond_batch failed with status %d", rc);
  val = PROTECT(Rf_allocVector(REALSXP, nid));
  grd = PROTECT(Rf_allocMatrix(REALSXP, nid, neta));
  gv = REAL(val);
  gg = REAL(grd);
  for (i = 0; i < nid; i++) {
    gv[i] = value[i];
    for (j = 0; j < neta; j++) {
      gg[i + (size_t) j * nid] = grad[(size_t) i * neta + j];
    }
  }
  nbad = PROTECT(Rf_ScalarInteger(rc));
  ret = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(ret, 0, val);
  SET_VECTOR_ELT(ret, 1, grd);
  SET_VECTOR_ELT(ret, 2, nbad);
  nm = PROTECT(Rf_allocVector(STRSXP, 3));
  SET_STRING_ELT(nm, 0, Rf_mkChar("value"));
  SET_STRING_ELT(nm, 1, Rf_mkChar("grad"));
  SET_STRING_ELT(nm, 2, Rf_mkChar("nBad"));
  Rf_setAttrib(ret, R_NamesSymbol, nm);
  UNPROTECT(5);
  return ret;
}
