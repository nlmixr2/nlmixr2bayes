/* Consume nlmixr2est's FOCEi conditional-likelihood C API through its
 * external-pointer table (inst/include/nlmixr2estFoceiPtr.h, installed at
 * .onLoad from nlmixr2est::.nlmixr2estFoceiPtrs()), and re-export ONE plain-C
 * symbol via R_RegisterCCallable for the injected Stan header
 * (inst/include/nlmixr2stan_lp.hpp) to find with R_GetCCallable.
 *
 * Deliberately NO Rcpp, NO Armadillo, NO Eigen, NO StanHeaders here: the
 * external function is compiled into rstan's own translation unit, and
 * everything crossing a shared-object boundary is a double pointer or an int,
 * so no C++ type ever appears in more than one DSO (see the design notes in
 * nlmixr2/nlmixr2stan#1). */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <nlmixr2estFoceiPtr.h>

/* instantiates the 7 function-pointer globals + SEXP iniNlmixr2estFocei(SEXP) */
iniNlmixr2estFoceiGlobals

/* Subject-parallel thread count for the batch entries, set from R before
 * sampling (rxode2/nlmixr2est parallelize over subjects inside each call;
 * Stan itself stays single-threaded). */
static int nlmixr2stanCores = 1;

SEXP _nlmixr2stan_setCores(SEXP n) {
  int c = Rf_asInteger(n);
  if (c == NA_INTEGER || c < 1) c = 1;
  nlmixr2stanCores = c;
  return Rf_ScalarInteger(nlmixr2stanCores);
}

SEXP _nlmixr2stan_apiVersion(void) {
  if (nlmixr2FoceiApiVersionP == NULL) return Rf_ScalarInteger(-1);
  return Rf_ScalarInteger(nlmixr2FoceiApiVersionP());
}

/* THE symbol the compiled Stan model resolves via
 * R_GetCCallable("nlmixr2stan", "nlmixr2stan_cond_batch").  Signature frozen
 * by the injected header; never change it without bumping the header's guard.
 * eta/grad are nid x neta ROW-MAJOR (the nlmixr2est ABI). */
int nlmixr2stan_cond_batch(const double *eta, int nid, int neta,
                           double *value, double *grad) {
  if (nlmixr2FoceiCondBatchP == NULL) return -100; /* table not installed */
  return nlmixr2FoceiCondBatchP(eta, nid, neta, nlmixr2stanCores, value, grad);
}

/* ---- tier-2 state: theta base vector + mu-reference map -------------------
 * Stan samples only the model thetas; the loaded problem's parameter vector
 * is length npars = ntheta + <omega block>.  The omega tail stays at the
 * link-time values (the conditional likelihood does not depend on it; the
 * gradient assembly subtracts the same Omega^-1 term nlmixr2est added, so
 * the pair stays exact).  muRef[k] holds the 1-based theta index that
 * mu-references eta k (0 = none): for such a theta,
 * d/dtheta_p log p(y_i|eta_i) equals the eta_k gradient. */
static double *nlmixr2stanThetaBase = NULL;
static int nlmixr2stanNpars = 0;
static int *nlmixr2stanMuRef = NULL;
static int nlmixr2stanNeta = 0;

SEXP _nlmixr2stan_setThetaBase(SEXP thetaS) {
  int n = (int) Rf_xlength(thetaS), i;
  if (TYPEOF(thetaS) != REALSXP || n < 1) Rf_error("'theta' must be a numeric vector");
  if (nlmixr2stanThetaBase != NULL) { R_Free(nlmixr2stanThetaBase); nlmixr2stanThetaBase = NULL; }
  nlmixr2stanThetaBase = R_Calloc((size_t) n, double);
  for (i = 0; i < n; i++) nlmixr2stanThetaBase[i] = REAL(thetaS)[i];
  nlmixr2stanNpars = n;
  return Rf_ScalarInteger(n);
}

SEXP _nlmixr2stan_setMuRef(SEXP idxS) {
  int n = (int) Rf_xlength(idxS), i;
  if (TYPEOF(idxS) != INTSXP) Rf_error("'muRef' must be an integer vector");
  if (nlmixr2stanMuRef != NULL) { R_Free(nlmixr2stanMuRef); nlmixr2stanMuRef = NULL; }
  nlmixr2stanMuRef = R_Calloc((size_t) (n > 0 ? n : 1), int);
  for (i = 0; i < n; i++) nlmixr2stanMuRef[i] = INTEGER(idxS)[i];
  nlmixr2stanNeta = n;
  return Rf_ScalarInteger(n);
}

SEXP _nlmixr2stan_clearThetaBase(void) {
  if (nlmixr2stanThetaBase != NULL) { R_Free(nlmixr2stanThetaBase); nlmixr2stanThetaBase = NULL; }
  if (nlmixr2stanMuRef != NULL) { R_Free(nlmixr2stanMuRef); nlmixr2stanMuRef = NULL; }
  nlmixr2stanNpars = 0;
  nlmixr2stanNeta = 0;
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
int nlmixr2stan_cond_batch_theta(const double *theta, int ntheta,
                                 const double *eta, int nid, int neta,
                                 double *value, double *gradEta,
                                 double *gradTheta) {
  int rc, rcT, i, k, p;
  size_t j;
  if (nlmixr2FoceiCondBatchP == NULL) return -100;
  if (nlmixr2stanThetaBase == NULL || nlmixr2stanMuRef == NULL ||
      ntheta > nlmixr2stanNpars || neta != nlmixr2stanNeta) return -101;
  for (j = 0; j < (size_t) ntheta; j++) nlmixr2stanThetaBase[j] = theta[j];
  rc = nlmixr2FoceiSetThetaP(nlmixr2stanThetaBase, nlmixr2stanNpars);
  if (rc != 0) return -102;
  rc = nlmixr2FoceiCondBatchP(eta, nid, neta, nlmixr2stanCores, value, gradEta);
  if (rc < 0) return rc;
  /* forward-sensitivity theta columns (zero-filled by the callee); -4 =
   * model not wired, tolerable only when no sensitivity thetas exist */
  rcT = nlmixr2FoceiCondThetaGradP(eta, nid, neta, nlmixr2stanCores, gradTheta);
  if (rcT == -4) {
    if (nlmixr2FoceiThetaSensIdxP(NULL, 0) != 0) return -103;
    for (j = 0; j < (size_t) nid * ntheta; j++) gradTheta[j] = 0.0;
  } else if (rcT < 0) {
    return rcT;
  }
  /* mu-referenced columns: d/dtheta_p = d/deta_k of the conditional */
  for (i = 0; i < nid; i++) {
    for (k = 0; k < neta; k++) {
      p = nlmixr2stanMuRef[k];
      if (p >= 1 && p <= ntheta) {
        gradTheta[(size_t) i * ntheta + (p - 1)] +=
          gradEta[(size_t) i * neta + k];
      }
    }
  }
  return rc;
}

/* R-callable mirror (col-major in/out) for tests and the FD oracle */
SEXP _nlmixr2stan_condBatchTheta(SEXP thetaS, SEXP etaS) {
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
  rc = nlmixr2stan_cond_batch_theta(REAL(thetaS), ntheta, eta, nid, neta,
                                    value, gradEta, gradTheta);
  if (rc < 0) Rf_error("nlmixr2stan_cond_batch_theta failed with status %d", rc);
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

SEXP _nlmixr2stan_dims(void) {
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

SEXP _nlmixr2stan_setTheta(SEXP thetaS) {
  if (nlmixr2FoceiSetThetaP == NULL) Rf_error("nlmixr2est FOCEi pointer table not installed");
  if (TYPEOF(thetaS) != REALSXP) Rf_error("'theta' must be a numeric vector");
  return Rf_ScalarInteger(nlmixr2FoceiSetThetaP(REAL(thetaS), (int)Rf_xlength(thetaS)));
}

SEXP _nlmixr2stan_setOmegaInv(SEXP m) {
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
SEXP _nlmixr2stan_condBatch(SEXP etaS) {
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
  rc = nlmixr2stan_cond_batch(eta, nid, neta, value, grad);
  if (rc < 0) Rf_error("nlmixr2stan_cond_batch failed with status %d", rc);
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
