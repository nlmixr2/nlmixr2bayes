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
