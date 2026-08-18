// Path A: re-drive a solve that rxode2 already built, instead of paying
// rxSolve()'s ~2.3 ms of R-level setup and teardown on every HMC gradient.
//
// rxSolve.default calls rxSolveFree() at its START, to clear the *previous*
// solve, so the rx_solve structure is still live once rxSolve() returns.  That
// is what makes this possible -- and also what makes it fragile: any other
// rxSolve() in the session replaces that structure.  Hence the invariant check
// before every fast solve, and the fall back to the slow R path when it fails.

// Rcpp and the C++ standard headers must come before rxode2ptr.h: that header
// pulls in R's own without R_NO_REMAP, and the resulting `length`/`error`
// macros break libc++.
#include <Rcpp.h>

#include <cstdio>
#include <cstring>
#include <vector>

#include <time.h>  // rxode2ptr.h uses clock_t without including it
#include <rxode2ptr.h>

#include "rxstan.h"

namespace {

struct FastState {
  bool active = false;
  int handle = -1;

  rx_solve *rx = NULL;
  double *solve0 = NULL;  // ind 0's solve buffer, part of the identity check

  int nsub = 0;
  int npars = 0;
  int neq = 0;
  int nobs = 0;

  std::vector<double> params;  // canonical full parameter vector, per subject
  std::vector<int> sensIdx;    // par_ptr index of each block parameter
  std::vector<int> outIdx;     // state index of each requested output
  std::vector<int> sensState;  // state index of d out[k] / d p_block[j], k-major
  std::vector<int> blockOf;    // parameter block each subject draws from
  int nBlocks = 1;
  bool quiet = false;

  // Observation time-rows, flattened subject-major to match rxSolve()'s output.
  std::vector<int> obsSub;
  std::vector<int> obsRow;
};

FastState &fast() {
  static FastState s;
  return s;
}

// Is the structure we cached still the one rxode2 holds?  Catches the common
// failure -- someone else called rxSolve, freeing and rebuilding it -- cheaply.
bool still_ours(const FastState &f) {
  if (!f.active || f.rx == NULL) return false;
  rx_solve *rx = getRxSolve_();
  if (rx != f.rx) return false;
  if (getRxNsub(rx) != f.nsub) return false;
  if (getRxNpars(rx) != f.npars) return false;

  rx_solving_options *op = getSolvingOptions(rx);
  if (getOpNeq(op) != f.neq) return false;
  if (getRxNobs(rx) != f.nobs) return false;

  rx_solving_options_ind *ind = getSolvingOptionsInd(rx, 0);
  if (getIndSolve(ind) != f.solve0) return false;
  return true;
}

}  // namespace

extern "C" {

int rxs_fast_available(int handle) {
  FastState &f = fast();
  return (f.active && f.handle == handle && still_ours(f)) ? 1 : 0;
}

void rxs_fast_invalidate(void) { fast().active = false; }

int rxs_fast_solve(int handle, const double *p, int np, double *y, double *dydp,
                   int ny, char *errbuf, int errlen) {
  FastState &f = fast();
  if (!rxs_fast_available(handle)) {
    if (errlen > 0) std::snprintf(errbuf, errlen, "fast path unavailable");
    return 1;
  }
  const int nBlock = (int)f.sensIdx.size();
  const int nOut = (int)f.outIdx.size();
  const int nRow = (int)f.obsRow.size();

  if (np != nBlock * f.nBlocks || ny != nRow * nOut) {
    if (errlen > 0) std::snprintf(errbuf, errlen, "fast path dimension mismatch");
    return 2;
  }

  rx_solve *rx = f.rx;
  rx_solving_options *op = getSolvingOptions(rx);

  // Re-assert every parameter, not just the differentiated ones, so a foreign
  // solve that happened to leave the shape intact cannot leak its values in.
  // Reset each subject's solve counter too: after a failed solve rxode2 leaves
  // ind dirty, and par_solve() on its own will keep returning that NaN state.
  for (int s = 0; s < f.nsub; ++s) {
    double *ps = &f.params[(size_t)s * f.npars];
    const double *pb = p + (size_t)f.blockOf[s] * nBlock;
    for (int j = 0; j < nBlock; ++j) ps[f.sensIdx[j]] = pb[j];

    rx_solving_options_ind *ind = getSolvingOptionsInd(rx, s);
    setIndSolve(ind, 0);
    for (int i = 0; i < f.npars; ++i) setIndParPtr(ind, i, ps[i]);
  }

  resetOpBadSolve(op);
  if (f.quiet) rxSetSilentErr(1);
  par_solve(rx);
  if (f.quiet) rxSetSilentErr(0);
  if (hasOpBadSolve(op)) {
    if (errlen > 0) std::snprintf(errbuf, errlen, "rxode2 reported a bad solve");
    return 3;
  }

  for (int r = 0; r < nRow; ++r) {
    rx_solving_options_ind *ind = getSolvingOptionsInd(rx, f.obsSub[r]);
    const double *row = getOpIndSolve(op, ind, f.obsRow[r]);
    for (int k = 0; k < nOut; ++k) {
      const int i = k * nRow + r;
      const double v = row[f.outIdx[k]];
      if (!R_finite(v)) {
        if (errlen > 0) {
          std::snprintf(errbuf, errlen, "non-finite solve at output %d", i);
        }
        return 4;
      }
      y[i] = v;
      for (int j = 0; j < nBlock; ++j) {
        const double g = row[f.sensState[k * nBlock + j]];
        if (!R_finite(g)) {
          if (errlen > 0) {
            std::snprintf(errbuf, errlen,
                          "non-finite sensitivity at output %d, parameter %d", i,
                          j);
          }
          return 5;
        }
        dydp[i + ny * j] = g;
      }
    }
  }
  return 0;
}

// Captures the live solve immediately after rxsRegister()'s probe rxSolve().
SEXP rxstanFastSetup(SEXP handleSXP, SEXP sensIdxSXP, SEXP outIdxSXP,
                     SEXP sensStateSXP, SEXP nobsSXP, SEXP blockOfSXP,
                     SEXP quietSXP) {
  FastState &f = fast();
  f.active = false;
  f.quiet = (Rf_asLogical(quietSXP) == TRUE);

  rx_solve *rx = getRxSolve_();
  if (rx == NULL) return Rf_ScalarLogical(FALSE);
  if (getRxNsub(rx) < 1) return Rf_ScalarLogical(FALSE);

  rx_solving_options *op = getSolvingOptions(rx);
  f.rx = rx;
  f.nsub = getRxNsub(rx);
  f.npars = getRxNpars(rx);
  f.neq = getOpNeq(op);
  f.nobs = getRxNobs(rx);
  f.solve0 = getIndSolve(getSolvingOptionsInd(rx, 0));

  if (Rf_length(blockOfSXP) != f.nsub) return Rf_ScalarLogical(FALSE);
  f.blockOf.assign(INTEGER(blockOfSXP), INTEGER(blockOfSXP) + f.nsub);
  f.nBlocks = 0;
  for (int s = 0; s < f.nsub; ++s) {
    if (f.blockOf[s] < 0) return Rf_ScalarLogical(FALSE);
    if (f.blockOf[s] + 1 > f.nBlocks) f.nBlocks = f.blockOf[s] + 1;
  }

  // Each subject keeps its own copy: only the block parameters differ, but the
  // fast path writes the whole vector back through setIndParPtr anyway.
  f.params.assign((size_t)f.nsub * f.npars, 0.0);
  for (int s = 0; s < f.nsub; ++s) {
    rx_solving_options_ind *ind = getSolvingOptionsInd(rx, s);
    for (int i = 0; i < f.npars; ++i) {
      f.params[(size_t)s * f.npars + i] = getIndParPtr(ind, i);
    }
  }

  f.sensIdx.assign(INTEGER(sensIdxSXP),
                   INTEGER(sensIdxSXP) + Rf_length(sensIdxSXP));
  f.outIdx.assign(INTEGER(outIdxSXP),
                  INTEGER(outIdxSXP) + Rf_length(outIdxSXP));
  f.sensState.assign(INTEGER(sensStateSXP),
                     INTEGER(sensStateSXP) + Rf_length(sensStateSXP));

  for (size_t i = 0; i < f.sensIdx.size(); ++i) {
    if (f.sensIdx[i] < 0 || f.sensIdx[i] >= f.npars) {
      return Rf_ScalarLogical(FALSE);
    }
  }
  for (size_t i = 0; i < f.outIdx.size(); ++i) {
    if (f.outIdx[i] < 0 || f.outIdx[i] >= f.neq) return Rf_ScalarLogical(FALSE);
  }
  for (size_t i = 0; i < f.sensState.size(); ++i) {
    if (f.sensState[i] < 0 || f.sensState[i] >= f.neq) {
      return Rf_ScalarLogical(FALSE);
    }
  }

  // Observation rows, subject-major, matching how rxSolve() stacks its output.
  f.obsSub.clear();
  f.obsRow.clear();
  for (int s = 0; s < f.nsub; ++s) {
    rx_solving_options_ind *ind = getSolvingOptionsInd(rx, s);
    const int nt = getIndNallTimes(ind);
    for (int j = 0; j < nt; ++j) {
      if (getIndEvid(ind, getIndIx(ind, j)) == 0) {
        f.obsSub.push_back(s);
        f.obsRow.push_back(j);
      }
    }
  }

  // If our notion of "observation" disagrees with rxSolve's row count then the
  // layout assumption is wrong and the fast path must not be used at all.
  if ((int)f.obsRow.size() != Rf_asInteger(nobsSXP)) {
    return Rf_ScalarLogical(FALSE);
  }

  f.handle = Rf_asInteger(handleSXP);
  f.active = true;
  return Rf_ScalarLogical(TRUE);
}

SEXP rxstanFastInvalidate(void) {
  rxs_fast_invalidate();
  return R_NilValue;
}

SEXP rxstanFastAvailable(SEXP handleSXP) {
  return Rf_ScalarLogical(rxs_fast_available(Rf_asInteger(handleSXP)));
}

// Toggles rxode2's own solver diagnostics.  Only the repeated per-solve
// chatter is affected: a failed solve is still reported through the return
// code and still becomes a domain_error, and rxsSolveStats() keeps counting.
SEXP rxstanSetSilent(SEXP silentSXP) {
  rxSetSilentErr(Rf_asLogical(silentSXP) == TRUE ? 1 : 0);
  return R_NilValue;
}

// --- layout probes, kept because the fast path's correctness rests on them ---

SEXP rxstanProbeSolve(void) {
  rx_solve *rx = getRxSolve_();
  if (rx == NULL) return R_NilValue;

  rx_solving_options *op = getSolvingOptions(rx);
  const int nsub = getRxNsub(rx);

  const char *nms[] = {"nsub", "npars", "neq", "nall", "nobs", "nlhs",
                       "nAllTimes0", "nDoses0", "nevid2_0"};

  SEXP out = PROTECT(Rf_allocVector(INTSXP, 9));
  INTEGER(out)[0] = nsub;
  INTEGER(out)[1] = getRxNpars(rx);
  INTEGER(out)[2] = getOpNeq(op);
  INTEGER(out)[3] = getRxNall(rx);
  INTEGER(out)[4] = getRxNobs(rx);
  INTEGER(out)[5] = getOpNlhs(op);

  if (nsub > 0) {
    rx_solving_options_ind *ind = getSolvingOptionsInd(rx, 0);
    INTEGER(out)[6] = getIndNallTimes(ind);
    INTEGER(out)[7] = getIndNdoses(ind);
    INTEGER(out)[8] = getIndNevid2(ind);
  } else {
    INTEGER(out)[6] = INTEGER(out)[7] = INTEGER(out)[8] = NA_INTEGER;
  }

  SEXP nm = PROTECT(Rf_allocVector(STRSXP, 9));
  for (int i = 0; i < 9; ++i) SET_STRING_ELT(nm, i, Rf_mkChar(nms[i]));
  Rf_setAttrib(out, R_NamesSymbol, nm);
  UNPROTECT(2);
  return out;
}

SEXP rxstanProbeParams(void) {
  rx_solve *rx = getRxSolve_();
  if (rx == NULL || getRxNsub(rx) < 1) return R_NilValue;

  const int npars = getRxNpars(rx);
  rx_solving_options_ind *ind = getSolvingOptionsInd(rx, 0);

  SEXP out = PROTECT(Rf_allocVector(REALSXP, npars));
  for (int i = 0; i < npars; ++i) REAL(out)[i] = getIndParPtr(ind, i);
  UNPROTECT(1);
  return out;
}

SEXP rxstanProbeStates(void) {
  rx_solve *rx = getRxSolve_();
  if (rx == NULL || getRxNsub(rx) < 1) return R_NilValue;

  rx_solving_options *op = getSolvingOptions(rx);
  const int neq = getOpNeq(op);
  rx_solving_options_ind *ind = getSolvingOptionsInd(rx, 0);
  const int nt = getIndNallTimes(ind);

  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, nt, neq + 2));
  double *o = REAL(out);
  for (int j = 0; j < nt; ++j) {
    double *yj = getOpIndSolve(op, ind, j);
    for (int k = 0; k < neq; ++k) o[j + nt * k] = yj[k];
    o[j + nt * neq] = getTime(getIndIx(ind, j), ind);
    o[j + nt * (neq + 1)] = (double)getIndEvid(ind, getIndIx(ind, j));
  }
  UNPROTECT(1);
  return out;
}

}  // extern "C"
