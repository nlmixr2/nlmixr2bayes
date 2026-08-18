// The no-longjmp bridge between Stan's C++ and rxode2's R-level solve.
//
// Everything Stan calls goes through rxs_solve_sens(), which evaluates
// rxstan:::.rxsSolveOne() under R_tryCatchError().  That guard is the whole
// point of this file: without it an Rf_error() inside rxode2 would longjmp
// straight out of Stan's log_prob(), skipping the destructors that release the
// autodiff arena.

#include "rxstan.h"

#include <cstdio>
#include <cstring>
#include <map>
#include <vector>

namespace {

struct Dims {
  int ny;
  int nBlock;   // parameters per block
  int nBlocks;  // 1 when parameters are shared, nsub for a population model
  std::vector<int> outBlock;  // length ny
  double nSolve;
  double nFail;
};

std::map<int, Dims>& dims_registry() {
  static std::map<int, Dims> reg;
  return reg;
}

// rxstan:::.rxsSolveOne, looked up once and kept alive for the session.
// Only ever called from solve_body(), i.e. under R_tryCatchError, so an
// unbound symbol surfaces as a caught R error rather than a longjmp.
SEXP solver_closure() {
  static SEXP fn = R_NilValue;
  if (fn == R_NilValue) {
    SEXP ns = PROTECT(R_FindNamespace(Rf_mkString("nlmixr2bayes")));
    SEXP found = PROTECT(Rf_eval(Rf_install(".rxsSolveOne"), ns));
    R_PreserveObject(found);
    fn = found;
    UNPROTECT(2);
  }
  return fn;
}

struct BodyData {
  int handle;
  const double* p;
  int np;
};

struct ErrData {
  char* buf;
  int len;
  bool failed;
};

SEXP solve_body(void* bdata) {
  BodyData* b = static_cast<BodyData*>(bdata);
  SEXP fn = solver_closure();

  SEXP hv = PROTECT(Rf_ScalarInteger(b->handle));
  SEXP pv = PROTECT(Rf_allocVector(REALSXP, b->np));
  std::memcpy(REAL(pv), b->p, sizeof(double) * static_cast<size_t>(b->np));
  SEXP call = PROTECT(Rf_lang3(fn, hv, pv));
  SEXP res = Rf_eval(call, R_GlobalEnv);
  UNPROTECT(3);
  return res;
}

SEXP solve_handler(SEXP cond, void* hdata) {
  ErrData* e = static_cast<ErrData*>(hdata);
  e->failed = true;
  if (e->len > 0) std::snprintf(e->buf, e->len, "unspecified rxode2 failure");

  // conditionMessage() can itself fail; a bare message is better than a crash.
  SEXP call = PROTECT(Rf_lang2(Rf_install("conditionMessage"), cond));
  SEXP msg = PROTECT(Rf_eval(call, R_BaseEnv));
  if (TYPEOF(msg) == STRSXP && Rf_length(msg) > 0 && e->len > 0) {
    std::snprintf(e->buf, e->len, "%s", CHAR(STRING_ELT(msg, 0)));
  }
  UNPROTECT(2);
  return R_NilValue;
}

}  // namespace

extern "C" {

int rxs_layout(int handle, int* ny, int* np, int* nBlock, int* nBlocks) {
  std::map<int, Dims>::const_iterator it = dims_registry().find(handle);
  if (it == dims_registry().end()) return 1;
  if (ny != NULL) *ny = it->second.ny;
  if (np != NULL) *np = it->second.nBlock * it->second.nBlocks;
  if (nBlock != NULL) *nBlock = it->second.nBlock;
  if (nBlocks != NULL) *nBlocks = it->second.nBlocks;
  return 0;
}

int rxs_dims(int handle, int* ny, int* np) {
  return rxs_layout(handle, ny, np, NULL, NULL);
}

int rxs_out_block(int handle, int* out, int n) {
  std::map<int, Dims>::const_iterator it = dims_registry().find(handle);
  if (it == dims_registry().end()) return 1;
  if (n != it->second.ny) return 2;
  for (int i = 0; i < n; ++i) out[i] = it->second.outBlock[i];
  return 0;
}

int rxs_stats(int handle, double* nSolve, double* nFail) {
  std::map<int, Dims>::const_iterator it = dims_registry().find(handle);
  if (it == dims_registry().end()) return 1;
  if (nSolve != NULL) *nSolve = it->second.nSolve;
  if (nFail != NULL) *nFail = it->second.nFail;
  return 0;
}

// The slow, always-correct route: evaluate .rxsSolveOne() in R.  It returns the
// value column followed by the block-diagonal Jacobian, ny x (nBlock + 1),
// which is already the shape rxode2's rx__sens_* columns come in.
static int solve_via_r(int handle, const double* p, int np, double* y,
                       double* dydp, int ny, int nBlock, char* errbuf,
                       int errlen) {
  BodyData b;
  b.handle = handle;
  b.p = p;
  b.np = np;

  ErrData e;
  e.buf = errbuf;
  e.len = errlen;
  e.failed = false;

  SEXP res = PROTECT(R_tryCatchError(solve_body, &b, solve_handler, &e));
  if (e.failed) {
    UNPROTECT(1);
    return 3;
  }

  if (TYPEOF(res) != REALSXP || Rf_length(res) != ny * (nBlock + 1)) {
    if (errlen > 0) {
      std::snprintf(errbuf, errlen,
                    ".rxsSolveOne returned %d doubles, expected %d",
                    static_cast<int>(Rf_length(res)), ny * (nBlock + 1));
    }
    UNPROTECT(1);
    return 4;
  }

  const double* m = REAL(res);
  for (int i = 0; i < ny; ++i) {
    if (!R_finite(m[i])) {
      if (errlen > 0) {
        std::snprintf(errbuf, errlen, "non-finite solve at output %d", i);
      }
      UNPROTECT(1);
      return 5;
    }
    y[i] = m[i];
  }
  for (int j = 0; j < nBlock; ++j) {
    for (int i = 0; i < ny; ++i) {
      const double g = m[i + ny * (j + 1)];
      if (!R_finite(g)) {
        if (errlen > 0) {
          std::snprintf(errbuf, errlen,
                        "non-finite sensitivity at output %d, parameter %d", i,
                        j);
        }
        UNPROTECT(1);
        return 6;
      }
      dydp[i + ny * j] = g;
    }
  }

  UNPROTECT(1);
  return 0;
}

static int check_dims(int handle, int np, int ny, int* nBlock, char* errbuf,
                      int errlen) {
  int rny = 0, rnp = 0;
  if (rxs_layout(handle, &rny, &rnp, nBlock, NULL) != 0) {
    if (errlen > 0) std::snprintf(errbuf, errlen, "unknown handle %d", handle);
    return 1;
  }
  if (rny != ny || rnp != np) {
    if (errlen > 0) {
      std::snprintf(errbuf, errlen,
                    "handle %d dimension mismatch: caller %dx%d, registry %dx%d",
                    handle, ny, np, rny, rnp);
    }
    return 2;
  }
  return 0;
}

int rxs_solve_sens(int handle, const double* p, int np, double* y, double* dydp,
                   int ny, char* errbuf, int errlen) {
  if (errlen > 0) errbuf[0] = '\0';

  int nBlock = 0;
  const int bad = check_dims(handle, np, ny, &nBlock, errbuf, errlen);
  if (bad != 0) return bad;

  Dims& d = dims_registry()[handle];
  d.nSolve += 1;

  // Path A when rxode2 still holds the solve we set up, otherwise rebuild it
  // through R.  ANY fast-path failure disarms it and retries slowly: a failed
  // solve leaves rxode2's per-subject state dirty, and only a full rxSolve()
  // reliably clears it.  A successful slow solve re-arms the fast path, so the
  // cost of a rejected HMC proposal is a couple of slow solves rather than a
  // silently poisoned cache.
  if (rxs_fast_available(handle)) {
    if (rxs_fast_solve(handle, p, np, y, dydp, ny, errbuf, errlen) == 0) {
      return 0;
    }
    rxs_fast_invalidate();
    if (errlen > 0) errbuf[0] = '\0';
  }

  const int rc =
      solve_via_r(handle, p, np, y, dydp, ny, nBlock, errbuf, errlen);
  if (rc != 0) d.nFail += 1;
  return rc;
}

SEXP rxstanStats(SEXP handleSXP) {
  double nSolve = 0, nFail = 0;
  if (rxs_stats(Rf_asInteger(handleSXP), &nSolve, &nFail) != 0) {
    Rf_error("unknown rxstan handle %d", Rf_asInteger(handleSXP));
  }
  SEXP out = PROTECT(Rf_allocVector(REALSXP, 2));
  REAL(out)[0] = nSolve;
  REAL(out)[1] = nFail;
  SEXP nm = PROTECT(Rf_allocVector(STRSXP, 2));
  SET_STRING_ELT(nm, 0, Rf_mkChar("solves"));
  SET_STRING_ELT(nm, 1, Rf_mkChar("failures"));
  Rf_setAttrib(out, R_NamesSymbol, nm);
  UNPROTECT(2);
  return out;
}

SEXP rxstanResetStats(SEXP handleSXP) {
  std::map<int, Dims>::iterator it = dims_registry().find(Rf_asInteger(handleSXP));
  if (it == dims_registry().end()) {
    Rf_error("unknown rxstan handle %d", Rf_asInteger(handleSXP));
  }
  it->second.nSolve = 0;
  it->second.nFail = 0;
  return R_NilValue;
}

SEXP rxstanSetDims(SEXP handleSXP, SEXP nySXP, SEXP nBlockSXP, SEXP nBlocksSXP,
                   SEXP outBlockSXP) {
  Dims d;
  d.ny = Rf_asInteger(nySXP);
  d.nBlock = Rf_asInteger(nBlockSXP);
  d.nBlocks = Rf_asInteger(nBlocksSXP);
  d.outBlock.assign(INTEGER(outBlockSXP),
                    INTEGER(outBlockSXP) + Rf_length(outBlockSXP));
  d.nSolve = 0;
  d.nFail = 0;
  if ((int)d.outBlock.size() != d.ny) {
    Rf_error("outBlock has %d entries, expected %d",
             (int)d.outBlock.size(), d.ny);
  }
  dims_registry()[Rf_asInteger(handleSXP)] = d;
  return R_NilValue;
}

SEXP rxstanClearDims(SEXP handleSXP) {
  dims_registry().erase(Rf_asInteger(handleSXP));
  return R_NilValue;
}

// Exercises the exact path Stan takes, from R.  Used by the test suite so the
// C ABI is covered without compiling a Stan model.
SEXP rxstanSolve(SEXP handleSXP, SEXP pSXP) {
  const int handle = Rf_asInteger(handleSXP);
  int ny = 0, np = 0, nBlock = 0;
  if (rxs_layout(handle, &ny, &np, &nBlock, NULL) != 0) {
    Rf_error("unknown rxstan handle %d", handle);
  }
  if (Rf_length(pSXP) != np) {
    Rf_error("expected %d parameters, got %d", np, (int)Rf_length(pSXP));
  }

  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, ny, nBlock + 1));
  char err[512];
  const int rc =
      rxs_solve_sens(handle, REAL(pSXP), np, REAL(out), REAL(out) + ny, ny, err,
                     static_cast<int>(sizeof(err)));
  if (rc != 0) {
    UNPROTECT(1);
    Rf_error("rxs_solve_sens failed (%d): %s", rc, err);
  }
  UNPROTECT(1);
  return out;
}

// Same, but forced down the slow R route, so tests can assert the two agree.
SEXP rxstanSolveSlow(SEXP handleSXP, SEXP pSXP) {
  const int handle = Rf_asInteger(handleSXP);
  int ny = 0, np = 0, nBlock = 0;
  if (rxs_layout(handle, &ny, &np, &nBlock, NULL) != 0) {
    Rf_error("unknown rxstan handle %d", handle);
  }
  if (Rf_length(pSXP) != np) {
    Rf_error("expected %d parameters, got %d", np, (int)Rf_length(pSXP));
  }

  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, ny, nBlock + 1));
  char err[512];
  err[0] = '\0';
  const int rc = solve_via_r(handle, REAL(pSXP), np, REAL(out), REAL(out) + ny,
                             ny, nBlock, err, static_cast<int>(sizeof(err)));
  if (rc != 0) {
    UNPROTECT(1);
    Rf_error("rxs_solve_sens (slow) failed (%d): %s", rc, err);
  }
  UNPROTECT(1);
  return out;
}

}  // extern "C"
