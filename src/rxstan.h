#ifndef RXSTAN_SRC_RXSTAN_H
#define RXSTAN_SRC_RXSTAN_H

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif

// --- C ABI consumed by the Stan translation unit via R_GetCCallable ---------
//
// All of these are strictly no-longjmp: any R-level error raised while solving
// is caught and reported through the return code, so a C++ caller (Stan's
// log_prob, mid-reverse-sweep) is never longjmp'd past its destructors.
//
// Parameters are organised in BLOCKS.  A population model gives each subject
// its own block of `nBlock` parameters, and because subjects are independent
// the Jacobian is block diagonal: output i depends only on block
// `outBlock[i]`.  `dydp` therefore carries just that block, ny x nBlock, not
// the mostly-zero ny x np dense matrix.  Parameters shared across subjects are
// the same thing with a single block.

// Writes the registered layout for `handle`; any pointer may be NULL.
// np == nBlock * nBlocks.  Returns 0 on success, non-zero if unknown.
int rxs_layout(int handle, int *ny, int *np, int *nBlock, int *nBlocks);

// Shorthand for rxs_layout(handle, ny, np, NULL, NULL).
int rxs_dims(int handle, int *ny, int *np);

// Fills `out` (length ny) with the parameter block each output belongs to.
int rxs_out_block(int handle, int *out, int n);

// Solves `handle` at parameter vector `p` (length np, block-major).
// On success writes `ny` values into `y` and the ny x nBlock block-diagonal
// Jacobian into `dydp`, COLUMN-MAJOR (dydp[i + ny*j] = d y_i / d p_block[j]),
// and returns 0.  On failure returns non-zero and writes a NUL-terminated
// message to `errbuf`.
//
// The outputs are ODE STATES at the observation times, never derived `lhs`
// quantities: rxode2 only emits rx__sens_<state>_BY_<param>__, so anything
// derived has to be composed in Stan where its chain rule is autodiffed.
int rxs_solve_sens(int handle, const double *p, int np, double *y, double *dydp,
                   int ny, char *errbuf, int errlen);

// Solve counters for `handle`.  Kept even when diagnostics are silenced, so a
// quiet run still reports how often the solver failed.
int rxs_stats(int handle, double *nSolve, double *nFail);

// --- .Call entry points ----------------------------------------------------
SEXP rxstanProbeRxode2(void);
SEXP rxstanSetDims(SEXP handleSXP, SEXP nySXP, SEXP nBlockSXP, SEXP nBlocksSXP,
                   SEXP outBlockSXP);
SEXP rxstanClearDims(SEXP handleSXP);
SEXP rxstanSolve(SEXP handleSXP, SEXP pSXP);
SEXP rxstanSolveSlow(SEXP handleSXP, SEXP pSXP);
SEXP rxstanStats(SEXP handleSXP);
SEXP rxstanResetStats(SEXP handleSXP);
SEXP rxstanSetSilent(SEXP silentSXP);

// --- Path A: re-drive a solve rxode2 already built (src/fast.cpp) ----------
int rxs_fast_available(int handle);
void rxs_fast_invalidate(void);
int rxs_fast_solve(int handle, const double *p, int np, double *y, double *dydp,
                   int ny, char *errbuf, int errlen);

SEXP rxstanFastSetup(SEXP handleSXP, SEXP sensIdxSXP, SEXP outIdxSXP,
                     SEXP sensStateSXP, SEXP nobsSXP, SEXP blockOfSXP,
                     SEXP quietSXP);
SEXP rxstanFastInvalidate(void);
SEXP rxstanFastAvailable(SEXP handleSXP);
SEXP rxstanProbeSolve(void);
SEXP rxstanProbeParams(void);
SEXP rxstanProbeStates(void);

#ifdef __cplusplus
}
#endif

#endif
