#ifndef __SUMR_RWRAPPERS_H__
#define __SUMR_RWRAPPERS_H__

#include <Rinternals.h>
#include "sumR.h"

// Function evaluation
static inline double feval(SEXP lF, SEXP rho)
{return REAL(eval(lF, rho))[0];}

// Function return macro
static inline SEXP retFun(double res, R_xlen_t mI)
{
  SEXP out = PROTECT(allocVector(VECSXP,2));
  SET_VECTOR_ELT(out, 0, Rf_ScalarReal(res));
  SET_VECTOR_ELT(out, 1, Rf_ScalarInteger(mI));

  SEXP nms = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(nms, 0, mkChar("sum"));
  SET_STRING_ELT(nms, 1, mkChar("n"));

  /* assign names to list */
  setAttrib(out, R_NamesSymbol, nms);

  UNPROTECT(2);
  return out;
}

// Needed global var
SEXP envir, lF;

// Wrapping sums for functions defined at the R level
static inline long double translator(R_xlen_t k, double *Theta)
{
  defineVar(install("k"), Rf_ScalarInteger(k), envir);
  return (long double)feval(lF, envir);
}

static inline long double algorithm_selector(
    long double logF(R_xlen_t k, double *T), double* params, double eps,
    R_xlen_t mI, double lL, R_xlen_t n0, int selector, R_xlen_t* n)
{
  if (selector == 2)
    return infiniteAdaptive_(logF, params, eps, mI, lL, n0, n);
  else if (selector == 1)
    return infiniteSumToThreshold_(logF, params, eps, mI, n0, n);
  else
    return infiniteSum_(logF, params, eps, mI, lL, n0, n);
}

SEXP inf_sum(SEXP logFun, SEXP params, SEXP eps, SEXP maxIter,
             SEXP logL, SEXP n0, SEXP rho, SEXP forceAlgo);

// Wrapper for C pre-compiled code
SEXP infinite_sum_callPrecomp(SEXP lF, SEXP params, SEXP epsilon, SEXP maxIter,
                           SEXP logL, SEXP n0, SEXP forceAlgo);

// Wrappers for the c-folding algorithm
SEXP inf_c_folding(SEXP logFun, SEXP params, SEXP eps, SEXP maxIter,
                   SEXP n0, SEXP rho, SEXP c, SEXP N_start);

SEXP infinite_c_folding_precomp(SEXP lF, SEXP params, SEXP epsilon, SEXP maxIter,
                        SEXP n0, SEXP c, SEXP N_start);

#endif
