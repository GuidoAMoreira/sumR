#ifndef __SUMR_RWRAPPERS_H__
#define __SUMR_RWRAPPERS_H__

#include <Rinternals.h>
#include "sumR.h"
#include "precompiled.h"

// Function evaluation
static inline double feval(SEXP lF, SEXP rho)
{return REAL(eval(lF, rho))[0];}

// Function return macro
static inline SEXP retFun(double res, long mI)
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
static inline long double translator(long k, double *Theta)
{
  defineVar(install("k"), Rf_ScalarInteger(k), envir);
  return (long double)feval(lF, envir);
}

// Selectors

static inline long double algorithm_selector(
    long double logF(long k, double *T), double* params, double eps,
    long mI, double lL, long n0, int selector, long* n)
{
  if (selector == 2)
    return infiniteAdaptive_(logF, params, eps, mI, lL, n0, n);
  else if (selector == 1)
    return infiniteSumToThreshold_(logF, params, eps, mI, n0, n);
  else
    return infiniteSum_(logF, params, eps, mI, lL, n0, n);
}

typedef long double (*lFptr)(long, double*);
static inline lFptr precompiled_selector(unsigned int funS){
  switch (funS){
  case 1:
    return negbin_marginal;
    break;
  case 2:
    return noObs;
    break;
  case 3:
    return COMP;
    break;
  case 4:
    return dR0;
    break;
  case 5:
    return powerLawDiff;
    break;
  case 6:
    return negbin_sentinel;
    break;
  case 7:
    return poisson_sentinel;
    break;
  case 8:
    return weird_series_constL;
    break;
  case 9:
    return weird_series;
    break;
  default:
    error("Compiled function not found.");
  }
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

// Wrappers for the fixed iterations algorithm
SEXP sum_n_times_precomp(SEXP lF, SEXP params, SEXP N, SEXP n0);

#endif
