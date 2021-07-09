#ifndef __SUMR_RWRAPPERS_H__
#define __SUMR_RWRAPPERS_H__

#include <Rinternals.h>
#include "sumR.h"
#include "precompiled.h"

double Rf_logspace_add(double, double);

// Quality of life macro
#define compareStr(s) (!strcmp(CHAR(STRING_PTR(funS)[0]), s))
#define arraysize(a) sizeof(a) / sizeof(double)
#define checkSize(s, n) if (s != n) error("Wrong number of parameters.\n")

typedef long double (*lFptr)(long, double*);

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

static inline long double algorithm_selector(lFptr logF, double *params,
                                             double eps, long mI, double lL,
                                             long n0, int selector, long *n)
{
  if (selector == 2)
    return infiniteAdaptive_(logF, params, eps, mI, lL, n0, n);
  else if (selector == 1)
    return infiniteSumToThreshold_(logF, params, eps, mI, n0, n);
  else
    return infiniteSum_(logF, params, eps, mI, lL, n0, n);
}

static inline lFptr precompiled_selector(SEXP funS, double *logL,
                                         double *params, R_xlen_t size){

  if compareStr("negbin_marginal"){
    checkSize(size, 4);
    *logL = log(params[0]) - Rf_logspace_add(log(params[0]), log(params[1])) +
      log1p(-params[2]);
    return negbin_marginal;
  }
  if compareStr("noObs"){
    checkSize(size, 1);
    *logL = log1p(-params[0]);
    return noObs;
  }
  if compareStr("COMP"){
    checkSize(size, 2);
    *logL = -INFINITY;
    return COMP;
  }
  if compareStr("dR0"){
    checkSize(size, 4);
    *logL = log(params[0]) + log1p(-params[3]) +
      (1 + params[1]) * (log1p(params[1]) - log(params[0] + params[1]));
    return dR0;
  }
  if compareStr("powerLawDiff"){
    checkSize(size, 3);
    *logL = log(0.9999);
    return powerLawDiff;
  }
  if compareStr("negbin_sentinel"){
    checkSize(size, 3);
    *logL = log(params[0]) - Rf_logspace_add(log(params[0]), log(params[1])) +
      log1p(-params[2]);
    return negbin_sentinel;
  }
  if compareStr("poisson_sentinel"){
    checkSize(size, 2);
    *logL = -INFINITY;
    return poisson_sentinel;
  }
  if compareStr("weird_series_constL"){
    checkSize(size, 1);
    *logL = log(params[0]);
    return weird_series_constL;
  }
  if compareStr("weird_series"){
    checkSize(size, 1);
    *logL = -1;
    return weird_series;
  }
  if compareStr("double_poisson"){
    checkSize(size, 2);
    *logL = -INFINITY;
    return dbl_poisson;
  }
  error("Compiled function not found.");
}

SEXP inf_sum(SEXP logFun, SEXP params, SEXP eps, SEXP maxIter,
             SEXP logL, SEXP n0, SEXP rho, SEXP forceAlgo);

// Wrapper for C pre-compiled code
SEXP infinite_sum_callPrecomp(SEXP lF, SEXP params, SEXP epsilon, SEXP maxIter,
                              SEXP n0, SEXP forceAlgo);

// Wrappers for the c-folding algorithm
SEXP inf_c_folding(SEXP logFun, SEXP params, SEXP eps, SEXP maxIter,
                   SEXP n0, SEXP rho, SEXP c, SEXP N_start);

SEXP infinite_c_folding_precomp(SEXP lF, SEXP params, SEXP epsilon, SEXP maxIter,
                                SEXP n0, SEXP c, SEXP N_start);

// Wrappers for the fixed iterations algorithm
SEXP sum_n_times_precomp(SEXP lF, SEXP params, SEXP N, SEXP n0);

#endif
