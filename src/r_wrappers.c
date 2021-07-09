#include <Rinternals.h>
#include "r_wrappers.h"

lFptr precompiled_selector(SEXP funS, double *logL,
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
    if (params[0] <= 0 || params[1] <= 0)
      error("Parameters must be positive.");
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
    if (params[0] <= 0 || params[1] <= 0)
      error("Parameters must be positive.");
    *logL = -INFINITY;
    return dbl_poisson;
  }
  error("Compiled function not found.");
}

SEXP inf_sum(SEXP logFun, SEXP params, SEXP eps, SEXP maxIter,
             SEXP logL, SEXP n0, SEXP rho, SEXP forceAlgo)
{
  defineVar(install("Theta"), params, rho);
  long double out;
  long n;

  // Global variables declared in r_wrappers.h. They are used in translator
  envir = rho;
  lF = logFun;

  out = algorithm_selector(translator, REAL(params), REAL(eps)[0],
                           INTEGER(maxIter)[0], REAL(logL)[0], INTEGER(n0)[0],
                           INTEGER(forceAlgo)[0], &n);

  return retFun((double)out, n);
}

SEXP infinite_sum_callPrecomp(SEXP lF, SEXP params, SEXP epsilon, SEXP maxIter,
                              SEXP n0, SEXP forceAlgo)
{
  long double out;
  long n;
  double lL;
  lFptr logFunction = precompiled_selector(lF, &lL, REAL(params),
                                           Rf_xlength(params));

  out = algorithm_selector(logFunction, REAL(params), REAL(epsilon)[0],
                           INTEGER(maxIter)[0], lL,
                           INTEGER(n0)[0], INTEGER(forceAlgo)[0], &n);

  return retFun((double)out, n);
}

//////// c-folding wrappers

SEXP inf_c_folding(SEXP logFun, SEXP params, SEXP eps, SEXP maxIter,
                   SEXP n0, SEXP rho, SEXP c, SEXP N_start)
{
  defineVar(install("Theta"), params, rho);
  long double out;
  long n;

  // Global variables declared in r_wrappers.h. They are used in translator
  envir = rho;
  lF = logFun;

  out = infiniteCFolding_(translator, REAL(params), REAL(eps)[0],
                          INTEGER(maxIter)[0], INTEGER(n0)[0], &n,
                          INTEGER(c)[0], INTEGER(N_start)[0]);

  return retFun((double)out, n);
}

SEXP infinite_c_folding_precomp(SEXP lF, SEXP params, SEXP epsilon,
                                SEXP maxIter, SEXP n0, SEXP c, SEXP N_start)
{
  long double out;
  long n;
  double logL;

  out = infiniteCFolding_(precompiled_selector(lF, &logL, REAL(params),
                                               Rf_xlength(params)),
                          REAL(params), REAL(epsilon)[0],
                          INTEGER(maxIter)[0], INTEGER(n0)[0], &n,
                          INTEGER(c)[0], INTEGER(N_start)[0]);

  return retFun((double)out, n);
}

//////// fixed iterations wrappers

SEXP sum_n_times_precomp(SEXP lF, SEXP params, SEXP N, SEXP n0)
{
  long double out;
  long n = INTEGER(N)[0];
  double logL;

  out = sumNTimes_(precompiled_selector(lF, &logL, REAL(params),
                                        Rf_xlength(params)),
                   REAL(params), n, INTEGER(n0)[0]);

  return retFun((double)out, n);
}
