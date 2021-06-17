#include <Rinternals.h>
#include "r_wrappers.h"

SEXP inf_sum(SEXP logFun, SEXP params, SEXP eps, SEXP maxIter,
             SEXP logL, SEXP n0, SEXP rho, SEXP forceAlgo, SEXP forceMax)
{
  defineVar(install("Theta"), params, rho);
  long double out;
  R_xlen_t n;

  // Global variables declared in r_wrappers.h. They are used in translator
  envir = rho;
  lF = logFun;

  out = algorithm_selector(translator, REAL(params), REAL(eps)[0],
                           INTEGER(maxIter)[0], REAL(logL)[0], INTEGER(n0)[0],
                           INTEGER(forceAlgo)[0], &n, INTEGER(forceMax)[0]);

  return retFun((double)out, n);
}

SEXP infinite_sum_callPrecomp(SEXP lF, SEXP params, SEXP epsilon, SEXP maxIter,
                              SEXP logL, SEXP n0, SEXP forceAlgo, SEXP forceMax)
{
  long double out;
  R_xlen_t n;

  out = algorithm_selector(precompiled_selector(INTEGER(lF)[0]),
                           REAL(params), REAL(epsilon)[0],
                           INTEGER(maxIter)[0], REAL(logL)[0],
                           INTEGER(n0)[0], INTEGER(forceAlgo)[0], &n,
                           INTEGER(forceMax)[0]);

  return retFun((double)out, n);


}

//////// c-folding wrappers

SEXP inf_c_folding(SEXP logFun, SEXP params, SEXP eps, SEXP maxIter,
                   SEXP n0, SEXP rho, SEXP c, SEXP N_start, SEXP forceMax)
{
  defineVar(install("Theta"), params, rho);
  long double out;
  R_xlen_t n;

  // Global variables declared in r_wrappers.h. They are used in translator
  envir = rho;
  lF = logFun;

  out = infiniteCFolding_(translator, REAL(params), REAL(eps)[0],
                          INTEGER(maxIter)[0], INTEGER(n0)[0], &n,
                          INTEGER(c)[0], INTEGER(N_start)[0],
                          INTEGER(forceMax)[0]);

  return retFun((double)out, n);
}

SEXP infinite_c_folding_precomp(SEXP lF, SEXP params, SEXP epsilon,
                                SEXP maxIter, SEXP n0, SEXP c, SEXP N_start,
                                SEXP forceMax)
{
  long double out;
  R_xlen_t n;


  out = infiniteCFolding_(precompiled_selector(INTEGER(lF)[0]),
                          REAL(params), REAL(epsilon)[0],
                          INTEGER(maxIter)[0], INTEGER(n0)[0], &n,
                          INTEGER(c)[0], INTEGER(N_start)[0],
                          INTEGER(forceMax)[0]);

  return retFun((double)out, n);
}
