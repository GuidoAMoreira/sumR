#include <Rinternals.h>
#include "precompiled.h"
#include "r_wrappers.h"

SEXP inf_sum(SEXP logFun, SEXP params, SEXP eps, SEXP maxIter,
                           SEXP logL, SEXP n0, SEXP rho, SEXP forceAlgo)
{
  defineVar(install("Theta"), params, rho);
  long double out;
  R_xlen_t n;

  // Global variables declared in r_wrappers.h. They are used in translator
  envir = rho;
  lF = logFun;

  out = algorithm_selector(translator, REAL(params), REAL(eps)[0],
                           INTEGER(maxIter)[0], REAL(logL)[0], INTEGER(n0)[0],
                           INTEGER(forceAlgo)[0], &n);

  return retFun((double)out, n);
}

SEXP infinite_sum_callPrecomp(SEXP lF, SEXP params, SEXP epsilon, SEXP maxIter,
                              SEXP logL, SEXP n0, SEXP forceAlgo)
{
  unsigned int funSelect = INTEGER(lF)[0];
  long double out;
  R_xlen_t n;

  switch (funSelect)
  {
    case 1:
      out = algorithm_selector(negbin_marginal, REAL(params), REAL(epsilon)[0],
                               INTEGER(maxIter)[0], REAL(logL)[0],
                               INTEGER(n0)[0], INTEGER(forceAlgo)[0], &n);
      break;
    case 2:
      out = algorithm_selector(noObs, REAL(params), REAL(epsilon)[0],
                               INTEGER(maxIter)[0], REAL(logL)[0],
                               INTEGER(n0)[0], INTEGER(forceAlgo)[0], &n);
      break;
    case 3:
      out = algorithm_selector(COMP, REAL(params), REAL(epsilon)[0],
                                INTEGER(maxIter)[0], REAL(logL)[0],
                                INTEGER(n0)[0], INTEGER(forceAlgo)[0], &n);
      break;
    case 4:
      out = algorithm_selector(dR0, REAL(params), REAL(epsilon)[0],
                                INTEGER(maxIter)[0], REAL(logL)[0],
                                INTEGER(n0)[0], INTEGER(forceAlgo)[0], &n);
      break;
    case 5:
      out = algorithm_selector(powerLawDiff, REAL(params), REAL(epsilon)[0],
                                INTEGER(maxIter)[0], REAL(logL)[0],
                                INTEGER(n0)[0], INTEGER(forceAlgo)[0], &n);
      break;
    case 6:
      out = algorithm_selector(negbin_sentinel, REAL(params), REAL(epsilon)[0],
                               INTEGER(maxIter)[0], REAL(logL)[0],
                               INTEGER(n0)[0], INTEGER(forceAlgo)[0], &n);
      break;
    case 7:
      out = algorithm_selector(poisson_sentinel, REAL(params), REAL(epsilon)[0],
                               INTEGER(maxIter)[0], REAL(logL)[0],
                               INTEGER(n0)[0], INTEGER(forceAlgo)[0], &n);
      break;
    case 8:
      out = algorithm_selector(weird_series_constL, REAL(params), REAL(epsilon)[0],
                               INTEGER(maxIter)[0], REAL(logL)[0],
                               INTEGER(n0)[0], INTEGER(forceAlgo)[0], &n);
      break;
    case 9:
      out = algorithm_selector(weird_series, REAL(params), REAL(epsilon)[0],
                               INTEGER(maxIter)[0], REAL(logL)[0],
                               INTEGER(n0)[0], INTEGER(forceAlgo)[0], &n);
      break;
    default:
      error("No implemented logFunction found.");
  }

  return retFun((double)out, n);


}

//////// c-folding wrappers

SEXP inf_c_folding(SEXP logFun, SEXP params, SEXP eps, SEXP maxIter,
                   SEXP n0, SEXP rho, SEXP c, SEXP N_start)
{
  defineVar(install("Theta"), params, rho);
  long double out;
  R_xlen_t n;

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
  unsigned int funSelect = INTEGER(lF)[0];
  long double out;
  R_xlen_t n;

  switch (funSelect)
  {
    case 1:
      out = infiniteCFolding_(negbin_marginal, REAL(params), REAL(epsilon)[0],
                               INTEGER(maxIter)[0], INTEGER(n0)[0], &n,
                               INTEGER(c)[0], INTEGER(N_start)[0]);
      break;
    case 2:
      out = infiniteCFolding_(noObs, REAL(params), REAL(epsilon)[0],
                               INTEGER(maxIter)[0], INTEGER(n0)[0], &n,
                               INTEGER(c)[0], INTEGER(N_start)[0]);
      break;
    case 3:
      out = infiniteCFolding_(COMP, REAL(params), REAL(epsilon)[0],
                               INTEGER(maxIter)[0], INTEGER(n0)[0], &n,
                               INTEGER(c)[0], INTEGER(N_start)[0]);
      break;
    case 4:
      out = infiniteCFolding_(dR0, REAL(params), REAL(epsilon)[0],
                              INTEGER(maxIter)[0], INTEGER(n0)[0], &n,
                              INTEGER(c)[0], INTEGER(N_start)[0]);
      break;
    case 5:
      out = infiniteCFolding_(powerLawDiff, REAL(params), REAL(epsilon)[0],
                              INTEGER(maxIter)[0], INTEGER(n0)[0], &n,
                              INTEGER(c)[0], INTEGER(N_start)[0]);
      break;
    case 6:
      out = infiniteCFolding_(negbin_sentinel, REAL(params), REAL(epsilon)[0],
                              INTEGER(maxIter)[0], INTEGER(n0)[0], &n,
                              INTEGER(c)[0], INTEGER(N_start)[0]);
      break;
    case 7:
      out = infiniteCFolding_(poisson_sentinel, REAL(params), REAL(epsilon)[0],
                              INTEGER(maxIter)[0], INTEGER(n0)[0], &n,
                              INTEGER(c)[0], INTEGER(N_start)[0]);
      break;
    case 8:
      out = infiniteCFolding_(weird_series_constL, REAL(params), REAL(epsilon)[0],
                              INTEGER(maxIter)[0], INTEGER(n0)[0], &n,
                              INTEGER(c)[0], INTEGER(N_start)[0]);
      break;
    case 9:
      out = infiniteCFolding_(weird_series, REAL(params), REAL(epsilon)[0],
                              INTEGER(maxIter)[0], INTEGER(n0)[0], &n,
                              INTEGER(c)[0], INTEGER(N_start)[0]);
      break;
    default:
      error("No implemented logFunction found.");
  }

  return retFun((double)out, n);
}
