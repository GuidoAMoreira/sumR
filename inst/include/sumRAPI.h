#ifndef __SUMR_API_H__
#define __SUMR_API_H__

#include <sumR.h>		// also includes R.h and Rinternals.h

#include <Rconfig.h>
#include <R_ext/Rdynload.h>

#ifdef __cplusplus
extern "C" {
#endif

// infiniteSum auto select
long double infiniteSum(
    long double logFun(R_xlen_t k, double *Theta),
    double *params, double eps,
    R_xlen_t maxIter, double logL, R_xlen_t n0, R_xlen_t* n) {
  static long double (*fun)(long double (R_xlen_t k, double *Theta),
                      double*, double, R_xlen_t, double, R_xlen_t,
                      R_xlen_t*) = NULL;
  fun = (long double(*)(long double (R_xlen_t k, double *Theta),
                     double*, double, R_xlen_t, double, R_xlen_t,
                     R_xlen_t*))
    R_GetCCallable("sumR", "infiniteSum_");
  return fun(logFun, params, eps, maxIter, logL, n0, n);
}

// Sum-To-Threshold algorithm
long double infiniteSumToThreshold(
    long double logFun(R_xlen_t k, double *Theta),
    double *params, double eps,
    R_xlen_t maxIter, R_xlen_t n0, R_xlen_t* n) {
  static long double (*fun)(long double (R_xlen_t k, double *Theta),
                      double*, double, R_xlen_t, R_xlen_t,
                      R_xlen_t*) = NULL;
  fun = (long double(*)(long double (R_xlen_t k, double *Theta),
                     double*, double, R_xlen_t, R_xlen_t,
                     R_xlen_t*))
    R_GetCCallable("sumR", "infiniteSumToThreshold_");
  return fun(logFun, params, eps, maxIter, n0, n);
}

// Adaptive Truncation algorithm
long double infiniteAdaptive(
    long double logFun(R_xlen_t k, double *Theta),
    double *params, double eps,
    R_xlen_t maxIter, double logL, R_xlen_t n0, R_xlen_t* n) {
  static long double (*fun)(long double (R_xlen_t k, double *Theta),
                      double*, double, R_xlen_t, double, R_xlen_t,
                      R_xlen_t*) = NULL;
  fun = (long double(*)(long double (R_xlen_t k, double *Theta),
                     double*, double, R_xlen_t, double, R_xlen_t,
                     R_xlen_t*))
    R_GetCCallable("sumR", "infiniteAdaptive_");
  return fun(logFun, params, eps, maxIter, logL, n0, n);
}

// c-folding algorithm
long double infiniteCFolding(
    long double logFun(R_xlen_t k, double *Theta),
    double *params, double eps,
    R_xlen_t maxIter, R_xlen_t n0, R_xlen_t* n, R_xlen_t c, R_xlen_t N_start) {
  static long double (*fun)(long double (R_xlen_t k, double *Theta),
                      double*, double, R_xlen_t, R_xlen_t,
                      R_xlen_t*, R_xlen_t, R_xlen_t) = NULL;
  fun = (long double(*)(long double (R_xlen_t k, double *Theta),
                     double*, double, R_xlen_t, R_xlen_t,
                     R_xlen_t*, R_xlen_t, R_xlen_t))
    R_GetCCallable("sumR", "infiniteCFolding_");
  return fun(logFun, params, eps, maxIter, n0, n, c, N_start);
}

// sum N times
long double sumNTimes(
    long double logFun(R_xlen_t k, double *Theta),
    double *params, R_xlen_t n, R_xlen_t n0) {
  static long double (*fun)(long double (R_xlen_t k, double *Theta),
                      double*, R_xlen_t, R_xlen_t) = NULL;
  fun = (long double(*)(long double (R_xlen_t k, double *Theta),
                     double*, R_xlen_t, R_xlen_t))
    R_GetCCallable("sumR", "sumNTimes_");
  return fun(logFun, params, n, n0);
}

#ifdef __cplusplus
}
#endif

#endif
