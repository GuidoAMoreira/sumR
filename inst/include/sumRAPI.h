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
    long double logFun(long k, double *Theta),
    double *params, double eps,
    long maxIter, double logL, long n0, long* n) {
  static long double (*fun)(long double (long k, double *Theta),
                      double*, double, long, double, long,
                      long*) = NULL;
  fun = (long double(*)(long double (long k, double *Theta),
                     double*, double, long, double, long,
                     long*))
    R_GetCCallable("sumR", "infiniteSum_");
  return fun(logFun, params, eps, maxIter, logL, n0, n);
}

// Sum-To-Threshold algorithm
long double infiniteSumToThreshold(
    long double logFun(long k, double *Theta),
    double *params, double eps,
    long maxIter, long n0, long* n) {
  static long double (*fun)(long double (long k, double *Theta),
                      double*, double, long, long,
                      long*) = NULL;
  fun = (long double(*)(long double (long k, double *Theta),
                     double*, double, long, long,
                     long*))
    R_GetCCallable("sumR", "infiniteSumToThreshold_");
  return fun(logFun, params, eps, maxIter, n0, n);
}

// Adaptive Truncation algorithm
long double infiniteAdaptive(
    long double logFun(long k, double *Theta),
    double *params, double eps,
    long maxIter, double logL, long n0, long* n) {
  static long double (*fun)(long double (long k, double *Theta),
                      double*, double, long, double, long,
                      long*) = NULL;
  fun = (long double(*)(long double (long k, double *Theta),
                     double*, double, long, double, long,
                     long*))
    R_GetCCallable("sumR", "infiniteAdaptive_");
  return fun(logFun, params, eps, maxIter, logL, n0, n);
}

// c-folding algorithm
long double infiniteCFolding(
    long double logFun(long k, double *Theta),
    double *params, double eps,
    long maxIter, long n0, long* n, long c, long N_start) {
  static long double (*fun)(long double (long k, double *Theta),
                      double*, double, long, long,
                      long*, long, long) = NULL;
  fun = (long double(*)(long double (long k, double *Theta),
                     double*, double, long, long,
                     long*, long, long))
    R_GetCCallable("sumR", "infiniteCFolding_");
  return fun(logFun, params, eps, maxIter, n0, n, c, N_start);
}

// sum N times
long double sumNTimes(
    long double logFun(long k, double *Theta),
    double *params, long n, long n0) {
  static long double (*fun)(long double (long k, double *Theta),
                      double*, long, long) = NULL;
  fun = (long double(*)(long double (long k, double *Theta),
                     double*, long, long))
    R_GetCCallable("sumR", "sumNTimes_");
  return fun(logFun, params, n, n0);
}

#ifdef __cplusplus
}
#endif

#endif
