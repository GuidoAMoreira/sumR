#ifndef __SUMR_H__
#define __SUMR_H__

#include <R.h>
#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif

// Approximates infinite sum with an adaptive truncation
long double infiniteAdaptive_(long double logFun(long k, double *Theta), double *params,
                        double eps, long maxIter, double logL, long n0, long *n);

// Approximates infinite sum by summing until added value is smaller than threshold
long double infiniteSumToThreshold_(long double logFun(long k, double *Theta),
                        double*, double, long, long, long*);

// Approximates infinite sum by summing in batches until the batch addes up to
// less than the desired margin
long double infiniteCFolding_(long double logFun(long k, double *Theta),
                              double*, double, long, long, long*,
                              long, long);

// Dispatches infiniteAdaptive or inviniteSumToThreshold based on logL value
long double infiniteSum_(
    long double logFun(long k, double *Theta),
    double *params, double eps,
    long maxIter, double logL, long n0, long* n);

// Sum N times with no convergence checking
long double sumNTimes_(
    long double logFun(long k, double *Theta),
    double *params, long N, long n0);

#ifdef __cplusplus
}
#endif

#endif
