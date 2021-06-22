#ifndef __SUMR_H__
#define __SUMR_H__

#include <R.h>
#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif

// Approximates infinite sum with an adaptive truncation
long double infiniteAdaptive_(long double (R_xlen_t k, double *Theta), double*,
                        double, R_xlen_t, double, R_xlen_t, R_xlen_t*);

// Approximates infinite sum by summing until added value is smaller than threshold
long double infiniteSumToThreshold_(long double (R_xlen_t k, double *Theta),
                        double*, double, R_xlen_t, R_xlen_t, R_xlen_t*);

// Approximates infinite sum by summing in batches until the batch addes up to
// less than the desired margin
long double infiniteCFolding_(long double (R_xlen_t k, double *Theta),
                              double*, double, R_xlen_t, R_xlen_t, R_xlen_t*,
                              R_xlen_t, R_xlen_t);

// Dispatches infiniteAdaptive or inviniteSumToThreshold based on logL value
long double infiniteSum_(
    long double logFun(R_xlen_t k, double *Theta),
    double *params, double eps,
    R_xlen_t maxIter, double logL, R_xlen_t n0, R_xlen_t* n);

// Sum N times with no convergence checking
long double sumNTimes_(
    long double logFun(R_xlen_t k, double *Theta),
    double *params, R_xlen_t N, R_xlen_t n0);

#ifdef __cplusplus
}
#endif

#endif
