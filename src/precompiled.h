#ifndef __SUMR_PRECOMPILED_H__
#define __SUMR_PRECOMPILED_H__
#include <Rinternals.h>

long double negbin_marginal(R_xlen_t k, double *Theta);

long double noObs(R_xlen_t k, double *Theta);

long double COMP(R_xlen_t k, double *Theta);

long double dR0(R_xlen_t k, double *Theta);

long double powerLawDiff(R_xlen_t k, double *Theta);

long double negbin_sentinel(R_xlen_t k, double *Theta);

long double poisson_sentinel(R_xlen_t k, double *Theta);

long double weird_series_constL(R_xlen_t k, double *Theta);

long double weird_series(R_xlen_t k, double *Theta);

#endif
