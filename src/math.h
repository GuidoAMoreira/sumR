#ifndef __SUMR_MATH_H__
#define __SUMR_MATH_H__

#include <Rinternals.h>

#ifndef LOG_2
#define LOG_2 0.69314718055994530941723212145817656831
#endif

// R internal functions that require definition
double Rf_logspace_sub(double, double);
double Rf_log1pexp(double);

// Auxiliary function. NOT A FULL logSumExp IMPLEMENTATION
void partial_logSumExp(long double*, long, long double,
                       long double*, int, long double*);

// This function helps reduce the floating point rounding error.
static inline void KahanSum(long double* tot, long double x, long double* c)
{
  long double t, y;
  y = x - *c;
  t = *tot + y;
  *c = t - *tot - y;
  *tot = t;
}

// Required in the adaptive algorithm
static inline long double delta(double logz, double loga, double log1ml)
{
  double ls = loga - log1ml;
  return (logz > ls ? Rf_logspace_sub(logz, ls) : Rf_logspace_sub(ls, logz));
}

#endif
