#include "sumR_internal.h"
#include "math.h"

long double infiniteErrorBoundingPairs_(long double logFun(long k, double *Theta),
                       double *params, double logL, double eps,
                       long maxIter, long n0, long* n)
{
  // Declaration
  long nMax;
  long double maxA, lEps = logl(eps) + LOG_2, total = 0.,
    totalBack = 0., log1mL = Rf_logspace_sub(0, logL), c = 0., cb = 0.,
    *logFunVal = R_Calloc((size_t)(maxIter + 1), long double);
  *n = 0;

  // Find the maximum
  // Finding function max. Only check convergence after max is reached
  logFunVal[*n] = logFun(n0, params);
  while (!R_FINITE(logFunVal[*n])) // In case the series starts with inf values.
    logFunVal[++*n] = logFun(++n0, params);

  do
    logFunVal[++*n] = logFun(++n0, params);
  while (logFunVal[*n] >= logFunVal[*n - 1] && *n <= (maxIter - 1));

  // If too many iterations. Last iter is max.
  if (*n == maxIter)
  {
    partial_logSumExp(logFunVal, maxIter - 1, logFunVal[*n], &c, 0, &total);
    long double result = logFunVal[*n] + log1p(total);
    
    if (logFunVal != NULL) R_Free(logFunVal);
    return result;
  }

  // I know which is the max due to the stop criteria.
  // Assumed local max = global max.
  maxA = logFunVal[*n - 1];
  nMax = *n;
  if (*n > 1)
    partial_logSumExp(logFunVal, *n - 2, maxA, &c, 0, &total);

  do
    logFunVal[++*n] = logFun(++n0, params);
  while  ( (log1mL ? // if L = 0 there's a simpler convergence check
              delta(logFunVal[*n] + logFunVal[*n - 1] -
              Rf_logspace_sub(logFunVal[*n - 1], logFunVal[*n]),
              logFunVal[*n], log1mL) >= lEps :
              (logFunVal[*n - 1] - logFunVal[*n] < Rf_log1pexp(logFunVal[*n] -
                lEps))) &
                (*n < maxIter));
  // Braden bounds
  KahanSum(&totalBack, expl(logFunVal[*n] - log1mL - LOG_2 - maxA), &cb);
  KahanSum(&totalBack, expl(logFunVal[*n] + logFunVal[*n - 1] -
    Rf_logspace_sub(logFunVal[*n - 1], logFunVal[*n]) - LOG_2 -
    maxA), &cb);
  partial_logSumExp(&logFunVal[nMax], *n - nMax - 1, maxA, &cb, 1, &totalBack);
  
  if (logFunVal != NULL) R_Free(logFunVal);
  return maxA + log1pl(total + totalBack);
}
