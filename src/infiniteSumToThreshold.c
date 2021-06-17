#include "sumR.h"
#include "math.h"

long double infiniteSumToThreshold_(long double logFun(R_xlen_t k, double *Theta),
                              double *params, double eps,
                              R_xlen_t maxIter, R_xlen_t n0, R_xlen_t* n,
                              int forceMax)
{
  // Declaration
  R_xlen_t nMax;
  long double maxA, lEps = log(eps), logFunVal[maxIter + 1], total = 0.,
    totalBack = 0., c = 0., cb = 0.;
  *n = 0;

  // Finding function max. Only check convergence after max is reached
  logFunVal[*n] = logFun(n0, params);
  while (!R_FINITE(logFunVal[*n]))
    logFunVal[++*n] = logFun(++n0, params);

  do
    logFunVal[++*n] = logFun(++n0, params);
  while (logFunVal[*n] >= logFunVal[*n - 1] && *n <= (maxIter - 1));

  // If too many iterations. Last iter is max.
  if (*n == maxIter)
  {
    partial_logSumExp(logFunVal, maxIter - 1, logFunVal[*n], &c, 0, &total);
    return logFunVal[*n] + log1p(total);
  }

  // I know which is the max due to the stop criteria.
  // Assumed local max = global max.
  maxA = logFunVal[*n - 1];
  nMax = *n;
  if (*n > 1)
    partial_logSumExp(logFunVal, *n - 2, maxA, &c, 0, &total);

  // Calculate the tail. Only loop once.
  do
    logFunVal[++*n] = logFun(++n0, params);
  while ((logFunVal[*n] >= lEps || forceMax) && (*n <= (maxIter - 1)));
  partial_logSumExp(&logFunVal[nMax], *n - nMax, maxA, &cb, 1, &totalBack);

  return maxA + log1pl(total + totalBack);
}
