#include "sumR.h"
#include "math.h"

long double sumNTimes_(long double logFun(long, double *Theta),
                      double *params, long N, long n0)
{
  // Declaration
  long double maxA, logFunVal[N + 1], total = 0., totalBack = 0., c = 0., cb = 0.;
  long n = 0, nMax;

  // Finding function max.
  logFunVal[n] = logFun(n0, params);
  maxA = logFunVal[n]; nMax = n;
  while (!R_FINITE(logFunVal[n]))
    logFunVal[++n] = logFun(++n0, params);

  do
    logFunVal[++n] = logFun(++n0, params);
  while (logFunVal[n] >= logFunVal[n - 1] && n < N);

  if (n == N)
  {
    partial_logSumExp(logFunVal, N - 1, logFunVal[n], &c, 0, &total);
    return logFunVal[n] + log1pl(total);
  }

  // I know which is the max due to the stop criteria.
  // Assumed local max = global max.
  maxA = logFunVal[n - 1];
  nMax = n;
  if (n > 1)
    partial_logSumExp(logFunVal, n - 2, maxA, &c, 0, &total);

  // Calculate the tail. Only loop once.
  do
    logFunVal[++n] = logFun(++n0, params);
  while (n < N);
  partial_logSumExp(&logFunVal[nMax], n - nMax, maxA, &cb, 1, &totalBack);

  return maxA + log1pl(total + totalBack);
}
