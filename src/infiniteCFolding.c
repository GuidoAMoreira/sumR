#include "sumR.h"
#include "math.h"

long double infiniteCFolding_(long double logFun(R_xlen_t k, double *Theta),
                                 double *params, double eps,
                                 R_xlen_t maxIter, R_xlen_t n0, R_xlen_t* n,
                                 R_xlen_t c, R_xlen_t N_start, int forceMax)
{
  // Declaration
  R_xlen_t N, N_inc = N_start * c;
  long double maxA, lEps = logl(eps), logFunVal[maxIter + 1],
          partial = 0., *checkStart = logFunVal, S = 0., cc = 0., total = 0;
  *n = 0;

  logFunVal[*n] = logFun(n0, params);
  while (!R_FINITE(logFunVal[*n]) && *n < (maxIter - 1))
    logFunVal[++*n] = logFun(++n0, params);

  // Find the maximum
  do
    logFunVal[++*n] = logFun(++n0, params);
  while (logFunVal[*n] >= logFunVal[*n - 1] && *n <= (maxIter - 1));

  // If too many iterations. Last iter is max.
  if (*n == maxIter)
  {
    partial_logSumExp(logFunVal, maxIter - 1, logFunVal[*n], &cc, 0, &total);
    return logFunVal[*n] + log1pl(total);
  }

  // I know which is the max due to the stop criteria.
  // Assumed local max = global max.
  maxA = logFunVal[*n - 1];
  lEps -= maxA; // For the convergence checking
  while (*n % N_inc && *n < maxIter) // Complete the next checkpoint.
    logFunVal[++*n] = logFun(++n0, params);
  N = *n == N_inc ? N_start : ((*n - N_inc) / N_inc) * N_inc; // Second to last completed checkpoint.
  partial_logSumExp(logFunVal, N, maxA, &cc, 0, &partial);
  checkStart += N + 1; // pointer displacement
  N = *n == N_inc ? N_inc - N_start - 1 : N_inc - 1; // How many iterations to the last checkpoints.
  cc = 0.;
  partial_logSumExp(checkStart, N, maxA, &cc, 0, &S);

  // Calculate the tail. Only loop once.
  while ((logl(S) >= lEps || forceMax) && *n < maxIter)
  {
    partial += S;
    for (N = 0, S = 0.; N < N_inc; N++)
      KahanSum(&S, exp(logFun(++n0, params) - maxA), &cc);
    *n += N_inc;
  }

  return maxA + logl(partial + S);
}
