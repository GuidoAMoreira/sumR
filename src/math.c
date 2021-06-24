#include <Rinternals.h>
#include "math.h"

void partial_logSumExp(long double* fun, long evals, long double maxA,
                       long double* c, int backwards, long double* res)
{
  // Go backwards if the series is decreasing and forward otherwise.
  if (backwards)
    for (R_xlen_t i = evals; i >= 0; i--)
      KahanSum(res, expl(fun[i] - maxA), c);
  else
    for (R_xlen_t i = 0; i <= evals; i++)
      KahanSum(res, expl(fun[i] - maxA), c);
}
