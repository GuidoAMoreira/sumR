#include "sumR.h"
#include "math.h"

long double infiniteSum_(
    long double logFun(long k, double *Theta),
    double *params, int alternating, double eps,
    long maxIter, double logL, long n0, long* n){
  if (logL < - LOG_2 || alternating)
    return infiniteSumToThreshold_(logFun, params, alternating, eps, maxIter, n0, n);
  else if (logL < 0)
    return infiniteAdaptive_(logFun, params, eps, maxIter, logL, n0, n);
  else
    return infiniteCFolding_(logFun, params, eps, maxIter, n0, n, 2, 20);
}
