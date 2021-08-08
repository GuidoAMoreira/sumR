#include "sumR.h"
#include "math.h"

long double infiniteSum_(
    long double logFun(long k, double *Theta),
    double *params, int alternating, double eps,
    long maxIter, double logL, long n0, long* n)
{return logL < - LOG_2 || alternating ?
  infiniteSumToThreshold_(logFun, params, alternating, eps, maxIter, n0, n) :
  infiniteAdaptive_(logFun, params, eps, maxIter, logL, n0, n);}
