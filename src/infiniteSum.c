#include "sumR.h"
#include "math.h"

long double infiniteSum_(
    long double logFun(R_xlen_t k, double *Theta),
    double *params, double eps,
    R_xlen_t maxIter, double logL, R_xlen_t n0, R_xlen_t* n, int forceMax)
{return logL < - LOG_2 ?
  infiniteSumToThreshold_(logFun, params, eps, maxIter, n0, n, forceMax) :
  infiniteAdaptive_(logFun, params, eps, maxIter, logL, n0, n, forceMax);}
