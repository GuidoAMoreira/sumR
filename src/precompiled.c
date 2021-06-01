#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "precompiled.h"
#include "math.h"

double Rf_dnbinom(double x, double size, double prob, int give_log);
double Rf_dbinom(double x, double n, double p, int give_log);
double Rf_dpois(double x, double lambda, int give_log);

long double negbin_marginal(R_xlen_t k, double *Theta)
{
  double x = Theta[3], s = Theta[1];
  return (k<x ? -INFINITY :
            Rf_dnbinom(k, s, s / (s + Theta[0]), 1) +
              Rf_dbinom(x, k, Theta[2], 1));
}

long double noObs(R_xlen_t k, double *Theta)
{return k * log1p(-Theta[0]);}

long double COMP(R_xlen_t k, double *Theta)
{return k * log(Theta[0]) - Theta[1] * lgamma1p(k);}

long double dR0(R_xlen_t k, double *Theta)
{
  double x = Theta[2];
  if (k == 0 || k < x) {return -INFINITY;}
  else
  {
    long double wy = Theta[1] * k, wyym1 = wy + k - 1;
    return lgamma(wyym1) - (lgamma(wy) + Rf_lgamma1p(k)) + ((k-1) * (log(Theta[0])-
                  log(Theta[1])) - wyym1 * log1p(Theta[0]/Theta[1])) +
                  Rf_dbinom(x, k, Theta[3], 1);
  }
}

long double powerLawDiff(R_xlen_t k, double *Theta)
{return (k<Theta[1] ? -INFINITY :
           -Theta[0] * log(k) + Rf_logspace_sub(0,
                           -Theta[2] - Theta[3] * (k - Theta[1])));}

                           long double negbin_sentinel(R_xlen_t k, double *Theta)
                           {
                             return Rf_dnbinom(k, Theta[1], Theta[1]/ (Theta[1] + Theta[0]), 1) +
                               k * log1p(-Theta[2]);
                           }

long double poisson_sentinel(R_xlen_t k, double *Theta)
{
  return Rf_dpois(k, Theta[0], 1) +
    k * log1p(-Theta[1]);
}

long double weird_series_constL(R_xlen_t k, double *Theta)
{return (k == 0 ? -INFINITY : -(2 * log(k) + k * log(Theta[0])));}

long double weird_series(R_xlen_t k, double *Theta)
{return (k == 0 ? -INFINITY : Rf_lgamma1p(k) - k * log(k));}


