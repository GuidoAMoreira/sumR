# sumR <img src="man/figures/logo.jpeg" align="right" height="139"/>

Approximate infinite series with a guaranteed error margin!

# Installation

Installing this package in R from this GitHub page is straightforward. Make sure you have the package `devtool` installed and that your Operating System has a C compiler that R can access. This is automatically true in a standard installation in Linux, requires Xcode or gcc (depends on MacOS version) in Mac, and Rtools in Windows. Then, just run:

``` r
devtools::install_github("GuidoAMoreira/sumR")
```

# Interfacing low-level summation with Rcpp

It is possible to define a series at the C level and use a `Rcpp` wrapper function. For this to be possible, all that is needed is the use of the *depends* Rcpp attribute and including the adequate header. The following code at the R level exports a function that sums a series.

``` r
library(Rcpp)

sourceCpp(code='
#include <Rcpp.h>

// [[Rcpp::depends(sumR)]]

#include <sumR.h>

long double some_series(long n, double *p)
{
  long double out = n * log1pl(-p[0]);
  return out;
}

// [[Rcpp::export]]
double sum_series(double param)
{
  double parameter[1];
  long double r;
  long n; // Number of iterations. Does not require initialization.

  parameter[0] = param;

  r = infiniteSum(some_series, parameter, log1p(-parameter[0]), 0, exp(-35), 100000, 0, &n);
  
  Rcpp::Rcout << "Summation took " << n << " iterations to converge.\\n";

  return (double)r;
}
')

sum_series(0.08)
```

The available functions and their arguments are listed below.

# Interfacing low-level summation with another package

Package `sumR` facilitates using its low-level C or C++ function in other packages. In order to use it, these steps are necessary:

1.  Make sure that the DESCRIPTION file in your package includes `sumR` in its **LinkingTo** and **Imports** fields.
2.  Make sure that the NAMESPACE file in your packages includes a line with `import(sumR)`. If you are using the roxygen2 documentation package, this can be achieved by adding `#' @import sumR` in one of your R files, such as mypackage-package.R, before running `roxygen2::roxygenize()`.
3.  Include sumR's API header file, sumRAPI.h, in your C or C++ file that will use the desired sumR low-level function.
4.  Use sumR's functions in your code.

The following code exemplifies a C file in a package after steps 1. and 2. above were taken. See [Writing R Extensions](https://cran.r-project.org/doc/manuals/r-release/R-exts.html) to learn about the `SEXP` type and related macros and functions.

``` c
#include <Rinternals.h>
#include <Rmath.h> // Required for the log1p and log1pl functions
#include <sumR.h>

long double some_series(long n, double *p)
{
  long double out = n * log1pl(-p[0]);
  return out;
}

SEXP sum_series(SEXP param)
{
  double parameter[1];
  long double r;
  long n; // Number of iterations. Doesn't require initialization.

  parameter[0] = REAL(param)[0];

  r = infiniteSum(some_series, parameter, log1p(-parameter[0]), 0, exp(-35), 100000, 0, &n);
  
  Rprintf("Summation took %d iterations to converge.\n", n);
  
  return Rf_ScalarReal((double)r);
}
```

Then your package can have an R wrapper function such as:

``` r
#' Wrapper function that sums a series for a given parameter
#' @export
sumSeries <- function(p) .Call("sum_series", param = p, PACKAGE = "mypackage")
```

The this function can be tested after the package has been installed with:

``` r
sumSeries(0.08)
```

# Available sumR functions

The interfaced functions from sumR are:

``` c
long double infiniteSum(long double logFun(long k, double *Theta), double *params, double logL, int alternating, double eps, long maxIter, long n0, long* n);
long double infiniteSumToThreshold(long double logFun(long k, double *Theta), double *params, int alternating, double eps, long maxIter, long n0, long* n);
long double infiniteErrorBoundingPairs(long double logFun(long k, double *Theta), double *params, double logL, double eps, long maxIter, long n0, long* n);
long double infiniteBatches(long double logFun(long k, double *Theta), double *params, long batch_size, double eps, long maxIter, long n0, long* n);
long double sumNTimes(long double logFun(long k, double *Theta), double *params, long n, long n0);
```

Function `infiniteSum` dispatches the arguments to `infiniteSumToThreshold`, `infiniteErrorBoundingPairs` or `infiniteBatches` depending on the value of `logL` and returns the result of the respectively chosen function. Namely, it is the first if `logL` < log(0.5), the second if `logL` < 0 and the third otherwise.

See the help documentation in the sumR package for information about the interfaced function arguments. `sumNTimes` is documented under `finiteSum`.

## Notes

When making a wrapper function, we have found that manually typecasting the result of the low-level C function to double before passing it to R is more stable in some systems than straight up using Rf_ScalarReal on the long double result.

Since all `sumR` code is in C, a C++ functor cannot be passed as the logFun argument, not even in a templated wrapper.

None of the C functions do error checking. It is the user's responsibility to pass adequate arguments. Pay particular mind to the C stack limit (the functions create an array of `maxIter` size).

If defining an R function with `cppFunction`, be mindful that argument `k` of `logFun` is a single `long`. If multiple values are expected, they should be dealt with in the wrapper.
