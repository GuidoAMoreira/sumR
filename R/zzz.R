# Selecting precompiled
#' @importFrom matrixStats logSumExp
precomp_select <- function(name, pars){
  if (name == "negbin_marginal"){
    stopifnot(length(pars) == 4)
    lL <- log(pars[1]) -
      matrixStats::logSumExp(c(log(pars[1]), log(pars[2]))) + log1p(- pars[3])
    nmbr = 1L
  } else if (name == "noObs"){
    stopifnot(length(pars) == 1)
    lL <- log1p(-pars[1])
    nmbr = 2L
  } else if (name == "COMP"){
    stopifnot(length(pars) == 2)
    lL <- -Inf
    nmbr = 3L
  } else if (name == "dR0"){
    stopifnot(length(pars) == 4)
    lL <- log(pars[1]) + log1p(-pars[4]) +
      (1 + pars[2]) * (log1p(pars[2]) - log(pars[1] + pars[2]))
    nmbr = 4L
  } else if (name == "PL_diff"){
    stopifnot(length(pars) == 3)
    lL <- log(.9999)
    nmbr = 5L
  } else if (name == "negbin_sentinel"){
    stopifnot(length(pars) == 3)
    lL <- log(pars[1]) -
      matrixStats::logSumExp(c(log(pars[1]), log(pars[2]))) +
      log1p(-pars[3])
    nmbr = 6L
  } else if (name == "poisson_sentinel"){
    stopifnot(length(pars) == 2)
    logL <- -Inf
    nmbr = 7L
  } else if (name == "weird_series_constL"){
    stopifnot(length(pars) == 1)
    lL <- log(pars[1])
    nmbr = 8L
  } else if (name == "weird_series"){
    stopifnot(length(pars) == 1)
    lL <- -1
    nmbr = 9L
  }

  list(selection = nmbr, lL = lL)
}

# Conway-Maxwell-Poisson normalizing constant
COMP <- function(k, Theta) k * log(Theta[1]) - Theta[2] * lfactorial(k)
