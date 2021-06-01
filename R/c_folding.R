#' Approximates the sum of a positive discrete infinite series with a single
#' maximum using the c-folding algorithm
#'
#' A simple method to perform the summation. It adds the values in batches and
#' stops when the accumulated batch is smaller than the desired threshold. There
#' is an implementation purely in \code{R} and one in \code{C}. The one in R is
#' usually slightly faster due to vectorized computing.
#' @param logFunction The function that returns the series value
#' \ifelse{html}{\out{a<sub>n</sub>}}{\eqn{a_n}} in
#' the log scale. Can either be an \code{R} function or a string indicating one
#' of the pre-coded functions. See \code{\link{precompiled}} for a list of
#' available functions. If defined in \code{R}, the function's definition must
#' have two arguments. The first argument must be the integer argument
#' equivalent to \eqn{n} in \ifelse{html}{\out{a<sub>n</sub>}}{\eqn{a_n}} and
#' the second must be a vector of numeric parameters.
#' @param parameters A numeric vector with parameters used in logFunction.
#' Vectorized summation over various parameter values sets is not implemented.
#' Use \code{\link{apply}} or their variants to achieve this.
#' @param epsilon The desired error margin for the approximation. See 'details'.
#' @param maxIter The maximum number of iterations for the approximation. In
#' most cases, this number will not be reached unless it is very small.
#' @param n0 The sum will be approximated for the series starting at this value.
#' @param c The fold by which \code{N_start} is multiplied to calculate the
#' checkpoints. See 'details'.
#' @param N_start The size of the first batch. It is multiplied by \code{c} to
#' find the checkpoints. See 'details'.
#' @return A list with two named members, \code{sum} and \code{n}. \code{sum} is
#' the approximated value in the log scale and \code{n} is the total number of
#' iterations, that is, the number of times the function was evaluated.
#' @seealso \code{\link{precompiled}} provides a list with precompiled functions
#' that can be used for the summation. \code{\link{infiniteSum}} is a more
#' efficient algorithm.
#' @details The series \ifelse{html}{\out{a<sub>n</sub>}}{\eqn{a_n}} must pass
#' the ratio convergence test, meaning that the ratio
#' \ifelse{html}{\out{a<sub>n+1</sub>/a<sub>n</sub>}}{\eqn{a_{n+1}/a_n}} must
#' converge to a number \eqn{L < 1} when \eqn{n} goes to infinity. The
#' approximation can be proven to overshoot the necessary number of function
#' evaluations when \eqn{L < 0.5}, in the sense that its result would be within
#' \code{epsilon} distance of the true value with much fewer evalutaions. If
#' \eqn{0.5 \le L < 1}, then there is no theoretical proof of overshooting, but
#' practical experimentation has shown that this is the case in many examples.
#'
#' The c-folding algorithm consists of evaluating the function a fixed number of
#' times for two checkpoints. If the difference between the sum at these
#' checkpoints is smaller than \code{epsilon}, the code stops and the later
#' checkpoint sum is returned. Else, continue summing until the next checkpoint.
#' The first checkpoint is \code{N_start}. The second is \code{N_start} *
#' \code{c}. Subsequent ones are \code{N_start} * 2\code{c}, \code{N_start} *
#' 3\code{c} and so on.
#'
#' This function's efficiency is reliant on the choice of \code{N_start} and
#' \code{c}. If they are set too large, the algorithm overshoots the necessary
#' number of function evaluations too much. If they are set too small, the
#' algorithm will need to process too many partial summations which slows it
#' down. However, if they are well calibrated for the series, they can
#' potentially be very efficient.
#'
#' Since the batch sizes are known before the calculations are made,
#' function evaluations can be vectorized. This is why there are two functions
#' available. \code{infiniteSum_cFolding} does the calculations at the \code{R}
#' level, while \code{infiniteSum_cFolding_C} interfaces the low level \code{C}
#' code. However, the \code{C} code does not use vectorization to reduce
#' dependency on third party libraries, and therefore the \code{R} level
#' function should be faster.
#'
#' Another requirement in the current installment of this function is that the
#' series must have only a single maximum. This is the case for most discrete
#' probability distributions and marginalization problems. This limitation
#' will be addressed in the future.
#' @examples
#' ## Define some function that is known to pass the ratio test.
#' param = 0.1
#' funfun <- function(k, p) return(k * log1p(-p[1]))
#' result <- infiniteSum_cFolding(funfun, parameters = param)
#'
#' ## This series is easy to verify analytically
#' TrueSum = -log(param)
#' TrueSum - result$sum
#' # Notice that it required 400 function evaluations for the approximation.
#' result$n
#'
#' ## A common problem is finding the normalizing constant for the
#' ## Conway-Maxwell-Poisson distribution. It has already been included
#' ## in the precompiled list of functions.
#' comp_params = c(lambda = 5, nu = 3)
#' result <- infiniteSum_cFolding("COMP", comp_params)
#' @importFrom matrixStats logSumExp
#' @export
infiniteSum_cFolding <- function(logFunction, parameters = numeric(),
                                 epsilon = 1e-15, maxIter = 1e5, n0 = 0, c = 2,
                                 N_start = 20){

  stopifnot(is.function(logFunction) || is.character(logFunction),
            length(logFunction) == 1,
            is.numeric(parameters),
            is.numeric(epsilon),
            epsilon > 0,
            length(epsilon) == 1,
            is.numeric(maxIter),
            maxIter > 0,
            length(maxIter) == 1,
            is.numeric(n0),
            n0 >= 0,
            length(n0) == 1,
            is.numeric(c),
            c > 1,
            length(c) == 1,
            is.numeric(N_start),
            N_start > 0,
            length(N_start) == 1)

  if (is.character(logFunction)){
    if (logFunction == "COMP"){
      stopifnot(length(parameters) == 2)
      logFunction <- COMP
    }
  }

  # Setup
  lEps = log(epsilon)
  N_inc = c * N_start
  nextCheckPoint <- (n0 + N_start - 1)
  funValues <- logFunction(n0:nextCheckPoint, parameters)
  lastCheckPoint <- nextCheckPoint + 1
  nextCheckPoint <- n0 + N_inc - 1
  increment <- logFunction(lastCheckPoint:nextCheckPoint, parameters)
  n = n0 + N_inc

  # Convergence checking
  while (matrixStats::logSumExp(increment) > lEps && n < maxIter){
    funValues <- c(funValues, increment)
    lastCheckPoint <- nextCheckPoint + 1
    nextCheckPoint <- nextCheckPoint + N_inc - 1
    increment <- logFunction(lastCheckPoint:nextCheckPoint, parameters)
    n <- n + N_inc
  }
  funValues <- c(funValues, increment)

  list(sum = matrixStats::logSumExp(funValues),
       n = n)
}

#' @name infiniteSum_cFolding
#' @export
infiniteSum_cFolding_C <- function(logFunction, parameters = numeric(),
                                   epsilon = 1e-15, maxIter = 1e5, n0 = 0,
                                   c = 2, N_start = 20){

  stopifnot(is.function(logFunction) || is.character(logFunction),
            length(logFunction) == 1,
            is.numeric(parameters),
            is.numeric(epsilon),
            epsilon > 0,
            length(epsilon) == 1,
            is.numeric(maxIter),
            maxIter > 0,
            length(maxIter) == 1,
            is.numeric(n0),
            n0 >= 0,
            length(n0) == 1,
            is.numeric(c),
            c > 1,
            length(c) == 1,
            is.numeric(N_start),
            N_start > 0,
            length(N_start) == 1)

  maxIter <- as.integer(maxIter); n0 <- as.integer(n0)
  c <- as.integer(c); N_start <- as.integer(N_start)

  if (is.character(logFunction)){
    if (logFunction == "negbin_marginal"){
      stopifnot(length(parameters) == 4)
      out <- .Call("infinite_c_folding_precomp",
                   1L, parameters, epsilon, maxIter, n0, c, N_start,
                   PACKAGE = "sumR")
    }
    else if (logFunction == "noObs"){
      stopifnot(length(parameters) == 1)
      out <- .Call("infinite_c_folding_precomp",
                   2L, parameters, epsilon, maxIter, n0, c, N_start,
                   PACKAGE = "sumR")
    }
    else if (logFunction == "COMP"){
      stopifnot(length(parameters) == 2)
      out <- .Call("infinite_c_folding_precomp",
                   3L, parameters, epsilon, maxIter, n0, c, N_start,
                   PACKAGE = "sumR")
    }
    else if (logFunction == "dR0"){
      stopifnot(length(parameters) == 4)
      out <- .Call("infinite_c_folding_precomp",
                   4L, parameters, epsilon, maxIter, n0, c, N_start,
                   PACKAGE = "sumR")
    }
    else if (logFunction == "PL_diff"){
      stopifnot(length(parameters) == 3)
      out <- .Call("infinite_c_folding_precomp",
                   5L, parameters, epsilon, maxIter, n0, c, N_start,
                   PACKAGE = "sumR")
    }
    else if (logFunction == "negbin_sentinel"){
      stopifnot(length(parameters) == 3)
      out <- .Call("infinite_c_folding_precomp",
                   6L, parameters, epsilon, maxIter, n0, c, N_start,
                   PACKAGE = "sumR")
    }
    else if (logFunction == "poisson_sentinel"){
      stopifnot(length(parameters) == 2)
      out <- .Call("infinite_c_folding_precomp",
                   7L, parameters, epsilon, maxIter, n0, c, N_start,
                   PACKAGE = "sumR")
    }
    else if (logFunction == "weird_series_constL"){
      stopifnot(length(parameters) == 1)
      out <- .Call("infinite_c_folding_precomp",
                   8L, parameters, epsilon, maxIter, n0, c, N_start,
                   PACKAGE = "sumR")
    }
    else if (logFunction == "weird_series"){
      stopifnot(length(parameters) == 1)
      out <- .Call("infinite_c_folding_precomp",
                   9L, parameters, epsilon, maxIter, n0, c, N_start,
                   PACKAGE = "sumR")
    }
  } else if(is.function(logFunction)) {
    stopifnot(logL < 0)

    f <- function(k, Theta) logFunction(k, Theta)

    out <- .Call("inf_c_folding",
                 body(f), parameters, epsilon,
                 maxIter, n0, new.env(), c, N_start,
                 PACKAGE = "sumR")
  } else {
    warning('Argument lFun must either be the name of a precompiled function
            or a function. See help("precompiled") to see which functions are
            available.')
    return(list(sum = 0, n = 0))
  }

  out
}
