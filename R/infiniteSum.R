#' Approximates the sum of a positive discrete infinite series with a single
#' maximum
#'
#' For series that pass the ratio test, the approximation is analytically
#' guaranteed to have an error that is smaller than epsilon. This can
#' occasionally not happen due to floating point arithmetic.
#' @param logFunction The function that returns the series value
#' \ifelse{html}{\out{a<sub>n</sub>}}{\eqn{a_n}} in
#' the log scale. Can either be an \code{R} function or a string indicating one
#' of the precompiled functions. See \code{\link{precompiled}} for a list of
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
#' @param logL The log of the limit value of
#' \ifelse{html}{\out{a<sub>n+1</sub>/a<sub>n</sub>}}{\eqn{a_{n+1}/a_n}} which
#' must be smaller than 1, or smaller than 0 in the log scale. See 'details'.
#' @param n0 The sum will be approximated for the series starting at this value.
#' @param forceAlgorithm A value to control which summation algorithm to use.
#' See 'details'.
#' @return A list with two named members, \code{sum} and \code{n}. \code{sum} is
#' the approximated value in the log scale and \code{n} is the total number of
#' iterations, that is, the number of times the function was evaluated.
#' @details The approximated sum is based on some theoretical results which,
#' analytically, guarantee that the approximation will be within \code{epsilon}
#' distance to the true value. It is possible that the numerical result fails
#' to fall in this distance due to floating point arithmetic. The \code{C} code
#' under the hood is being continuously reviewed to minimize this problem. They
#' seem to occur more often when the series decays very fast to zero or when the
#' total is a large number.
#'
#' For these theoretical results to work, the series must pass the ratio test,
#' which means that the ratio
#' \ifelse{html}{\out{a<sub>n+1</sub>/a<sub>n</sub>}}{\eqn{a_{n+1}/a_n}} must
#' converge to a number \eqn{L < 1} when \eqn{n} goes to infinity. The log of
#' \eqn{L} should be provided to the function for a better approximation.
#' Otherwise, a warning will be given and the program will try to estimate its
#' value, however, the theoretical result that guarantees the approximation does
#' not necessarily hold in this case.
#'
#' Another requirement in the current installment of this function is that the
#' series must have only a single maximum. This is the case for most discrete
#' probability distributions and marginalization problems. This limitation
#' will be addressed in the future.
#'
#' There are currently two implemented algorithms that perform the calculations.
#' The first, called Sum-To-Threshold, sums the series values until the series
#' values are smaller than \code{epsilon}. This is the fastest algorithm, but
#' it is only guaranteed to provide an approximation within the desired error
#' margin when \eqn{L < 0.5}.
#'
#' The second algorithm, called Adaptive Truncation is based on a more general
#' result which works for any \eqn{0 \le L < 1}. This algorithm sums the series
#' until
#'
#' \ifelse{html}{\out{<center>|a<sub>n+1</sub>/(1-L) - a<sub>n+1</sub> a<sub>n</sub>/(a<sub>n</sub> - a<sub>n+1</sub>)| < 2 ε.</center>}}{\deqn{\left|\frac{a_{n+1}}{1-L} - \frac{a_{n+1}a_n}{a_n-a_{n+1}}\right| < 2 \epsilon.}}
#'
#' Then the approximation is the added values of the sum plus
#'
#' \ifelse{html}{\out{<center> 0.5 (a<sub>n+1</sub>/(1-L) + a<sub>n+1</sub> a<sub>n</sub>/(a<sub>n</sub> - a<sub>n+1</sub>))</center>}}{\deqn{0.5 (\frac{a_{n+1}}{1-L} + \frac{a_{n+1}a_n}{a_n-a_{n+1}}).}}
#'
#' The Adaptive Truncation method usually requires less function evaluations
#' than the Sum-To-Threshold one, however the convergence checking is more
#' demanding, which means that it is typically slower, albeit slightly. If
#' \eqn{L = 0}, the convergence checking can be reduced and the Adaptive
#' Truncation becomes almost as fast as the Sum-To-Threshold method.
#'
#' The \code{forceAlgorithm} parameter can be used to
#' control which algorithm to use. When it is 0, the program automatically
#' selects the Sum-To-Threshold when \eqn{L < 0.5} and the Adaptive Truncation
#' otherwise. If it is set to 1, then the Sum-To-Threshold algorithm is forced.
#' If it is 2, then the Adaptive Truncation is forced. A small note, the
#' Adaptive Truncation algorithm can go up to \code{maxIter} + 1 function
#' evaluations. This is due to its convergence checking dependence on
#' \ifelse{html}{\out{a<sub>n+1</sub>}}{\eqn{a_{n+1}}}.
#'
#' The function to be summed can be an R function or a string naming the
#' precompiled function in the package. The list of precompiled functions can
#' be found in \code{\link{precompiled}}, and more functions will be added in
#' time. As is intuitive, using a precompiled function is much faster than using
#' an \code{R} function. In fact, it has been observed to be dozens time faster.
#'
#' The advanced user can program their own precompiled functions and use the
#' package's summation algorithms by linking the appropriate header file. See
#' the \href{https://github.com/GuidoAMoreira/sumR}{GitHub} readme for the names to use.
#' @seealso \code{\link{precompiled}} provides a list with precompiled functions
#' that can be used for the summation.
#' @examples
#' ## Define some function that is known to pass the ratio test.
#' param = 0.1
#' funfun <- function(k, p) return(k * log1p(-p[1]))
#' result <- infiniteSum(funfun, parameters = param, logL = log1p(-param))
#'
#' ## This series is easy to verify analytically
#' TrueSum = -log(param)
#' TrueSum - result$sum
#' # Since exp(logL) = 0.9, the Adaptive Truncation
#' # algorithm is used. Notice that it only required
#' # 2 function evaluations for the approximation, that is
#' result$n
#'
#' ## A common problem is finding the normalizing constant for the
#' ## Conway-Maxwell-Poisson distribution. It has already been included
#' ## in the precompiled list of functions.
#' comp_params = c(lambda = 5, nu = 3)
#' result <- infiniteSum("COMP", comp_params)
#' @export
infiniteSum <- function(logFunction, parameters = numeric(), epsilon = 1e-15,
                        maxIter = 1e5, logL = NULL, n0 = 0, forceAlgorithm = 0){

  stopifnot(is.function(logFunction) || is.character(logFunction),
            length(logFunction) == 1,
            is.numeric(parameters),
            is.numeric(epsilon),
            epsilon > 0,
            length(epsilon) == 1,
            is.numeric(maxIter),
            maxIter > 0,
            length(maxIter) == 1,
            length(logL) %in% 0:1,
            is.numeric(n0),
            n0 >= 0,
            length(n0) == 1,
            forceAlgorithm %in% 0:2)

  test_logL <- logL
  if (forceAlgorithm == 1){
    if (!is.null(logL))
      warning("Sum-To-Threshold algorithm doesn't use parameter logL. It will be ignored.")
    logL <- -1 # Any negative value to pass on to C.
  }
  maxIter <- as.integer(maxIter); n0 <- as.integer(n0)
  forceAlgorithm <- as.integer(forceAlgorithm)

  if (is.character(logFunction)){
    if (!is.null(test_logL)) warning("Summation over precompiled functions uses pre-determined logL. Inputted value ignored.")
    selection <- precomp_select(logFunction, parameters)
    out <- .Call("infinite_sum_callPrecomp", selection$selection, parameters,
                 epsilon, maxIter, selection$lL, n0, forceAlgorithm,
                 PACKAGE = "sumR")
  } else if(is.function(logFunction)) {
    if (is.null(logL) && forceAlgorithm != 1)
      stop('Parameter logL is NULL. Please provide its value. See help("infiniteSum") for details.')
    stopifnot(logL < 0)

    f <- function(k, Theta) logFunction(k, Theta)

    out <- .Call("inf_sum",
                body(f), parameters, epsilon,
                maxIter, logL,
                n0, new.env(), forceAlgorithm,
                PACKAGE = "sumR")
  } else {
    warning('Argument lFun must either be the name of a precompiled function or a function. See help("precompiled") to see which functions are available.')
    return(list(sum = 0, n = 0))
  }

  out
}
