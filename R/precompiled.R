#' List of precompiled functions in the sumR package
#'
#' More functions are periodically added to this list for convenience and speed.
#' These functions are all evaluated in the log scale and pass the ratio test,
#' that is, the limit of
#' \ifelse{html}{\out{a<sub>n+1</sub>/a<sub>n</sub>}}{\eqn{a_{n+1}/a_n}} as n
#' goes to infinity is a value \eqn{0 \le L < 1}. The value of \eqn{L} is
#' indicated in each entry. It is calculated automatically when the precompiled
#' functions are used in the summation.
#'
#' @note In some cases, the sum of the series is known in closed form for some
#' values of the parameters. The package function does not check for these cases
#' and just performs the approximation. If the exact value is desired by the
#' user when it is known, they must take responsibility for checking and
#' providing these values.
#'
#' Another important thing to note is that the precompiled functions perform all
#' calculations with twice the numerical precision than R. Therefore, in some
#' cases, there might be very small differences in the sum when comparing the
#' results of the function using the precompiled function and the same function
#' defined at the R level.
#'
#' @section Conway-Maxwell-Poisson normalizing constant:
#' This series is the kernel of the Conway-Maxwell-Poisson distribution, which
#' generalizes the Poisson and Geometric distributions. Its form is
#'
#' \ifelse{html}{\out{<center> a<sub>n</sub> = λ<sup>n</sup> / (n!)<sup>ν</sup> </center>}}{\deqn{a_n = \frac{\lambda^n}{(n!)^\nu},}}
#'
#' \deqn{L = 0, log(L) = -\infty,}
#'
#' for \eqn{\lambda > 0} and \eqn{\nu > 0}.
#'
#' When \eqn{\nu = 1}, this series reduces to the Poisson distribution kernel
#' and the sum (in the log scale) is known to be \eqn{\lambda}. When
#' \eqn{\nu = 0} and \eqn{0 < \lambda < 1}, the series reduces to the Geometric
#' distribution kernel with parameter \eqn{1 - \lambda}. The series is known to
#' sum to \eqn{1}. Finally, as \eqn{\nu} goes to \eqn{\infty} the distribution
#' approaches a Bernoulli distribution with parameter
#' \eqn{\lambda / (1 - \lambda)}.
#'
#' Another known result is when \eqn{\nu = 2}, in which case the sum is the
#' modified Bessel function of the first kind of order 0 evaluated at
#' \eqn{2\sqrt{\lambda}}.
#'
#' - String to access the precompiled function: \code{"COMP"}.
#'
#' - parameter vector: \code{c(lambda, nu)}.
#'
#' @section Double Poisson normalizing constant:
#' This series is the kernel of the double Poisson distribution, which is a
#' special case of the double exponential family, which extends it. Its form is
#'
#' \ifelse{html}{\out{<center> a<sub>n</sub> = φ<sup>0.5</sup> e<sup>-φ µ</sup> (e<sup>-n</sup>n<sup>n</sup> / n!)(eµ / n)<sup>φ n</sup> </center>}}{\deqn{a_n = \sqrt{\phi}e^{-\phi\mu}\left(\frac{e^{-n}n^n}{n!}\right)\left(\frac{e\mu}{n}\right)^{\phi n},}}
#'
#' \deqn{L = 0, log(L) = -\infty,}
#'
#' for \eqn{\lambda > 0} and \eqn{\phi > 0}.
#'
#' When \eqn{\phi = 1}, this series reduces to the Poisson distribution kernel
#' and the sum (in the log scale) is known to be 0.
#'
#' - String to acccess the precompiled function: \code{"double_poisson"}.
#'
#' - parameter vector: \code{c(mu, phi)}
#'
#' @name precompiled
#' @seealso \code{\link{infiniteSum}}, \code{\link{finiteSum}} and
#' \code{\link{infiniteSum_cFolding}}
NULL
