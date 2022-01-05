#' Class for objects containing iterated summations
#' 
#' Contains the summations in the log scale. The value can either be an
#' approximation to an infinite series or a finite sum.
#' @name summed-objects
#' @section Elements in the list:
#' \describe{
#' \item{\code{sum}}{The resulting sum in the log scale.}
#' \item{\code{n}}{The performed number of iterations. This value represents
#' the number of series elements evaluations performed during the summation.}
#' \item{\code{method}}{The method used for the summation.}
#' }
#' @seealso \code{\link{infiniteSum}}, \code{\link{infiniteSum_cFolding}} and
#' \code{\link{finiteSum}} for available methods.
NULL

#' @rdname summed-objects
#' @param x The \code{summed} object.
#' @param ... Currently unused.
#' @return The invisible object.
#' @method print summed
#' @export
print.summed <- function(x, ...) {
  cat("Method ", x$method, " performed ", x$n, " iterations and reached ",
      "the sum in the log scale: ", x$sum, sep = "")
  invisible(x)
}
