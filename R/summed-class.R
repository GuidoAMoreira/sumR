#' S3 Class for objects containing iterated summations
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
#' \item{\code{maxReached}}{\code{TRUE} or \code{FALSE}. Indicates whether the
#' maximum iterations was reached.}
#' }
#' @seealso \code{\link{infiniteSum}}, \code{\link{infiniteSum_batches}} and
#' \code{\link{finiteSum}} for available methods.
NULL

#' @rdname summed-objects
#' @param x The \code{summed} object.
#' @param ... Currently unused.
#' @return For \code{print}: The invisible object.
#' @method print summed
#' @export
print.summed <- function(x, ...) {
  cat("Method ", x$method, " performed ",
      ifelse(x$maxReached, "the maximum of ", ""),
      x$n, " iterations and reached ",
      "the sum in the log scale: ", x$sum, "\n", sep = "")
  invisible(x)
}

#' @rdname summed-objects
#' @param x The \code{summed} object.
#' @param ... Currently unused.
#' @return For \code{as.numeric}/\code{as.double}: The approximated sum.
#' @method as.double summed
#' @export
as.double.summed <- function(x, ...) x$sum
