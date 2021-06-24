#' sumR: Approximate series withing a desired error margin
#'
#' \if{html}{\figure{logo.jpeg}{options: width="120px"}}{\figure{logo.jpeg}{options: width=0.5in}}\cr
#' The \pkg{sumR} package provides some \code{R} functions that allow the
#' user to approximate the sum of series. It also interfaces the low level
#' \code{C} functions that perform the actual summation, so that other packages
#' can apply them to their own functions written in \code{C} or \code{C++}. The
#' \href{https://github.com/GuidoAMoreira/sumR}{GitHub} page provides a short
#' tutorial on how to achieve this.
#'
#' The underlying code is under frequent revision for improvements in speed,
#' memory use and precision. As the need to implement other summation algorithms
#' arise, they will be added.
#'
#' The theoretical foundations for the approximations will be submitted to a
#' peer-reviewed journal shortly.
#'
#' @author
#' Authors:
#' \itemize{
#' \item{Guido A. Moreira}
#' \item{Luiz Max Carvalho}
#' }
#'
#' @docType package
#' @name sumR
#'
#' @useDynLib sumR, .registration = TRUE
NULL

.onUnload <- function (libpath) {
  library.dynam.unload("sumR", libpath)
}
