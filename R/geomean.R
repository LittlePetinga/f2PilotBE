#' Geometric Mean
#'
#' @description
#' `geomean()` calculates geometric mean, i.e., the \emph{N}th root of the product of the \emph{N}
#' observations, equivalent to `exp(mean(log(x)))`.
#'
#' @param x An `R` object.
#' @param na.rm A logical value indicating whether `NA` values should be stripped
#'   before the computation proceeds. Defaults to `FALSE`.
#'
#' @return Value of the geometric mean.
#'
#'
#' @examples
#' x <- c(49.66, 53.62, 52.25, 50.44, 46.84, 50.51, 47.26, 48.90, 49.83, 49.67)
#' geomean(x, TRUE)
#'
#' @export
geomean <- function(x, na.rm = FALSE) {
  return(exp(mean(log(x), na.rm = na.rm)))
}
