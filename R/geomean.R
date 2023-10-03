#' Geometric Mean
#'
#' @description
#' `geomean()` calculates geometric mean, i.e., the \emph{N}th root of the product of the \emph{N}
#' observations, equivalent to exp(mean(log(x))).
#'
#' @param x An `R` object.
#' @param na.rm A logical value indicating whether `NA` values should be stripped
#'   before the computation proceeds. Defaults to `FALSE`.
#'
#' @return Value of the geometric mean.
#'
#' @author Sara Carolina Henriques
#'
#' @examples
#' geomean(x, TRUE)
#' geomean(x, FALSE)
#'
#' # Can be used for the calculation of geometric mean of concentration data,
#' # for Test and Reference treatments, over time.
#' # Calculate the geometric mean for each Treatment, by timepoint
#' dta <- data.frame(t(tapply(dta_id$Conc,
#'                            dta_id[,c('Treatment','Time')],
#'                            geomean)))
#' dta$Time <- as.numeric(row.names(dta))
#' row.names(dta) <- c()
#'
#' # To stack treatment information
#' dta <- cbind(dta[3],stack(dta[1:2]))
#' names(dta)[names(dta) %in% c('values','ind')] <- c('Conc',
#'                                                    'Treatment')
#'
#' @export
geomean <- function(x, na.rm = FALSE) {
  return(exp(mean(log(x), na.rm = na.rm)))
}
