#' Maximum Concentration (Cmax)
#'
#' @description
#' `Cmax()` calculates the maximum observed concentration post-dose \eqn{C_{\text{max}}}
#' directly obtained from the observed concentration-time profile.
#'
#' @param dta Dataframe with concentration-time data.
#' @param Conc Name of the column with concentration data.
#'
#' @return Value of the maximum observed concentration.
#'
#'
#' @examples
#' dta <- data.frame(
#'   Time = c(0, 0.25, 0.5, 0.75, 1, 1.5, 1.75, 2, 2.25, 2.5,
#'            2.75, 3, 3.25, 3.5, 3.75, 4, 6, 8, 12, 24),
#'   Conc = c(0.00, 221.23, 377.19, 494.73, 555.74, 623.86, 615.45, 663.38, 660.29, 621.71,
#'            650.33, 622.28, 626.72, 574.94, 610.51, 554.02, 409.14, 299.76, 162.85, 27.01))
#' Cmax(dta, 'Conc')
#'
#' @export
Cmax <- function(dta, Conc = 'Conc') {

  # Rename columns
  names(dta)[names(dta) == Conc] <- "Conc"

  # Cmax
  return(max(dta$Conc))
}
