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
#' Cmax(dta, 'Conc')
#' Cmax(x, 'Concentration')
#'
#' @export
Cmax <- function(dta, Conc = 'Conc') {

  # Rename columns
  names(dta)[names(dta) == Conc] <- "Conc"

  # Cmax
  return(max(dta$Conc))
}
