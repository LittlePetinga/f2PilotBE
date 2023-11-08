#' Time of Cmax (Tmax)
#'
#' @description
#' `Tmax()` calculates the time of the maximum observed concentration post-dose (\eqn{t_{\text{max}}}).
#'
#' @param dta Dataframe with concentration-time data.
#' @param Time Name of the column with time data.
#' @param Conc Name of the column with concentration data.
#'
#' @return Value of the time of the maximum observed concentration.
#'
#'
#' @examples
#' Tmax(dta, 't', 'Conc')
#' Tmax(x, 'Time', 'Concentration')
#'
#' @export
Tmax <- function(dta, Time = 'Time', Conc = 'Conc') {

  # Rename columns
  names(dta)[names(dta) == Time] <- "Time"
  names(dta)[names(dta) == Conc] <- "Conc"

  # tmax
  return(dta$Time[which.max(dta$Conc)])
}
