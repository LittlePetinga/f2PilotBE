#' AUC to Last Timepoint
#'
#' @description `AUClast()` calculates area under the concentration-time curve (AUC)
#' from time zero to the time of the last measurable concentration (\eqn{AUC_{\text{0-t}}}).
#'
#' @param dta Dataframe with concentration-time data.
#' @param Time Name of the column with time data.
#' @param Conc Name of the column with concentration data.
#' @param method Method for the AUC calculation:
#'   * `Linear-Up/Log-Down` (default method): Linear Trapezoidal Rule for for Increasing Values, Log Trapezoidal Rule for Decreasing Values.
#'   * `Linear/Log`: Linear Trapezoidal Rule for for Values below \eqn{C_{\text{max}}}, Log Trapezoidal Rule for Values above \eqn{C_{\text{max}}}.
#'   * `Linear`: Linear Trapezoidal Rule.
#'
#' @return Value of \eqn{AUC_{\text{0-t}}}.
#'
#'
#' @examples
#' AUClast(dta, 't', 'Conc')
#' AUClast(x, 'Time', 'Concentration', 'Linear/Log')
#'
#' @export
AUClast <- function(dta, Time = 'Time', Conc = 'Conc', method = 'Linear-Up/Log-Down') {
  AUC.cum <- AUC(dta, Time = Time, Conc = Conc, method = method)
  AUClast <- AUC.cum[which.max(AUC.cum$Time),'AUC']
  return(AUClast)
}
