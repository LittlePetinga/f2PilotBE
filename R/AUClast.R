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
#' dta <- data.frame(
#'   Time = c(0, 0.25, 0.5, 0.75, 1, 1.5, 1.75, 2, 2.25, 2.5,
#'            2.75, 3, 3.25, 3.5, 3.75, 4, 6, 8, 12, 24),
#'   Conc = c(0.00, 221.23, 377.19, 494.73, 555.74, 623.86, 615.45, 663.38, 660.29, 621.71,
#'            650.33, 622.28, 626.72, 574.94, 610.51, 554.02, 409.14, 299.76, 162.85, 27.01))
#' AUClast(dta, 't', 'Conc')
#'
#' @export
AUClast <- function(dta, Time = 'Time', Conc = 'Conc', method = 'Linear-Up/Log-Down') {
  AUC.cum <- AUC(dta, Time = Time, Conc = Conc, method = method)
  AUClast <- AUC.cum[which.max(AUC.cum$Time),'AUC']
  return(AUClast)
}
