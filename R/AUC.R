#' Area Under the Curve (AUC)
#'
#' @description `AUC()` calculates the cumulative area under the concentration-time curve (AUC) over time.
#'
#' @param dta Dataframe with concentration-time data.
#' @param Time Name of the column with time data.
#' @param Conc Name of the column with concentration data.
#' @param method Method for the AUC calculation:
#'   * `Linear-Up/Log-Down` (default method): Linear Trapezoidal Rule for for Increasing Values, Log Trapezoidal Rule for Decreasing Values.
#'   * `Linear/Log`: Linear Trapezoidal Rule for for Values below Cmax, Log Trapezoidal Rule for Values above Cmax.
#'   * `Linear`: Linear Trapezoidal Rule.
#'
#' @return Dataframe with cumulative AUC over time.
#'
#'
#' @examples
#' dta <- data.frame(
#'   Time = c(0, 0.25, 0.5, 0.75, 1, 1.5, 1.75, 2, 2.25, 2.5,
#'            2.75, 3, 3.25, 3.5, 3.75, 4, 6, 8, 12, 24),
#'   Conc = c(0.00, 221.23, 377.19, 494.73, 555.74, 623.86, 615.45, 663.38, 660.29, 621.71,
#'            650.33, 622.28, 626.72, 574.94, 610.51, 554.02, 409.14, 299.76, 162.85, 27.01))
#' AUC(dta, Time = 'Time', Conc = 'Conc')
#'
#' @export
AUC <- function(dta, Time = 'Time', Conc = 'Conc', method = 'Linear-Up/Log-Down') {

  # Rename columns
  names(dta)[names(dta) == Time] <- "Time"
  names(dta)[names(dta) == Conc] <- "Conc"

  # Order data by time
  dta <- dta[order(dta$Time),]


  # Methods for AUC calculation

  Linear.AUC <- function(C1, C2, t1, t2) {
    # Linear Trapezoidal Method for AUC Calculation
    # The linear trapezoidal method uses linear interpolation, between a given time interval (t1 - t2),
    # to calculate the AUC
    return(((C2 + C1)*(t2 - t1))/2)
  }

  Log.AUC <- function(C1, C2, t1, t2) {
    # Logarithmic Trapezoidal Method for AUC Calculation
    # The logarithmic trapezoidal method uses logarithmic interpolation, between a given time interval
    # (t1 - t2), to calculate the AUC. This method is more accurate when concentrations are decreasing,
    # since drug follows an exponential elimination process, i.e. linear on a logarithmic scale
    return(((C1 - C2)/(log(C1, base = exp(1))-log(C2, base = exp(1))))*(t2 - t1))
  }

  # For the first timepoint, AUC corresponds to the concentration
  dta.AUC <- data.frame(Time = min(dta$Time),
                        AUC  = dta$Conc[which.min(dta$Time)])

  # For the subsequent timepoints, a trapezoidal method is applied
  for (t in 2:length(dta$Time)) {
    if (method == 'Linear') {
      dta.AUC.i <- data.frame(Time = dta$Time[t],
                              AUC  = Linear.AUC(C1 = dta$Conc[which(dta$Time == dta$Time[t-1])],
                                                C2 = dta$Conc[which(dta$Time == dta$Time[t])],
                                                t1 = dta$Time[t-1],
                                                t2 = dta$Time[t]))

    } else if (method == 'Linear-Up/Log-Down') {
      if (dta$Conc[which(dta$Time == dta$Time[t])] >= dta$Conc[which(dta$Time == dta$Time[t-1])]) {
        # Linear for increasing concentration
        dta.AUC.i <- data.frame(Time = dta$Time[t],
                                AUC  = Linear.AUC(C1 = dta$Conc[which(dta$Time == dta$Time[t-1])],
                                                  C2 = dta$Conc[which(dta$Time == dta$Time[t])],
                                                  t1 = dta$Time[t-1],
                                                  t2 = dta$Time[t]))
      } else {
        # Log for decreasing concentratios
        dta.AUC.i <- data.frame(Time = dta$Time[t],
                                AUC  = Log.AUC(C1 = dta$Conc[which(dta$Time == dta$Time[t-1])],
                                               C2 = dta$Conc[which(dta$Time == dta$Time[t])],
                                               t1 = dta$Time[t-1],
                                               t2 = dta$Time[t]))

        if (dta.AUC.i$AUC == Inf) {
          dta.AUC.i$AUC  <- Linear.AUC(C1 = dta$Conc[which(dta$Time == dta$Time[t-1])],
                                       C2 = dta$Conc[which(dta$Time == dta$Time[t])],
                                       t1 = dta$Time[t-1],
                                       t2 = dta$Time[t])
        }
      }

    } else if (method == 'Linear/Log') {
      if (dta$Conc[which(dta$Time == dta$Time[t])] <= max(dta$Conc)) {
        # Linear for concentrations below Cmax
        dta.AUC.i <- data.frame(Time = dta$Time[t],
                                AUC  = Linear.AUC(C1 = dta$Conc[which(dta$Time == dta$Time[t-1])],
                                                  C2 = dta$Conc[which(dta$Time == dta$Time[t])],
                                                  t1 = dta$Time[t-1],
                                                  t2 = dta$Time[t]))
      } else {
        # Log for concentrations above Cmax
        dta.AUC.i <- data.frame(Time = dta$Time[t],
                                AUC  = Log.AUC(C1 = dta$Conc[which(dta$Time == dta$Time[t-1])],
                                               C2 = dta$Conc[which(dta$Time == dta$Time[t])],
                                               t1 = dta$Time[t-1],
                                               t2 = dta$Time[t]))
      }

    }
    # Append all AUC data
    dta.AUC <- rbind(dta.AUC,dta.AUC.i)
  }
  # Cumulative sum of all AUCs
  dta.AUC$AUC <- cumsum(dta.AUC$AUC)

  # Output
  return(dta.AUC)
}
