#' AUC Similarity f2 Factor
#'
#' @description `f2.AUC()` calculates the \eqn{f_2} Similarity factor for AUC, from concentration data.
#'
#' @param dta Dataframe with concentration data from Test and Reference products.
#' @param Time Name of the column with time data.
#' @param Conc Name of the column with concentration data
#'   (only for stacked data, `Trt.cols = FALSE`).
#' @param Trt Name of the column with treatment/formulation information
#'   (only for stacked data, `Trt.cols = FALSE`).
#' @param Ref Nomenclature in `dta` for the Reference product, defaults to `R`.
#' @param Test Nomenclature in `dta` for the Test product, defaults to `T`.
#' @param method Method for the AUC calculation:
#'   * `Linear-Up/Log-Down` (default method): Linear Trapezoidal Rule for for Increasing Values, Log Trapezoidal Rule for Decreasing Values.
#'   * `Linear/Log`: Linear Trapezoidal Rule for for Values below \eqn{C_{\text{max}}}, Log Trapezoidal Rule for Values above \eqn{C_{\text{max}}}.
#'   * `Linear`: Linear Trapezoidal Rule.
#' @param Trt.cols A logical value indicating whether treatment is presented in
#'   columns/pivoted (`TRUE`), or in rows/stacked (`FALSE`).
#' @param details A logical value indicating whether detailed results will be
#'   presented (`TRUE`), or if the function only returns the \eqn{f_2} factor
#'   (`FALSE`). Defaults to `FALSE`.
#' @param plot A logical value indicating whether graphical representation will be
#'   returned. Defaults to `TRUE`.
#'
#' @return Returns a list with the following elements:
#'   * `Raw Concentration Data`: Dataframe of input concentration data, with
#'     treatment information in columns (pivoted).
#'   * `Cumulative AUC`: Dataframe with cumulative AUC over time, for Test and Reference product, calculated from [AUC()] function.
#'   * `Reference AUClast`: Vector of Reference \eqn{AUC_{\text{0-t}}}.
#'   * `Normalized AUC`: dataframe with normalize Test and Reference AUC
#'     concentrations over time, by Reference \eqn{AUC_{\text{0-t}}}, until \eqn{t_{\text{last}}}
#'   * `AUC f2 Factor`: dataframe from [f2()] function.
#'
#' @author Sara Carolina Henriques
#'
#' @references
#' Henriques, S.C.; Albuquerque, J.; Paixão, P.; Almeida, L.; Silva, N.E. (2023).
#' Alternative Analysis Approaches for the Assessment of Pilot Bioavailability/Bioequivalence Studies.
#' *Pharmaceutics*. *15*(5), 1430.
#' [10.3390/pharmaceutics15051430](https://doi.org/10.3390/pharmaceutics15051430).
#'
#' Henriques, S.C.; Paixão, P.; Almeida, L.; Silva, N.E. (2023).
#' Predictive Potential of C~max~ Bioequivalence in Pilot Bioavailability/Bioequivalence Studies,
#' through the Alternative ƒ~2~ Similarity Factor Method. *Pharmaceutics*. *15*(10), 2498.
#' [10.3390/pharmaceutics15102498](https://doi.org/10.3390/pharmaceutics15102498).
#'
#'
#' @examples
#' # Calculate f2 Factor for AUC
#' # when treatment data is pivoted:
#' f2.AUC(dta, Time = 'Time', Reference = 'R', Test = 'T')
#'
#' # when treatment data is stacked:
#' f2.AUC(dta, Time = 'Time', Conc = 'Conc',
#'        Treatment = 'Treatment', Reference = 'R', Test = 'T')
#'
#' @export
f2.AUC <- function(dta, Time = 'Time', Conc = 'Conc',
                   Trt = 'Treatment', Ref = 'Reference', Test = 'Test',
                   method = 'Linear-Up/Log-Down',
                   Trt.cols = TRUE, details = FALSE, plot = TRUE) {

  # If data is pivoted
  if (Trt.cols) {

    # Rename columns
    cols <- c(Time, Ref, Test)
    dta <- dta[,cols]
    names(dta)[names(dta) == Time] <- "Time"
    names(dta)[names(dta) == Ref]  <- "Reference"
    names(dta)[names(dta) == Test] <- "Test"

  } else { # If data is stacked

    # Rename columns
    cols <- c(Trt, Time, Conc)
    dta <- dta[,cols]
    names(dta)[names(dta) == Trt] <- "Treatment"
    names(dta)[names(dta) == Time] <- "Time"
    names(dta)[names(dta) == Conc] <- "Conc"

    # Data from Test
    dta_T <- dta[dta$Treatment == Test, c('Time','Conc')]
    names(dta_T)[names(dta_T) == 'Conc'] <- 'Test'

    # Data from Reference
    dta_R <- dta[dta$Treatment == Ref, c('Time','Conc')]
    names(dta_R)[names(dta_R) == 'Conc'] <- 'Reference'

    dta <- merge(dta_R,dta_T)
  }


  # Calculate cumulative AUC for Reference
  AUC.R <- AUC(dta, Time = 'Time', Conc = 'Reference', method = method)
  names(AUC.R)[names(AUC.R) == 'AUC'] <- 'Reference'

  # Calculate cumulative AUC for Test
  AUC.T <- AUC(dta, Time = 'Time', Conc = 'Test', method = method)
  names(AUC.T)[names(AUC.T) == 'AUC'] <- 'Test'

  # Merge Reference and Test AUC data
  AUC.dta <- merge(AUC.R,AUC.T)


  # Calculate Reference AUClast
  AUClast.R <- AUClast(dta, Time = 'Time', Conc = 'Reference', method = method)


  # Normalize Test and Reference cumulative AUC over time by Reference AUClast
  # Results presented as percentage (%)
  Normalized <- AUC.dta
  Normalized$Reference <- 100*(Normalized$Reference/AUClast.R)
  Normalized$Test <- 100*(Normalized$Test/AUClast.R)
  Normalized$Difference <- Normalized$Reference - Normalized$Test
  Normalized$Sq.Diff <- (Normalized$Reference - Normalized$Test)^2



  # Calculate f2 Factor for AUC
  f2.AUC <- f2(Normalized)


  # Output
  out <- list(
    'Raw Concentration Data' = dta,
    'Cumulative AUC' = AUC.dta,
    'Reference AUClast' = AUClast.R,
    'Normalized AUC' = Normalized,
    'AUC f2 Factor' = f2.AUC
  )


  if (plot) {

    # Plot of normalized concentration over time, until Reference Tmax
    Norm.plot <- (ggplot(data=Normalized)
                  + geom_line(aes(x = Time,
                                  y = Test,
                                  colour = "Test"))
                  + geom_line(aes(x = Time,
                                  y = Reference,
                                  colour = "Reference"))
                  + scale_colour_manual(name="",
                                        values = c("Reference"= "black",
                                                   "Test" = "blue"))
                  + labs(x = 'Time (h)',
                         y = 'Normalized AUC (%)')
                  + annotate('text', x = tail(Normalized$Time,1)*(3/4), y = 40,
                             label = bquote('AUC'['0-t']~italic(f)[2] == .(round(f2.AUC$f2,2))))
                  + theme_bw()
                  + theme(panel.background = element_blank(),
                          # panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          legend.position = "bottom"))

  }


  # Defined output
  if (details & plot) {

    suppressWarnings(print(Norm.plot))
    return(out)

  } else if (!details & plot) {

    suppressWarnings(print(Norm.plot))
    cat("AUC f2 Factor:", round(f2.AUC$f2,2))
    return(invisible(out))

  } else if (details & !plot) {

    return(out)

  } else if (!details & !plot) {

    cat("AUC f2 Factor:", round(f2.AUC$f2,2))
    return(invisible(out))

  }
}
