#' Cmax Similarity f2 Factor
#'
#' @description `f2.Cmax()` calculates the \eqn{f_2} Similarity factor for \eqn{C_{\text{max}}},
#' from concentration data.
#'
#' @param dta Dataframe with concentration data from Test and Reference products.
#' @param Time Name of the column with time data.
#' @param Conc Name of the column with concentration data
#'   (only for stacked data, `Trt.cols = FALSE`).
#' @param Treatment Name of the column with treatment/formulation information
#'   (only for stacked data, `Trt.cols = FALSE`).
#' @param Reference Nomenclature in `dta` for the Reference product, defaults to `R`.
#' @param Test Nomenclature in `dta` for the Test product, defaults to `T`.
#' @param Trt.cols A logical value indicating whether treatment is presented in
#'   columns/pivoted (`TRUE`), or in rows/stacked (`FALSE`).
#' @param details A logical value indicating whether detailed results will be
#'   presented (`TRUE`), or if the function only returns the f2 factor
#'   (`FALSE`). Defaults to `FALSE`.
#' @param plot A logical value indicating whether graphical representation will be
#'   returned. Defaults to `TRUE`.
#'
#' @return Returns a list with the following elements:
#'   * `Raw Concentration Data`: Dataframe of input concentration data, with
#'     treatment information in columns (pivoted).
#'   * `Reference Tmax`: Vector of Reference \eqn{t_{\text{max}}}
#'   * `Reference Cmax`: Vector of Reference \eqn{C_{\text{max}}}
#'   * `Normalized Concentrations`: dataframe with normalize Test and Reference
#'     concentrations over time, by Reference \eqn{C_{\text{max}}}, until Reference \eqn{t_{\text{max}}}
#'   * `Cmax f2 Factor`: dataframe from [f2()] function.
#'
#'
#' @details Henriques *et al* (2023) proposed to used the \eqn{f_2} factor to assess the similarity on
#' the rate of drug absorption by normalizing Test and Reference mean concentration-time profiles to the
#' maximum plasma concentration (\eqn{C_{\text{max}}}) derived from the mean Reference profile, until
#' Reference Cmax is observed (Reference \eqn{t_{\text{max}}}):
#' \deqn{
#'   C_{t}^{N} = 100 \cdot \frac{\overline{C}_{t}}{C_{\text{max R}}}, \text{ where } 0 \leq t \leq
#'   t_{\text{max R}}
#' }
#' where \eqn{C_{t}^{N}} is the normalized concentration at time \eqn{t}, \eqn{\overline{C}_{t}}
#' is the mean (Test or Reference) concentration at time \eqn{t}, \eqn{C_{\text{max R}}} is the
#' \eqn{C_{\text{max}}} of the Reference mean concentration-time profile, and \eqn{t_{\text{max R}}}
#' the time of observation of \eqn{C_{\text{max R}}}. The similarity \eqn{f_2} factor is calculated as
#' \deqn{
#'   C_{\text{max}} {f_2} = 50 \cdot \log \left( 100 \cdot \left[ 1 + \frac{1}{n} \sum_{t=1}^{t=n}
#'   \left(\overline{R}_t^N - \overline{T}_t^N\right)^2 \right]^{-0.5} \right)
#' }
#' where \eqn{C_{\text{max}} {f_2}} is the similarity factor calculated for \eqn{C_{\text{max}}},
#' \eqn{n} is the number of time points until Reference \eqn{t_{\text{max}}}, and \eqn{{\overline{R}_t^N}}
#' and \eqn{{\overline{T}_t^N}} are the are the normalized concentration at time \eqn{t} for Reference and
#' Test products, respectively.
#'
#' @author Sara Carolina Henriques
#'
#' @references
#' Henriques, S.C.; Albuquerque, J.; Paixão, P.; Almeida, L.; Silva, N.E. (2023).
#' Alternative Analysis Approaches for the Assessment of Pilot Bioavailability/Bioequivalence Studies.
#' *Pharmaceutics*. *15*(5), 1430.
#' [10.3390/pharmaceutics15051430](https://doi.org/10.3390/pharmaceutics15051430).
#'
#'
#'
#' @examples
#' # Calculate f2 Factor for Cmax
#' # when treatment data is pivoted:
#' f2.Cmax(dta, Time = 'Time', Reference = 'R', Test = 'T'
#'
#' # when treatment data is stacked:
#' f2.Cmax(dta, Time = 'Time', Conc = 'Conc',
#'         Treatment = 'Treatment', Reference = 'R', Test = 'T')
#'
#' @export
f2.Cmax <- function(dta, Time = 'Time', Conc = 'Conc',
                    Treatment = 'Treatment', Reference = 'R', Test = 'T',
                    Trt.cols = TRUE, details = FALSE, plot = TRUE) {

  # If data is pivoted
  if (Trt.cols) {

    # Rename columns
    cols <- c(Time, Reference, Test)
    dta <- dta[,cols]
    names(dta)[names(dta) == Time] <- "Time"
    names(dta)[names(dta) == Reference] <- "Reference"
    names(dta)[names(dta) == Test] <- "Test"

  } else { # If data is stacked

    # Rename columns
    cols <- c(Treatment, Time, Conc)
    dta <- dta[,cols]
    names(dta)[names(dta) == Treatment] <- "Treatment"
    names(dta)[names(dta) == Time] <- "Time"
    names(dta)[names(dta) == Conc] <- "Conc"

    # Data from Test
    dta_T <- dta[dta$Treatment==Test,c('Time','Conc')]
    names(dta_T)[names(dta_T) == 'Conc'] <- 'Test'

    # Data from Reference
    dta_R <- dta[dta$Treatment==Reference,c('Time','Conc')]
    names(dta_R)[names(dta_R) == 'Conc'] <- 'Reference'

    dta <- merge(dta_R,dta_T)
  }


  # Identify Reference maximum observed concentration (Cmax)
  Cmax.R <- Cmax(dta, Conc = 'Reference')
  # and time of occurrence of Reference Cmax (Tmax)
  Tmax.R <- Tmax(dta, Conc = 'Reference')



  # Normalize Test and Reference concentrations over time by Reference Cmax,
  # until Reference Tmax.
  # Results presented as percentage (%)
  Normalized <- NULL
  for (t in 1:which(dta$Time == Tmax.R)) {
    Normalized.t <- data.frame(
      Time      = dta$Time[t],
      Reference = 100*(dta$Reference[t]/Cmax.R),
      Test      = 100*(dta$Test[t]/Cmax.R)
    )
    Normalized.t <- within(Normalized.t, {
      Sq.Diff    = (Reference - Test)^2
      Difference = Reference - Test
    })
    Normalized <- rbind(Normalized,Normalized.t)
  }


  # Calculate f2 Factor for Cmax
  f2.Cmax <- f2(Normalized)


  # Output
  out <- list(
    'Raw Concentration Data' = dta,
    'Reference Tmax' = Tmax.R,
    'Reference Cmax' = Cmax.R,
    'Normalized Concentrations' = Normalized,
    'Cmax f2 Factor' = f2.Cmax
  )


  if (plot) {

    # Plot of normalized concentration over time, until Reference Tmax
    require(ggplot2)
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
                         y = 'Normalized Concentration (%)')
                  + annotate('text', x = tail(Normalized$Time,1)*(3/4), y = 40,
                             label = bquote(C['max']~italic(f)[2] == .(round(f2.Cmax$f2,2))))
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
    cat("Cmax f2 Factor:", round(f2.Cmax$f2,2))
    return(invisible(out))

  } else if (details & !plot) {

    return(out)

  } else if (!details & !plot) {

    cat("Cmax f2 Factor:", round(f2.Cmax$f2,2))
    return(invisible(out))

  }
}
