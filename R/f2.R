#' Similarity f2 Factor
#'
#' @description `f2()` calculates the \eqn{f_2} Similarity Factor, as proposed by Moore and Flanner in 1996.
#'
#' @param dta Dataframe with data from Test and Reference products.
#' @param Ref Name of the column with data from Reference product.
#' @param Test Name of the column with data from Test product.
#'
#' @return Dataframe with:
#'   * `n`: number of timepoints used for \eqn{f_2} computation.
#'   * `MSE`: mean squared error.
#'   * `Difference`: mean Test-to-Reference difference.
#'   * `f2`: Similarity \eqn{f_2} Factor
#'
#' @details The similarity \eqn{f_2} factor is a mathematical index widely used to compare dissolution
#' profiles, evaluating their similarity, using the percentage of drug dissolved per unit of time.
#'
#' The similarity \eqn{f_2} factor, proposed by Moore and Flanner in 1996, is derived from the
#' mean squared difference, and can be calculated as a function of the reciprocal of mean squared-root
#' transformation of the sum of square differences at all points:
#' \deqn{
#'   {f_2} = 50 \cdot \log \left( 100 \cdot \left[ 1 + \frac{1}{n} \sum_{t=1}^{t=n} \left(\overline{R}_t
#'   - \overline{T}_t\right)^2 \right]^{-0.5} \right)
#' }
#' where \eqn{f_2} is the similarity factor, \eqn{n} is the number of time points, and \eqn{{\overline{R}_t}}
#' and \eqn{{\overline{T}_t}} are the mean percentage of drug dissolved at time \eqn{t} after initiation of
#' the study, for Reference and Test products, respectively.
#'
#' The \eqn{f_2} similarity factor ranges from 0 (when \eqn{{\overline{R}_t} - {\overline{T}_t} = 100%}, at
#' all \eqn{t}) to 100 (when \eqn{{\overline{R}_t} - {\overline{T}_t} = 0%}, at all \eqn{t}).
#'
#' @author Sara Carolina Henriques
#'
#' @references
#' Henriques, S.C.; Albuquerque, J.; Paixão, P.; Almeida, L.; Silva, N.E. (2023).
#' Alternative Analysis Approaches for the Assessment of Pilot Bioavailability/Bioequivalence Studies.
#' *Pharmaceutics*. *15*(5), 1430.
#' [10.3390/pharmaceutics15051430](https://doi.org/10.3390/pharmaceutics15051430).
#'
#' Moore, J.W.; Flanner, H.H. (1996). Mathematical Comparison of Curves with an Emphasis on in Vitro
#' Dissolution Profiles. *Pharm. Technol*. *20*, 64–74.
#'
#'
#' @examples
#' f2(dta, 'R', 'T')
#'
#' @export
f2 <- function(dta, Ref = 'Reference', Test = 'Test') {

  # Rename columns
  names(dta)[names(dta) == Ref] <- "Reference"
  names(dta)[names(dta) == Test] <- "Test"

  # Number of timepoints
  n <- nrow(dta)-1

  # Mean Squared Error
  MSE <- (1/n)*(sum((dta$Reference-dta$Test)^2))

  # Mean Test-to-Reference Difference
  Difference <- ((MSE)^(0.5))

  # f2 Similarity Factor
  f2  <- 50*log10(100*((1+MSE)^(-0.5)))

  # Output
  return(data.frame(n = n, MSE = MSE, Difference = Difference, f2 = f2))
}
