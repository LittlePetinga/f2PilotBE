f2PilotBE
================

- [Background](#background)
  - [The Similarity f2 Factor](#the-similarity-f2-factor)
  - [Cmax Similarity f2 Factor](#cmax-similarity-f2-factor)
  - [Bioequivalence Evaluation of Pilot
    Studies](#bioequivalence-evaluation-of-pilot-studies)
- [The f2PilotBE Package](#the-f2pilotbe-package)
  - [Installation](#installation)
    - [From GitHub](#from-github)
  - [Examples](#examples)
    - [Load Package](#load-package)
    - [Calculate Geometric Mean
      Concentration](#calculate-geometric-mean-concentration)
    - [Calculate Cmax Similarity f2
      Factor](#calculate-cmax-similarity-f2-factor)
- [References](#references)

<!-- badges: start -->
<!-- ![GitHub Actions Status](https://github.com/LittlePetinga/f2PilotBE/workflows/R-CMD-check/badge.svg) -->

[![License:
GPL-3.0-or-later](https://img.shields.io/badge/License-GPL--3.0--or--later-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.html)
[![repo](https://img.shields.io/badge/Repo%20Since-Nov%202023-yellow)](https://github.com/LittlePetinga/f2PilotBE)
[![Version](https://img.shields.io/badge/Version-0.0.1-blue)](https://github.com/LittlePetinga/f2PilotBE)
[![Last
Update](https://img.shields.io/badge/Last%20Update-11%20Nov%202023-green)](https://github.com/LittlePetinga/f2PilotBE)
[![Project Status:
Active](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Paper f2
BE](https://img.shields.io/badge/Scientific%20Paper-Part%201-9D1781)](https://doi.org/10.3390/pharmaceutics15051430)
[![Paper f2
BE](https://img.shields.io/badge/Scientific%20Paper-Part%202-9D172A)](https://doi.org/10.3390/pharmaceutics15102498)
[![R
badge](https://img.shields.io/badge/Build%20with-♥%20and%20R-E0218A)](https://cdn.dribbble.com/users/6620596/screenshots/14792345/media/af61fa935b055891cb800a9e41ebb747.gif)
<!-- [![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/f2PilotBE)](https://cran.r-project.org/package=f2PilotBE) 
[![f2PilotBE badge](https://img.shields.io/badge/f2PilotBE-Ready%20to%20Use-brightgreen)](https://github.com/LittlePetinga/f2PilotBE) -->
<!-- badges: end -->

<!-- Hex Sticker -->

<img src="man/figures/logo.png" align="right" height="240">

The `f2PilotBE` package is designed to calculate the similarity factor
*f*<sub>2</sub> based on pharmacokinetic profiles from pilot
bioavailability/bioequivalence (BA/BE) studies, as an alternative
approach to assess the potential of bioequivalence of a Test product in
comparison to a Reference product in terms of maximum observed
concentration (C<sub>max</sub>) and area under the concentration-time
curve (AUC).

This package is based on the publications of Henriques *et al.* (2023)
\[<a href="#references">1</a>,<a href="#references">2</a>\], which
propose the use of the geometric mean (G<sub>mean</sub>) *f*<sub>2</sub>
factor, for the comparison of the absorption rate (given by
C<sub>max</sub>) of Test and Reference formulations. According to the
articles, G<sub>mean</sub> *f*<sub>2</sub> factor using a cut-off of 35
showed a good relationship between avoiding type I and type II errors
\[<a href="#references">1</a>,<a href="#references">2</a>\].

## Background

### The Similarity f2 Factor

The similarity *f*<sub>2</sub> factor is a mathematical index widely
used to compare dissolution profiles, evaluating their similarity, using
the percentage of drug dissolved per unit of time.

The similarity *f*<sub>2</sub> factor, proposed by Moore and Flanner in
1996 \[<a href="#references">3</a>\], is derived from the mean squared
difference, and can be calculated as a function of the reciprocal of
mean squared-root transformation of the sum of square differences at all
points:

$$f_{2}=50\cdot\log\biggl(100\cdot\biggl[1+\frac{1}{n}\sum_{t=1}^{t=n}{(\bar{R}_{t}-\bar{T}_{t})^2} \biggr]^{-0.5}\biggr)$$

where $f_{2}$ is the similarity factor, $n$ is the number of time
points, and $\bar{R}_{t}$ and $\bar{T}_{t}$ are the mean percentage of
drug dissolved at time $t$, for Reference and Test products respectively
\[<a href="#references">3</a>\].

The *f*<sub>2</sub> similarity factor ranges from 0 (when
$\bar{R}_{t}-\bar{T}_{t}=100\%$, at all $t$) to 100 (when
$\bar{R}_{t}-\bar{T}_{t}= 0\%$, at all $t$)
\[<a href="#references">3</a>\].

<a href="#fig1">Figure 1</a> (from Henriques *et al.* (2023)
\[<a href="#references">1</a>\]) presents the distribution of
*f*<sub>2</sub> similarity factor as a function of mean difference. Form
the *f*<sub>2</sub> equation, an average difference of
<span style="color: red;">10%</span>,
<span style="color: green;">15%</span>, and
<span style="color: blue;">20%</span> from all measured time points
results in a *f*<sub>2</sub> value of
<span style="color: red;">50</span> (red dotted lines),
<span style="color: green;">41</span> (green dotted lines) and
<span style="color: blue;">35</span> (blue dotted lines), respectively.

<!-- Figure from article (Part 1) -->
<figure>
<img src="https://pub.mdpi-res.com/pharmaceutics/pharmaceutics-15-01430/article_deploy/html/images/pharmaceutics-15-01430-g002.png?1683533911" alt="Image Description" height="400">
<figcaption>
Figure 1. <a name="fig1"></a> Distribution of ƒ<sub>2</sub> similarity
factor as a function of mean difference. ƒ<sub>2</sub> similarity factor
is derived from the mean squared difference and can be calculated as a
function of the reciprocal of the mean squared-root transformation of
the sum of square differences at all points. An average difference of
<span style="color: red;">10%</span>,
<span style="color: green;">15%</span>, and
<span style="color: blue;">20%</span> from all measured time points
results in a ƒ<sub>2</sub> value of <span style="color: red;">50</span>
(red dotted lines), <span style="color: green;">41</span> (green dotted
lines) and <span style="color: blue;">35</span> (blue dotted lines),
respectively.
</figcaption>
</figure>
<!-- Code for plot of distribution of ƒ~2~ similarity factor as a function of mean difference -->
<!-- Not run -->

As proposed by Henriques *et al.* (2023), the concept of similarity
factor *f*<sub>2</sub> can be applied as an alternative to the average
bioequivalence analysis, for pilot BA/BE studies
\[<a href="#references">1</a>,<a href="#references">2</a>\].

### Cmax Similarity f2 Factor

*f*<sub>2</sub> can be used to assess the similarity on the rate of drug
absorption by normalizing Test and Reference mean concentration-time
profiles to the maximum plasma concentration (C<sub>max</sub>) derived
from the mean Reference profile, until Reference C<sub>max</sub> is
observed (Reference t<sub>max</sub>)
\[<a href="#references">1</a>,<a href="#references">2</a>\]:

$$C_{t}^{N}=100\cdot\frac{\bar{C}_{t}}{C_{max,R}}$$, where
$0 \le t \le t_{max,R}$,

where $C_{t}^{N}$ is the normalized concentration at time $t$,
$\bar{C}_{t}$ is the mean (Test or Reference) concentration at time $t$,
$C_{max,R}$ is the C<sub>max</sub> of the Reference mean
concentration-time profile, and $t_{max,R}$ the time of observation of
$C_{max,R}$. The similarity *f*<sub>2</sub> factor is calculated as

$$C_{max}f_{2}=50\cdot\ log \biggl(100\cdot \biggl[1+\frac{1}{n} \sum_{t=1}^{t=n}{(R_{t}^{N} - T_{t}^{N})^2} \biggr]^{-0.5} \biggr)$$

where $n$ is the number of time points until Reference $t_{max}$, and
$R_{t}^{N}$ and $T_{t}^{N}$ are the normalized concentration at time
$t$, for Reference and Test products respectively.

### Bioequivalence Evaluation of Pilot Studies

For the planning of pilot BA/BE studies, a decision tree is proposed
(<a href="#fig2">Figure 2</a>, from Henriques *et al.* (2023)
\[<a href="#references">2</a>\]).

For drug products with a known Intra-Subject Coefficient of Variation
(ISCV%) below 20%, the authors propose the estimation of the sample size
for a pilot study assuming a Test-to-Reference Geometric Least Square
Mean Ratio (GMR) of 100%, a power of 80%, and an α of 0.05
\[<a href="#references">1</a>,<a href="#references">2</a>\]. However,
for cases of higher ISCV% or unknown variability, it is propose the use
of a fixed sample size of **20 subjects**, as the use of higher sample
sizes has not shown to increase the study power meaningfully, but was
sufficient to avoid substantial type I errors
\[<a href="#references">2</a>\].

Regarding the analysis of data from pilot studies, the authors propose
to initially analyze the data using the average bioequivalence approach.
For the case in which the calculated GMR and the corresponding 90% CI
are not within \[80.00–125.00\]%, the alternative G<sub>mean</sub>
*f*<sub>2</sub> factor method should be used with a cut off of 35, as it
was shown to be a valuable indicator of the potentiality of the Test
formulation to be bioequivalent in terms of C<sub>max</sub>:

1.  If the *f*<sub>2</sub> factor is above or equal to **35**
    (corresponding to a difference of 20% between Test and Reference
    concentration–time profiles until the Reference t<sub>max</sub>),
    the confidence to proceed to a pivotal study is higher than 90% when
    ISCV% is lower or equal to 20%; the confidence is higher than 80%
    when ISCV% is within 20% and 30%; and the confidence is higher than
    60% when ISCV% is higher than 40% \[<a href="#references">2</a>\].

2.  If the *f*<sub>2</sub> factor is above or equal to **41**
    (corresponding to a difference of 15% between Test and Reference
    concentration–time profiles until the Reference t<sub>max</sub>),
    the confidence to proceed to a pivotal study is higher than 90% for
    ISCV% until 40%, and higher than 80% for ISCV% within 50% to 60%
    \[<a href="#references">2</a>\].

3.  If the *f*<sub>2</sub> factor is above or equal to **50**
    (corresponding to a difference of 10% between Test and Reference
    concentration–time profiles until the Reference t<sub>max</sub>),
    the probability of the Test product to be truly bioequivalent to the
    Reference product in terms of C<sub>max</sub>, i.e., the confidence
    to proceed to a pivotal study, is higher than 90%, irrespective of
    the ISCV% \[<a href="#references">2</a>\].

<figure>
<img src="https://www.mdpi.com/pharmaceutics/pharmaceutics-15-02498/article_deploy/html/images/pharmaceutics-15-02498-g013.png" alt="Image Description" height="800">
<figcaption>
Figure 2. <a name="fig2"></a> Proposed decision tree for planning and
analysis of pilot BA/BE studies.
</figcaption>
</figure>

## The f2PilotBE Package

The `f2PilotBE` is equipped with the following functions to aid in the
calculation of the similarity *f*<sub>2</sub> factor:

- `AUC()` – calculates the cumulative area under the concentration-time
  curve (AUC) over time.

- `AUClast()` – calculates area under the concentration-time curve (AUC)
  from time zero to the time of the last measurable concentration
  (AUC<sub>0-t</sub>).

- `Cmax()` – calculates the maximum observed concentration post-dose
  (C<sub>max</sub>) directly obtained from the observed
  concentration-time profile.

- `f2.AUC()` – calculates the *f*<sub>2</sub> similarity factor for AUC,
  from concentration data.

- `f2.Cmax()` – calculates the *f*<sub>2</sub> similarity factor for
  C<sub>max</sub>, from concentration data.

- `f2()` – calculates the *f*<sub>2</sub> similarity factor, as proposed
  by Moore and Flanner in 1996.

- `geomean()` – calculates geometric mean, i.e., the *N*th root of the
  product of the *N* observations, equivalent to `exp(mean(log(x)))`.

- `Tmax()` – calculates the time of the C<sub>max</sub>.

### Installation

#### From GitHub

To install the development version of `f2PilotBE`, start by installing
the `devtools` package from CRAN:

``` r
install.packages("devtools")
```

Then the development version of `f2PilotBE` can be installed from GitHub
as:

``` r
library(devtools)
install_github("LittlePetinga/f2PilotBE")
```

### Examples

#### Load Package

Following installation, load the `f2PilotBE` package and the `ggplot2`
package.

``` r
# Load packages
library("ggplot2")      # required for plotting
library("f2PilotBE")
```

#### Calculate Geometric Mean Concentration

As proposed by the authors
\[<a href="#references">1</a>,<a href="#references">2</a>\], the
*f*<sub>2</sub> similarity factor should be calculated from the
arithmetic (A<sub>mean</sub>) or geometric (G<sub>mean</sub>)
concentration-time profiles.

The `f2PilotBE` provides a function (`geomean()`) for the calculation of
the G<sub>mean</sub>. For the calculation of the G<sub>mean</sub>
profile for Test and Reference products, first, import the individual
concentration-time data, and calculate the G<sub>mean</sub>
concentration-time profiles for each Treatment, as follows.

``` r
# Import individual concentration data
dta_id <- read.csv('dta_id.csv')


# Calculate the geometric mean for each Treatment, by time point
dta <- data.frame(t(tapply(dta_id$Concentration, 
                           dta_id[,c('Treatment','Time')], 
                           geomean)))
dta$Time <- as.numeric(row.names(dta))
row.names(dta) <- c()
```

The code above, will apply (`tapply`) the `geomean` function for each
given treatment, and each timepoint.

If the user prefers to use the A<sub>mean</sub>, instead of the
G<sub>mean</sub> concentration-time data, simply replace the function
`geomean`, with the function `mean` from base `R`.

#### Calculate Cmax Similarity f2 Factor

The `f2.Cmax()` function allows to calculate the *f*<sub>2</sub>
similarity factor for C<sub>max</sub>, from Test and Reference mean
concentration-time profiles.

This function allows some flexibility on the structure of the mean
concentration-time dataframe. Treatment (Test and Reference) information
should either be provided in columns (pivoted) or rows (stacked). For
both cases a column with `Time` information should always be provided:

1.  For the case where Treatment information is provided in **columns**
    (pivoted), as in the following example:

``` r
# Create mean concentration-time data for Test and Reference product
# where treatment information is pivoted
dta_piv <- data.frame(
  Time = c(0, 0.25, 0.5, 0.75, 1, 1.5, 1.75, 2, 2.25, 2.5,
           2.75, 3, 3.25, 3.5, 3.75, 4, 6, 8, 12, 24),
  Ref = c(0.00, 221.23, 377.19, 494.73, 555.74,
          623.86, 615.45, 663.38, 660.29, 621.71,
          650.33, 622.28, 626.72, 574.94, 610.51,
          554.02, 409.14, 299.76, 162.85, 27.01),
  Test = c(0.00, 149.24, 253.05, 354.49, 412.49,
           530.07, 539.68, 566.30, 573.54, 598.33,
           612.63, 567.48, 561.10, 564.47, 541.50,
           536.92, 440.32, 338.78, 185.03, 31.13)
 )
head(dta_piv)
#>   Time    Ref   Test
#> 1 0.00   0.00   0.00
#> 2 0.25 221.23 149.24
#> 3 0.50 377.19 253.05
#> 4 0.75 494.73 354.49
#> 5 1.00 555.74 412.49
#> 6 1.50 623.86 530.07
```

The *f*<sub>2</sub> similarity factor for C<sub>max</sub> can be
calculated from the mean concentration-time profiles, applying the
`f2.Cmax()` function as follows:

``` r
# Calculate f2 Factor for Cmax when treatment data is pivoted:
f2.Cmax(dta_piv, Time = "Time", Ref = "Ref", Test = "Test", plot = FALSE)
#> Cmax f2 Factor: 38.97
```

As demonstrated above, for the case when treatment data is pivoted, the
`f2.Cmax()` function requires the indication of: (1) the dataframe with
mean concentration-time data (in this case `dta = dta_piv`), (2) the
name of the column in the dataframe with time information (in this case
`Time = "Time"`), (3) the name of the column in the dataframe with
concentration information for the Reference product (in this case
`Ref = "Reference"`), and (4) the name of the column in the dataframe
with concentration information for the Test product (in this case
`Test = "Test"`).

Please note that `Trt.cols` is not required, as it defaults to `TRUE`
(i.e, logical value indicating whether treatment is presented in
columns/pivoted).

2.  For the case where Treatment information is provided in **rows**
    (stacked), as in the following example:

``` r
# Create mean concentration-time data for Test and Reference product
# where treatment information is stacked
dta_stk <- data.frame(
  Time = rep(c(0, 0.25, 0.5, 0.75, 1, 1.5, 1.75, 2, 2.25, 2.5,
               2.75, 3, 3.25, 3.5, 3.75, 4, 6, 8, 12, 24),2), 
  Trt = c(rep('R',20), rep('T',20)),
  Conc = c(c(0.00, 221.23, 377.19, 494.73, 555.74,
             623.86, 615.45, 663.38, 660.29, 621.71,
             650.33, 622.28, 626.72, 574.94, 610.51,
             554.02, 409.14, 299.76, 162.85, 27.01),
           c(0.00, 149.24, 253.05, 354.49, 412.49,
             530.07, 539.68, 566.30, 573.54, 598.33,
             612.63, 567.48, 561.10, 564.47, 541.50,
             536.92, 440.32, 338.78, 185.03, 31.13))
)
head(dta_stk)
#>   Time Trt   Conc
#> 1 0.00   R   0.00
#> 2 0.25   R 221.23
#> 3 0.50   R 377.19
#> 4 0.75   R 494.73
#> 5 1.00   R 555.74
#> 6 1.50   R 623.86
tail(dta_stk)
#>     Time Trt   Conc
#> 35  3.75   T 541.50
#> 36  4.00   T 536.92
#> 37  6.00   T 440.32
#> 38  8.00   T 338.78
#> 39 12.00   T 185.03
#> 40 24.00   T  31.13
```

The *f*<sub>2</sub> similarity factor for C<sub>max</sub> can be
calculated from the mean concentration-time profiles, applying the
`f2.Cmax()` function as follows:

``` r
# Calculate f2 Factor for Cmax when treatment data is stacked:
f2.Cmax(dta_stk, Time = "Time", Conc = "Conc",
        Trt = "Trt", Ref = "R", Test = "T",
        Trt.cols = FALSE, plot = FALSE)
#> Cmax f2 Factor: 38.97
```

As demonstrated above, for the case when treatment data is stacked, the
`f2.Cmax()` function requires the indication of: (1) the dataframe with
mean concentration-time data (in this case `dta = dta_stk`), (2) the
name of the column in the dataframe with time information (in this case
`Time = "Time"`), (3) the name of the column in the dataframe with
concentration information (in this case `Conc = "Conc"`), (4) the name
of the column in the dataframe with treatment information (in this case
`Trt = "Trt"`), (5) the nomenclature in the dataframe for the Reference
product (in this case `Ref = "R"`), (6) the nomenclature in the
dataframe for the Test product (in this case `Test = "T"`), and (7)
indication that treatment is presented in rows/stacked
(`Trt.cols = FALSE`).

For both cases, a graphical representation of the normalized Test and
Reference concentrations over time, by Reference C<sub>max</sub>, until
Reference t<sub>max</sub>, can be generated, by simply turning `plot` to
`TRUE`:

``` r
# Calculate f2 Factor for Cmax, and plot normalized concentration by Reference Cmax
f2.Cmax(dta_piv, Time = "Time", Ref = "Ref", Test = "Test", plot = TRUE)
```

<img src="man/figures/README-unnamed-chunk-12-1.png" width="70%" />

    #> Cmax f2 Factor: 38.97

By default, `f2.Cmax()` plots the Reference product in `black`, and the
Test product in <span style="color: blue;">`blue`</span>, nevertheless,
the user can personalize the output plot with different colours
according to their preference, by using `col.R` and `col.T`, for
Reference and Test product respectively:

``` r
# Calculate f2 Factor for Cmax, and plot normalized concentration by Reference Cmax
f2.Cmax(dta_piv, Time = "Time", Ref = "Ref", Test = "Test", 
        plot = TRUE, col.R = "darkblue", col.T = "red")
```

<img src="man/figures/README-unnamed-chunk-13-1.png" width="70%" />

    #> Cmax f2 Factor: 38.97

## References

<a name="references"></a>

1.  Henriques, S.C.; Albuquerque, J.; Paixão, P.; Almeida, L.; Silva,
    N.E. (2023). Alternative Analysis Approaches for the Assessment of
    Pilot Bioavailability/Bioequivalence Studies. *Pharmaceutics*.
    *15*(5), 1430.
    [10.3390/pharmaceutics15051430](https://doi.org/10.3390/pharmaceutics15051430).

2.  Henriques, S.C.; Paixão, P.; Almeida, L.; Silva, N.E. (2023).
    Predictive Potential of C<sub>max</sub> Bioequivalence in Pilot
    Bioavailability/Bioequivalence Studies, through the Alternative
    *f*<sub>2</sub> Similarity Factor Method. *Pharmaceutics*.
    *15*(10), 2498.
    [10.3390/pharmaceutics15102498](https://doi.org/10.3390/pharmaceutics15102498).

3.  Moore, J.W.; Flanner, H.H. (1996). Mathematical Comparison of Curves
    with an Emphasis on in Vitro Dissolution Profiles. *Pharm.
    Technol.*. *20*, 64–74.
