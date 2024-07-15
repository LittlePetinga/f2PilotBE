f2PilotBE
================

- [1. Background](#1-background)
  - [1.1. The Similarity ƒ2 Factor](#11-the-similarity-ƒ2-factor)
  - [1.2. Cmax Similarity ƒ2 Factor](#12-cmax-similarity-ƒ2-factor)
  - [1.3. Bioequivalence Evaluation of Pilot
    Studies](#13-bioequivalence-evaluation-of-pilot-studies)
- [2. The f2PilotBE Package](#2-the-f2pilotbe-package)
  - [2.1. Installation](#21-installation)
    - [2.1.1. From GitHub](#211-from-github)
  - [2.2. Examples](#22-examples)
    - [2.2.1. Load Package](#221-load-package)
    - [2.2.2. Calculate Geometric Mean
      Concentration](#222-calculate-geometric-mean-concentration)
    - [2.2.3. Calculate Cmax Similarity ƒ2
      Factor](#223-calculate-cmax-similarity-ƒ2-factor)
- [References](#references)

<!-- badges: start -->
<!-- [![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/f2PilotBE)](https://cran.r-project.org/package=f2PilotBE) -->
<!-- [![Codecov test -->
<!-- coverage](https://codecov.io/gh/saracarolinahenriques/f2PilotBE/branch/master/graph/badge.svg)](https://app.codecov.io/gh/saracarolinahenriques/f2PilotBE?branch=master) -->

[![f2PilotBE
badge](https://img.shields.io/badge/f2PilotBE-Ready%20to%20Use-brightgreen)](https://github.com/saracarolinahenriques/f2PilotBE)
[![License:
GPL-3.0-or-later](https://img.shields.io/badge/License-GPL--3.0--or--later-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.html)
[![repo](https://img.shields.io/badge/Repo%20Since-Jan%202024-ffcc00)](https://github.com/saracarolinahenriques/f2PilotBE)
[![Version](https://img.shields.io/badge/Version-0.4.0-00b5a3)](https://github.com/saracarolinahenriques/f2PilotBE)
[![Last
Update](https://img.shields.io/badge/Last%20Update-15%20Jul%202024-3cb371)](https://github.com/saracarolinahenriques/f2PilotBE)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check](https://github.com/saracarolinahenriques/f2PilotBE/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/saracarolinahenriques/f2PilotBE/actions/workflows/R-CMD-check.yaml)
[![Paper f2
BE](https://img.shields.io/badge/Scientific%20Paper-Part%201-d93467)](https://doi.org/10.3390/pharmaceutics15051430)
[![Paper f2
BE](https://img.shields.io/badge/Scientific%20Paper-Part%202-0074cc)](https://doi.org/10.3390/pharmaceutics15102498)
[![ModelBio](https://img.shields.io/badge/Presentation-ModelBio%202024-8a4bd9)](https://www.researchgate.net/publication/378769367_Pilot_BioavailabilityBioequivalence_Studies_Revisited_With_An_Alternative_f2_Factor_Approach)
[![PAGE](https://img.shields.io/badge/Presentation-PAGE%202024-C80815)](https://www.page-meeting.org/?abstract=11251)
[![R
badge](https://img.shields.io/badge/Build%20with-♥%20and%20R-E0218A)](https://cdn.dribbble.com/users/6620596/screenshots/14792345/media/af61fa935b055891cb800a9e41ebb747.gif)

<!-- badges: end -->
<!-- Hex Sticker -->

<img src="man/figures/logo.png" align="right" height="240"/>

The `f2PilotBE` package is designed to calculate the similarity factor
*ƒ*<sub>2</sub> based on pharmacokinetic profiles from pilot
bioavailability/bioequivalence (BA/BE) studies, as an alternative
approach to assess the potential of bioequivalence of a Test product in
comparison to a Reference product in terms of maximum observed
concentration (C<sub>max</sub>) and area under the concentration-time
curve (AUC).

This package is based on the publications of Henriques *et al.* (2023)
\[<a href="#references">1</a>,<a href="#references">2</a>\], which
propose the use of the geometric mean (G<sub>mean</sub>) *ƒ*<sub>2</sub>
factor, for the comparison of the absorption rate (given by
C<sub>max</sub>) of Test and Reference formulations. According to the
articles, G<sub>mean</sub> *ƒ*<sub>2</sub> factor using a cut-off of 35
showed a good relationship between avoiding type I and type II errors
\[<a href="#references">1</a>,<a href="#references">2</a>\].

## 1. Background

### 1.1. The Similarity ƒ2 Factor

The similarity *ƒ*<sub>2</sub> factor is a mathematical index widely
used to compare dissolution profiles, evaluating their similarity, using
the percentage of drug dissolved per unit of time.

The similarity *ƒ*<sub>2</sub> factor, proposed by Moore and Flanner in
1996 \[<a href="#references">3</a>\], is derived from the mean squared
difference, and can be calculated as a function of the reciprocal of
mean squared-root transformation of the sum of square differences at all
points:

$$f_{2}=50\cdot\log\biggl(100\cdot\biggl[1+\frac{1}{n}\sum_{t=1}^{t=n}{(\bar{R}_{t}-\bar{T}_{t})^2} \biggr]^{-0.5}\biggr)$$

where $f_{2}$ is the similarity factor, $n$ is the number of time
points, and $\bar{R}_{t}$ and $\bar{T}_{t}$ are the mean percentage of
drug dissolved at time $t$, for Reference and Test products respectively
\[<a href="#references">3</a>\].

The *ƒ*<sub>2</sub> similarity factor ranges from 0 (when
$\bar{R}_{t}-\bar{T}_{t}=100\%$, at all $t$) to 100 (when
$\bar{R}_{t}-\bar{T}_{t}= 0\%$, at all $t$)
\[<a href="#references">3</a>\].

<a href="#fig1">Figure 1</a> (from Henriques *et al.* (2023)
\[<a href="#references">1</a>\]) presents the distribution of
*ƒ*<sub>2</sub> similarity factor as a function of mean difference. Form
the *ƒ*<sub>2</sub> equation, an average difference of
<span style="color: red;">10%</span>,
<span style="color: green;">15%</span>, and
<span style="color: blue;">20%</span> from all measured time points
results in a *ƒ*<sub>2</sub> value of
<span style="color: red;">50</span> (red dotted lines),
<span style="color: green;">41</span> (green dotted lines) and
<span style="color: blue;">35</span> (blue dotted lines), respectively.

<!-- Figure from article (Part 1) -->
<figure>

<img src="https://pub.mdpi-res.com/pharmaceutics/pharmaceutics-15-01430/article_deploy/html/images/pharmaceutics-15-01430-g002.png?1683533911" alt="Image Description" height="400"/>

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
respectively. From Henriques, S.C. et al. (2023). Pharmaceutics, 15(5),
1430 (DOI: 10.3390/pharmaceutics15051430).
</figcaption>
</figure>
<!-- Code to plot distribution of ƒ~2~ similarity factor as a function of mean difference -->
<!-- Not run -->

As proposed by Henriques *et al.* (2023), the concept of similarity
factor *ƒ*<sub>2</sub> can be applied as an alternative to the average
bioequivalence analysis, for pilot BA/BE studies
\[<a href="#references">1</a>,<a href="#references">2</a>\].

### 1.2. Cmax Similarity ƒ2 Factor

*ƒ*<sub>2</sub> can be used to assess the similarity on the rate of drug
absorption by normalizing Test and Reference mean concentration-time
profiles to the maximum plasma concentration (C<sub>max</sub>) derived
from the mean Reference profile, until Reference C<sub>max</sub> is
observed (Reference t<sub>max</sub>)
\[<a href="#references">1</a>,<a href="#references">2</a>\]
(<a href="#fig2">Figure 2</a>):

<!-- Figure Estimation Cmax f2 Factor -->
<figure>
<img src="man/figures/Cmaxf2_1.png"/>
<figcaption>
Figure 2. <a name="fig2"></a> Estimation of ƒ<sub>2</sub> similarity
factor from mean concentration-time profiles.
</figcaption>
</figure>
<!-- $$C_{t}^{N}=100\cdot\frac{\bar{C}_{t}}{C_{max,R}}$$, where -->
<!-- $0 \le t \le t_{max,R}$, -->

where $C_{t}^{N}$ is the normalized concentration at time $t$,
$\bar{C}_{t}$ is the mean (Test or Reference) concentration at time $t$,
$C_{max,R}$ is the C<sub>max</sub> of the Reference mean
concentration-time profile, and $t_{max,R}$ the time of observation of
$C_{max,R}$.

The similarity *ƒ*<sub>2</sub> factor is calculated as

$$C_{max}f_{2}=50\cdot\ log \biggl(100\cdot \biggl[1+\frac{1}{n} \sum_{t=1}^{t=n}{(R_{t}^{N} - T_{t}^{N})^2} \biggr]^{-0.5} \biggr)$$

where $n$ is the number of time points until Reference $t_{max}$, and
$R_{t}^{N}$ and $T_{t}^{N}$ are the normalized concentration at time
$t$, for Reference and Test products respectively.

### 1.3. Bioequivalence Evaluation of Pilot Studies

For the planning of pilot BA/BE studies, a decision tree is proposed
(<a href="#fig3">Figure 3</a>, from Henriques *et al.* (2023)
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
*ƒ*<sub>2</sub> factor method should be used with a cut off of 35, as it
was shown to be a valuable indicator of the potentiality of the Test
formulation to be bioequivalent in terms of C<sub>max</sub>:

1.  If the *ƒ*<sub>2</sub> factor is above or equal to **35**
    (corresponding to a difference of 20% between Test and Reference
    concentration–time profiles until the Reference t<sub>max</sub>),
    the confidence to proceed to a pivotal study is higher than 90% when
    ISCV% is lower or equal to 20%; the confidence is higher than 80%
    when ISCV% is within 20% and 30%; and the confidence is higher than
    60% when ISCV% is higher than 40% \[<a href="#references">2</a>\].

2.  If the *ƒ*<sub>2</sub> factor is above or equal to **41**
    (corresponding to a difference of 15% between Test and Reference
    concentration–time profiles until the Reference t<sub>max</sub>),
    the confidence to proceed to a pivotal study is higher than 90% for
    ISCV% until 40%, and higher than 80% for ISCV% within 50% to 60%
    \[<a href="#references">2</a>\].

3.  If the *ƒ*<sub>2</sub> factor is above or equal to **50**
    (corresponding to a difference of 10% between Test and Reference
    concentration–time profiles until the Reference t<sub>max</sub>),
    the probability of the Test product to be truly bioequivalent to the
    Reference product in terms of C<sub>max</sub>, i.e., the confidence
    to proceed to a pivotal study, is higher than 90%, irrespective of
    the ISCV% \[<a href="#references">2</a>\].

<figure>

<img src="https://www.mdpi.com/pharmaceutics/pharmaceutics-15-02498/article_deploy/html/images/pharmaceutics-15-02498-g013.png" alt="Image Description" height="800"/>

<figcaption>
Figure 3. <a name="fig3"></a> Proposed decision tree for planning and
analysis of pilot BA/BE studies. From Henriques, S.C. et al. (2023).
Pharmaceutics, 15(10), 2498 (DOI: 10.3390/pharmaceutics15102498).
</figcaption>
</figure>

## 2. The f2PilotBE Package

The `f2PilotBE` is equipped with the following functions to aid in the
calculation of the similarity *ƒ*<sub>2</sub> factor:

- `AUC()` – calculates the cumulative area under the concentration-time
  curve (AUC) over time.

- `AUClast()` – calculates area under the concentration-time curve (AUC)
  from time zero to the time of the last measurable concentration
  (AUC<sub>0-t</sub>).

- `Cmax()` – calculates the maximum observed concentration post-dose
  (C<sub>max</sub>) directly obtained from the observed
  concentration-time profile.

- `f2.AUC()` – calculates the *ƒ*<sub>2</sub> similarity factor for AUC,
  from concentration data.

- `f2.Cmax()` – calculates the *ƒ*<sub>2</sub> similarity factor for
  C<sub>max</sub>, from concentration data.

- `f2()` – calculates the *ƒ*<sub>2</sub> similarity factor, as proposed
  by Moore and Flanner in 1996.

- `geomean()` – calculates geometric mean, i.e., the *N*th root of the
  product of the *N* observations, equivalent to `exp(mean(log(x)))`.

- `Tmax()` – calculates the time of the C<sub>max</sub>.

### 2.1. Installation

#### 2.1.1. From GitHub

To install the development version of `f2PilotBE`, start by installing
the `devtools` package from CRAN:

``` r
install.packages("devtools")
```

Then the development version of `f2PilotBE` can be installed from GitHub
as:

``` r
library(devtools)
install_github("saracarolinahenriques/f2PilotBE")
```

### 2.2. Examples

#### 2.2.1. Load Package

Following installation, load the `f2PilotBE` package and the `ggplot2`
package.

``` r
# Load packages
library("ggplot2")      # required for plotting
library("f2PilotBE")
```

#### 2.2.2. Calculate Geometric Mean Concentration

As proposed by the authors
\[<a href="#references">1</a>,<a href="#references">2</a>\], the
*ƒ*<sub>2</sub> similarity factor should be calculated from the
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

#### 2.2.3. Calculate Cmax Similarity ƒ2 Factor

The `f2.Cmax()` function allows to calculate the *ƒ*<sub>2</sub>
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

The *ƒ*<sub>2</sub> similarity factor for C<sub>max</sub> can be
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
```

``` r
tail(dta_stk)
#>     Time Trt   Conc
#> 35  3.75   T 541.50
#> 36  4.00   T 536.92
#> 37  6.00   T 440.32
#> 38  8.00   T 338.78
#> 39 12.00   T 185.03
#> 40 24.00   T  31.13
```

The *ƒ*<sub>2</sub> similarity factor for C<sub>max</sub> can be
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
    *15*(5), 1430. DOI:
    [10.3390/pharmaceutics15051430](https://doi.org/10.3390/pharmaceutics15051430).

2.  Henriques, S.C.; Paixão, P.; Almeida, L.; Silva, N.E. (2023).
    Predictive Potential of C<sub>max</sub> Bioequivalence in Pilot
    Bioavailability/Bioequivalence Studies, through the Alternative
    *ƒ*<sub>2</sub> Similarity Factor Method. *Pharmaceutics*.
    *15*(10), 2498. DOI:
    [10.3390/pharmaceutics15102498](https://doi.org/10.3390/pharmaceutics15102498).

3.  Moore, J.W.; Flanner, H.H. (1996). Mathematical Comparison of Curves
    with an Emphasis on in Vitro Dissolution Profiles. *Pharm.
    Technol.*. *20*, 64–74.
