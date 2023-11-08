f2PilotBE
================

- [The Similarity f2 Factor](#the-similarity-f2-factor)
  - [For Cmax](#for-cmax)
- [Installation](#installation)
  - [From GitHub](#from-github)
- [Examples](#examples)
  - [Load Package](#load-package)
  - [Geometric Mean Concentration](#geometric-mean-concentration)
  - [Calculate f2 Factor for Cmax](#calculate-f2-factor-for-cmax)
- [References](#references)

The `f2PilotBE` package is designed to calculate the similarity factor
*f*<sub>2</sub> for C<sub>max</sub> and AUC, as an alternative analysis
approaches for the assessment of bioavailability/bioequivalence (BA/BE)
studies.

This package is based on the publications of [Henriques *et al.*
(2023)](https://doi.org/10.3390/pharmaceutics15051430). In this work,
the authors proposed the use of the geometric mean (G<sub>mean</sub>)
*f*<sub>2</sub> factor using a cut-off of 35, for the comparison of the
absorption rate (given by the maximum observed concentration
\[C<sub>max</sub>\]) of Test and Reference formulations. For the tested
simulated scenarios, G<sub>mean</sub> *f*<sub>2</sub> factor using a
cut-off of 35 showed a good relationship between avoiding type I and
type II errors.

## The Similarity f2 Factor

The similarity *f*<sub>2</sub> factor is a mathematical index widely
used to compare dissolution profiles, evaluating their similarity, using
the percentage of drug dissolved per unit of time.

The similarity *f*<sub>2</sub> factor, proposed by Moore and Flanner in
1996, is derived from the mean squared difference, and can be calculated
as a function of the reciprocal of mean squared-root transformation of
the sum of square differences at all points:

$f_{2}=50\cdot\ log \biggl(100\cdot \biggl[1+\frac{1}{n} \sum_{t=1}^{t=n}{(\bar{R}_{t} - \bar{T}_{t})^2} \biggr]^{-0.5} \biggr)$

where $f_{2}$ is the similarity factor, $n$ is the number of time
points, and $\bar{R}_{t}$ and $\bar{T}_{t}$ are the mean percentage of
drug dissolved at time $t$, for Reference and Test products
respectively.

The *f*<sub>2</sub> similarity factor ranges from 0 (when
$\bar{R}_{t} - \bar{T}_{t} = 100 \%$ at all $t$) to 100 (when
$\bar{R}_{t} - \bar{T}_{t} = 0 \%$ at all $t$).

The figure below presents the distribution of *f*<sub>2</sub> similarity
factor as a function of mean difference. *f*<sub>2</sub> similarity
factor is derived from the mean squared difference and can be calculated
as a function of the reciprocal of the mean squared-root transformation
of the sum of square differences at all points. An average difference of
<span style="color: red;">10%</span>,
<span style="color: green;">15%</span>, and
<span style="color: blue;">20%</span> from all measured time points
results in a *f*<sub>2</sub> value of 50
(<span style="color: red;">red</span> dotted lines), 41
(<span style="color: green;">green</span> dotted lines) and 35
(<span style="color: blue;">blue</span> dotted lines), respectively.

<center>
<img src="https://pub.mdpi-res.com/pharmaceutics/pharmaceutics-15-01430/article_deploy/html/images/pharmaceutics-15-01430-g002.png?1683533911">
</center>

(Figure from \[[Henriques *et al.*
(2023)](https://doi.org/10.3390/pharmaceutics15051430)\])

As previously proposed by the authors, the concept of similarity factor
*f*<sub>2</sub> can be applied as an alternative to the average
bioequivalence analysis \[[Henriques *et al.*
(2023)](https://doi.org/10.3390/pharmaceutics15051430)\]. The similarity
between Test and Reference products by means of *f*<sub>2</sub> was
evaluated through the comparison of arithmetic (A<sub>mean</sub>) and
geometric (G<sub>mean</sub>) means of plasma concentration-time profiles
derived from the simulated individual pharmacokinetic profiles.

### For Cmax

*f*<sub>2</sub> can be used to assess the similarity on the rate of drug
absorption by normalizing Test and Reference mean concentration-time
profiles to the maximum plasma concentration (C<sub>max</sub>) derived
from the mean Reference profile, until Reference C<sub>max</sub> is
observed (Reference t<sub>max</sub>) \[[Henriques *et al.*
(2023)](https://doi.org/10.3390/pharmaceutics15051430)\]:

$C_{t}^{N}=100\cdot\frac{\bar{C}_{t}}{C_{max,R}}$, where
$0 \le t \le t_{max,R}$,

where $C_{t}^{N}$ is the normalized concentration at time $t$,
$\bar{C}_{t}$ is the mean (Test or Reference) concentration at time $t$,
$C_{max,R}$ is the C<sub>max</sub> of the Reference mean
concentration-time profile, and $t_{max,R}$ the time of observation of
$C_{max,R}$. The similarity *f*<sub>2</sub> factor is calculated as

$C_{max}f_{2}=50\cdot\ log \biggl(100\cdot \biggl[1+\frac{1}{n} \sum_{t=1}^{t=n}{(R_{t}^{N} - T_{t}^{N})^2} \biggr]^{-0.5} \biggr)$,

where $n$ is the number of time points until Reference $t_{max}$, and
$R_{t}^{N}$ and $T_{t}^{N}$ are the normalized concentration at time
$t$, for Reference and Test products respectively.

## Installation

### From GitHub

To install the development version of `f2PilotBE`, start by installing
the devtools package from CRAN:

``` r
install.packages("devtools")
```

Then the development version of `f2PilotBE` can be installed from GitHub
as:

``` r
library(devtools)
install_github("LittlePetinga/f2PilotBE")
```

## Examples

### Load Package

Following installation, load the `f2PilotBE` package and the `ggplot2`
package.

``` r
# Load packages
library("ggplot2")      # required for plotting
library("f2PilotBE")
```

### Geometric Mean Concentration

Import individual concentration data and calculate the G<sub>mean</sub>
concentration-time profiles for each Treatment.

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

### Calculate f2 Factor for Cmax

Calculate the *f*<sub>2</sub> for C<sub>max</sub> from the
G<sub>mean</sub> concentration-time profiles:

``` r
# Calculate f2 Factor for Cmax
# when treatment data is pivoted:
f2.Cmax(dta, Time = 'Time', Reference = 'R', Test = 'T',
        Trt.cols = TRUE, details = FALSE, plot = TRUE)
```

In the example above, treatment data is pivoted. can also be calculated
when treatment data is stacked:

``` r
# To stack treatment information
dta <- cbind(dta[3],stack(dta[1:2]))
names(dta)[names(dta) %in% c('values','ind')] <- c('Concentration', 'Treatment')

# Calculate f2 Factor for Cmax
# when treatment data is pivoted:
f2.Cmax(dta, Time = 'Time', Concentration = 'Concentration', 
        Treatment = 'Treatment', Reference = 'R', Test = 'T',
        Trt.cols = FALSE)
```

## References

Henriques, S.C.; Albuquerque, J.; Paixão, P.; Almeida, L.; Silva, N.E.
(2023). Alternative Analysis Approaches for the Assessment of Pilot
Bioavailability/Bioequivalence Studies. *Pharmaceutics*. *15*(5), 1430.
[10.3390/pharmaceutics15051430](https://doi.org/10.3390/pharmaceutics15051430)

Henriques, S.C.; Paixão, P.; Almeida, L.; Silva, N.E. (2023). Predictive
Potential of C<sub>max</sub> Bioequivalence in Pilot
Bioavailability/Bioequivalence Studies, through the Alternative
*f*<sub>2</sub> Similarity Factor Method. *Pharmaceutics*. *15*(10),
2498.
[10.3390/pharmaceutics15102498](https://doi.org/10.3390/pharmaceutics15102498).

Moore, J.W.; Flanner, H.H. (1996). Mathematical Comparison of Curves
with an Emphasis on in Vitro Dissolution Profiles. *Pharm. Technol.*.
*20*, 64–74.
