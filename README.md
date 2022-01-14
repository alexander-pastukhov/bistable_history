# Cumulative History Analysis For Bistable Perception Time Series

[![DOI](https://zenodo.org/badge/299245172.svg)](https://zenodo.org/badge/latestdoi/299245172)

A package to compute a cumulative history for time-series of perceptual
dominance in bistable displays.

Estimates cumulative history, an estimate of accumulating
adaptation/prediction error for the dominant percept, for time-series
for continuously viewed bistable perceptual rivalry displays. Computes
cumulative history via a homogeneous first order differential process.
I.e., it assumes exponential growth/decay of the history as a function
of time and perceptually dominant state. Supports Gamma, log normal, and
normal distribution families.

For details on rationale please refer to
([Pastukhov & Braun, 2011](https://doi.org/10.1167/11.10.12)).

## Installation

For current stable version use

```r
install.packages("bistablehistory")
```

The master branch is the development version. To install it please use

```r
library("devtools")
install_github("alexander-pastukhov/bistablehistory", dependencies = TRUE)
```

### Note

This package uses [Stan](https://mc-stan.org), a "state-of-the-art
platform for statistical modeling and high-performance statistical
computation". Therefore, it depends on the package
[rstantools](https://cran.r-project.org/package=rstantools), which in
turn depends on the [rstan](https://cran.r-project.org/package=rstan)
package, which uses the [V8 JavaScript library](https://v8.dev), through
the [V8 R package](https://cran.r-project.org/package=V8).

Therefore, you will need to install the V8 JavaScript library on your
system, and it is recommended that you also install the V8 R package
beforehand. For detailed instructions, please see
https://github.com/jeroen/v8.

You will also need the [R package
curl](https://cran.r-project.org/package=curl), which depends on
`libcurl-*` in various operating systems. Please see the documentation
at https://cran.r-project.org/package=curl.

## Usage

The main function is `fit_cumhist` that takes a data frame with
time-series as the first argument. Minimally, you need to specify `state`
--- string with the column name that encodes perceptually dominant state
--- and either `duration` (column name with duration of individual
dominance phases) or `onset` (column name with onset times of individual
dominance phases). Thus, for a simplest case of a single subject and
single run/block measurement with all defaults (gamma distribution,
fitted cumulative history time constant but fixed mixed state value and
history mixing proportion) the call would be

```r
library(bistablehistory)
data(br_singleblock)
gamma_fit <- fit_cumhist(br_singleblock,
                         state = "State",
                         duration = "Duration")
```

or, equivalently

```r
library(bistablehistory)
data(br_singleblock)
gamma_fit <- fit_cumhist(br_singleblock,
                         state = "State",
                         onset = "Time")
```

Now you can look at the fitted value for history time constant via

```r
history_tau(gamma_fit)
```

and main effect of history for both parameters of gamma distribution

```r
coef(gamma_fit)
```

For further details please see vignettes on package usage ([Usage
examples](https://CRAN.R-project.org/package=bistablehistory/vignettes/usage-examples.html)
and [Cumulative
history](https://CRAN.R-project.org/package=bistablehistory/vignettes/cumulative-history.html))
and on an example of writing Stan code directly ([Writing Stan
code](https://CRAN.R-project.org/package=bistablehistory/vignettes/writing-stan-code.html)).
