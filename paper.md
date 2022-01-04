---
title: 'bistablehistory: an R package for history-dependent analysis of perceptual time series'
tags:
  - R
  - vision
  - time series
  - binocular riavlry
authors:
  - name: Alexander Pastukhov
    orcid: 0000-0002-8738-8591
    affiliation: "1, 2"
affiliations:
  - name: Department of General Psychology and Methodology, University of Bamberg, Bamberg, Bavaria, Germany.
    index: 1
  - name: Forschungsgruppe EPÆG (Ergonomics, Psychological Æsthetics, Gestalt), Bamberg, Bavaria, Germany.
    index: 2
date: 01 October 2021
bibliography: paper.bib
---

# Summary
Our perception is a subject to a process of adaptation that changes its operating properties [@Clifford2007]. This process manifests itself in a plethora of perceptual illusions and so-called aftereffects. It also plays an important role during perception of multistable stimuli, such as a Necker cube (Fig. 1A). These are compatible with several, typically two, comparably likely perceptual interpretations. In case of the Necker cube, one can perceive as "upwards" or "downwards" and during continuous viewing the perception switches between these alternatives (Fig. 1B).

![**Figure 1. Perception of Necker cube stimulus.** A) Necker cube. B) Schematic time series of perceptual alternations. C) Example distribution of perceptual dominance phases' duration fitting with Gamma distribution. \label{fig:fig1}](fig1.png){width=80%}

Typically, time series for such multistable stimuli are fitted using Gamma distribution (Fig. 1C). However, this approach ignores the effect prior perception through adaptation [@VanEe2009] and violates assumptions about independent and identically distributed samples. The `bistablehistory` package solves this problem by accounting for the slow accumulation process via a homogeneous first order process [@PastukhovBraun2011], providing tools for fitting time series using various distributions. It also allows for fitting while accounting for additional random or fixed factors. In addition, it provides a tool for extracting the estimated accumulated adaptation or for computing it directly for further usage. The aim of the package is to streamline the analysis of time series for perceptual multistability and experiments on continuous perception in general.

The package is built using Stan probabilistic programming language [@carpenter2017stan]. Thus, it provides posterior distributions, ability to compare models via information criteria [@loo], etc. In addition, the package provides Stan code for performing the estimation and an example that explains how to implement a custom Stan model that relies on it. The source code for `bistablehistory` has been archived to Zenodo with the linked DOI: [@zenodo]

# Statement of Need
Analysis of time series from psychophysical experiments on perceptual experiments, in particular on multistable perception, frequently requires taking into account slow accumulation of adaptation. A typical approach is to use an easy-to-compute approximation via n autocorrelation coefficient [@VanEe2009]. However, such estimates are both less accurate than a first-order process approach [@PastukhovBraun2011] and are harder to use as covariates for the time series analysis. 

Package `bistablehistory` addresses this problem by providing tools that allow to compute an estimate of this process while, optionally, fit its parameters. The estimate could be used for further analysis or directly as part of a statistical (generalized) linear model that is fitted by a package function.

# Usage and Features
The main function is `fit_cumhist()` that takes a data frame with time-series as a first argument. In addition, you need to specify the name of the column that codes the perceptual state (`state` argument) and a column that holds either dominance phase duration (`duration`) or its onset (`onset`). The code below is a simplest case scenario, fitting data using Gamma distribution (default family) for a single run of a single participant. By default, the function fits time constant for slowly accumulating adaptation but uses default values for other parameters that influence history computation.

```r
library(bistablehistory)
data(br_singleblock)
gamma_fit <- fit_cumhist(br_singleblock,
                        state="State",
                        duration="Duration",
                        refresh=0)
```

Alternatively, you specify _onset_ of individual dominance phases that will be used to compute their duration.
```r
gamma_fit <- fit_cumhist(br_singleblock,
                        state="State",
                        onset="Time")
```

Now you can look at the fitted value for history time constant using `history_tau()`
```r
history_tau(gamma_fit)
```

main effect of history for both parameters of gamma distribution
```r
historyef(gamma_fit)
```
or exract an estimate of perceptual history / adaptation via `extract_history()` function.
```r
H <- extract_history(gam_fit)
```

A package vignette provide details on including random and fixed factors, as well as on use of cumulative history computation parameters.  

# References
