
<!-- README.md is generated from README.Rmd. Please edit that file -->

# aberrance

<!-- badges: start -->

[![](https://www.r-pkg.org/badges/version/aberrance)](https://CRAN.R-project.org/package=aberrance)
[![](https://cranlogs.r-pkg.org/badges/aberrance)](https://CRAN.R-project.org/package=aberrance)
<!-- badges: end -->

The *aberrance* package contains a collection of functions for detecting
several types of aberrant behavior, including:

- **Answer copying**, using statistics such as the $\omega$ statistic
  (Wollack, 1997).
- **Answer similarity**, using statistics such as the $GBT$ statistic
  (van der Linden & Sotaridona, 2006) and the $M4$ statistic (Maynes,
  2014).
- **Nonparametric misfit**, using statistics such as the $ZU3$ statistic
  (van der Flier, 1982) and the $H^T$ statistic (Sijtsma, 1986).
- **Parametric misfit**, using statistics such as the standardized
  log-likelihood statistic (Drasgow et al., 1985) and its various
  corrections (Bedrick, 1997; Gorney et al., 2024; Molenaar & Hoijtink,
  1990; Snijders, 2001).
- **Preknowledge**, using statistics such as the signed likelihood ratio
  test statistic (Sinharay, 2017).
- **Rapid guessing**, using methods such as the custom threshold method
  (Wise et al., 2004; Wise & Kong, 2005), the normative threshold method
  (Wise & Ma, 2012), and the cumulative proportion correct method (Guo
  et al., 2016).
- **Test tampering**, using statistics such as the erasure detection
  index (Wollack et al., 2015; Wollack & Eckerly, 2017) and its
  corrected versions (Sinharay, 2018).

## Installation

Install the released version from CRAN:

``` r
install.packages("aberrance")
```

Alternatively, install the development version from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("kyliegorney/aberrance")
```
