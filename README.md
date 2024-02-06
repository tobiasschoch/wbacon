# wbacon: Weighted BACON algorithms for multivariate outlier nomination (detection) and robust linear regression <img src="inst/varia/logo.svg" align="right" width=120 height=139 alt="" />

[![DOI](https://joss.theoj.org/papers/10.21105/joss.03238/status.svg)](https://doi.org/10.21105/joss.03238) [![CRAN](https://www.r-pkg.org/badges/version/wbacon)](https://cran.r-project.org/package=wbacon)


## Summary

[Billor et al.](#References) (2000) proposed the BACON (blocked adaptive computationally-efficient outlier nominators) algorithms for multivariate outlier nomination and robust linear regression. [Béguin and Hulliger](#References) (2008) extended the outlier detection method to weighted and incomplete data problems.  Both methods are implemented in the R statistical software ([R Core Team](#References), 2024) in the packages, respectively, `robustX` ([Mächler et al.](#References), 2023) and `modi` ([Hulliger](#References), 2023).

Our package offers a computationally efficient implementation in the `C` language with [OpenMP](#References) support for parallelization. Efficiency is achieved by using a weighted quantile based on the Quicksort algorithm, partial sorting in place of full sorting, reuse of computed estimates, and most importantly an up-/downdating scheme for the Cholesky and QR factorizations. The computational costs of up-/downdating are far less than re-computing the entire decomposition repeatedly.

The details of the package are discussed in the accompanying paper:

> Schoch, T. (2021) wbacon: Weighted BACON algorithms for multivariate outlier nomination (detection) and robust linear regression, *Journal of Open Source Software* **6**, 323. [DOI 10.21105/joss.03238](https://doi.org/10.21105/joss.03238)



## Available methods

* `wBACON()` is for multivariate outlier nomination and robust estimation of location/ center and covariance matrix
* `wBACON_reg()` is for robust linear regression (the method is robust against outliers in the response variable and the model's design matrix)

### Assumptions

The BACON algorithms assume that the underlying model is an appropriate description of the non-outlying observations; see [Billor et al.](#References) (2000). More precisely,

* the outlier nomination method assumes that the "good" data have (roughly) an elliptically contoured distribution (this includes the Gaussian distribution as a special case);
* the regression method assumes that the non-outlying ("good") data are described by a linear (homoscedastic) regression model and that the independent variables (having removed the regression intercept/constant, if there is a constant) follow (roughly) an elliptically contoured distribution.

> "Although the algorithms will often do something reasonable even when these assumptions are violated, it is hard to say what the results mean." [Billor et al.](#References) (2000, p. 289)

It is strongly recommended that the structure of the data be examined and whether the assumptions made about the "good" observations are reasonable.

### The role of the data analyst

In line with [Billor et al.](#References) (2000, p. 290), we use the term outlier "nomination" rather than "detection" to highlight that algorithms should not go beyond nominating observations as *potential* outliers; see also [Béguin and Hulliger](#References) (2008). It is left to the analyst to finally label outlying observations as such.

The software provides the analyst with tools and measures to study potentially outlying observations. It is strongly recommended to use the tools. See the package folders `vignettes` and `doc` for a vignette (guide) and further documentation.

## Installation

The package can be installed from CRAN using
```
install.packages("wbacon")
```

## Building

Make sure that the R package `devtools` is installed. Then, the `wbacon` package can be pulled from this GitHub repository and installed by
```
devtools::install_github("tobiasschoch/wbacon")
```

The package contains `C` code that needs to be compiled. 

## Community guidelines

#### Submitting an issue

If you have any suggestions for feature additions or any problems with the software that you would like addressed with the development community, please submit an issue on the Issues tab of the project GitHub repository. You may want to search the existing issues before submitting, to avoid asking a question or requesting a feature that has already been discussed.

#### How to contribute

If you are interested in modifying the code, you may fork the project for your own use, as detailed in the GNU GPL License we have adopted for the project. In order to contribute, please contact the developer by Tobias Schoch at gmail dot com (the names are separated by a dot) after making the desired changes.

#### Asking for help

If you have questions about how to use the software, or would like to seek out collaborations related to this project, you may contact Tobias Schoch (see contact details above).

## References

Béguin, C., and Hulliger, B. (2008). The BACON-EEM Algorithm for Multivariate Outlier Detection in Incomplete Survey Data, *Survey Methodology* **34**, 91-103.

Billor, N., Hadi, A. S., and Velleman , P. F. (2000). BACON: Blocked adaptive computationally-efficient outlier nominators, *Computational Statistics and Data Analysis* **34**, 279-298.  [doi: 10.1016/S0167-9473(99)00101-2](https://doi.org/10.1016%2FS0167-9473%2899%2900101-2)

Hulliger, B. (2023). modi: Multivariate Outlier Detection and Imputation for Incomplete Survey Data, R package version 0.1-2. URL [https://CRAN.R-project.org/package=modi](https://CRAN.R-project.org/package=modi)

Mächler, M. and W. A. Stahel (2023).  robustX: ’eXtra’ / ’eXperimental’ Functionality for Robust Statistics, R package version 1.2-7. URL [https://CRAN.R-project.org/package=robustX)](https://CRAN.R-project.org/package=robustX)

OpenMP Architecture Review Board (2018). OpenMP Application Program Interface Version 5.0. URL [https://https://www.openmp.org](https://https://www.openmp.org)

R Core Team (2024). R. A language and environment for statistical computing.  R Foundation for Statistical Computing, Vienna, Austria. URL [https://www.R-project.org](https://www.R-project.org).
