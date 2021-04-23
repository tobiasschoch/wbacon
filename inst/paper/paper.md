---
title: 'wbacon: Weighted BACON algorithms for multivariate outlier nomination (detection) and robust linear regression'
tags:
 - R
 - outlier detection
 - robustness
 - survey
 - linear regression
 - bounded influence
authors:
 - name: Tobias Schoch
   orcid: 0000-0002-1640-3395
   affiliation: 1
affiliations:
 - name: University of Applied Sciences and Arts Northwestern Switzerland, School of Business, Riggenbachstrasse 16, CH-4600 Olten, Switzerland
   index: 1
date: 29 March 2021
bibliography: paper.bib
---

# Summary

Outlier nomination (detection) and robust regression are computationally hard problems. This is all the more true when the number of variables and observations grow rapidly. Among all candidate methods, the BACON (blocked adaptive computationally efficient outlier nominators) algorithms of @billoretal2000 have favorable computational characteristics as they require only a few model evaluations irrespective of the sample size. This makes them superior and popular algorithms (at the time of writing [Google Scholar](https://scholar.google.com) reports more than 500 citations of the @billoretal2000 paper).

`wbacon` is a package/library for the `R` statistical software [@r2021]. The library is aimed at medium to large datasets that can possibly have (sampling) weights (e.g., data from complex survey samples). The library has a user-friendly `R` interface (with plotting methods, etc.) and is written mainly in the `C` language (with OpenMP support for parallelization; see @openmp2018) for performance reasons.

# The BACON algorithms

Technically, the BACON algorithms consist of the application of series of simple statistical estimation methods such as coordinate-wise means/medians, covariance matrix, Mahalanobis distances, or least squares regression on subsets of the data. The algorithms start from an initial small subset of "good" data and keep adding those observations to the subset whose distances/ discrepancies are smaller than a predefined threshold value. The algorithms terminate if the subset cannot be increased further. The observations not in the final subset are nominated as outliers. We follow @billoretal2000 and use the term "nomination" of outliers instead of "detection" to emphasize that the algorithms should not go beyond nominating observations as potential outliers. It is left to the analyst to finally label outlying observations as such.

A naive implementation of the BACON algorithms would call the (simple) estimation methods iteratively on a sequence of growing subsets of the data without bothering too much with reusing or updating existing blocks of data. This leads to an excessively large number of copy/ modify operations and (unnecessary) recomputations. Overall, we would end up with a computationally inefficient implementation. With small datasets, the inefficiencies would probably not be noticed. With large amounts of data, however, the situation is quite different.

# Statement of need
The two BACON algorithms are available from the package `robustX` [@robustx] for the `R` statistical software. The BACON algorithm for multivariate outlier nomination has been extended  to weighted problems (in the context of survey sampling) and incomplete data by @beguinhulliger2008. The extended method is available from the `R` package `modi` [@modi]. Both implementations are not explicitly targeted at big data applications. Our library fills this gap.

# What the library offers

The implementation of the `wbacon` package is targeted at medium to large datasets and is mainly implemented in the `C` language. The code depends heavily on the subroutines in the libraries `BLAS` [@blas2002] and `LAPACK` [@lapack1999]. If computation time is of great importance to the user, we recommend replacing the reference implementation of the `BLAS` library that ships with `R` by a version that has been adapted to the user's hardware (see e.g., OpenBLAS). The non-`BLAS` components of `wbacon` use multiple threads (if the compiler supports `OpenMP`; see @openmp2018). But the major improvements of `wbacon` (in terms of computation time) over the naive implementation are achieved by using partial sorting (in place of a full sort), reusing computed estimates, and employing an up-/downdating scheme for the Cholesky and QR factorizations. The computational costs of the up-/downdating schemes are far less than recomputing the entire decomposition repeatedly.

# Diagnostic tools

The BACON algorithms assume that the underlying model is an appropriate description of the non-outlying observations. More precisely [@billoretal2000],

* the outlier nomination method assumes that the "good" data have (roughly) an elliptically contoured distribution (this includes the Gaussian distribution as a special case);
* the regression method assumes that the non-outlying ("good") data are described by a linear (homoscedastic) regression model and that the independent variables (having removed the regression intercept/constant, if there is a constant) follow (roughly) an elliptically contoured distribution.

It is strongly recommended that the structure of the data be examined and whether the assumptions made about the "good" observations are reasonable.

> "Although the algorithms will often do something reasonable even when these assumptions are violated, it is hard to say what the results mean." [@billoretal2000, p. 289]

The `wbacon` library provides the analyst with tools to study potentially outlying observations. For multivariate outlier nomination, the package implements several diagnostic plots. Worth mentioning is the graph which plots the robust (Mahalanobis) distances against the univariate projection of the data that maximizes the separation criterion of @qiujoe2006. This kind of diagnostic graph attempts to separate outlying from non-outlying observations as much as possible; see @willemsetal2009. It is particularly helpful when the outliers are clustered or show patterns. For robust linear regression, the package offers the standard plotting methods that are available for objects of the class `lm`. In addition, it implements the plot of the robust distances of the (non-constant) design variables against the standardized residuals. This diagnostic plot been proposed by @rousseeuwvanzomeren1990. All plotting methods can be displayed as hexagonally binned scatter plots, using the functionality of the `hexbin` [@hexbin] package. This option is recommended for large datasets.

# Acknowledgements
I would like to acknowledge many fruitful discussions with Beat Hulliger. This research did not receive any special grant from funding agencies in the public, commercial, or not-for-profit sectors.

# References
