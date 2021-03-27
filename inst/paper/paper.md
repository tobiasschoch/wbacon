---
title: 'wbacon: Weighted BACON algorithms for multivariate outlier nomination (detection) and robust linear regression'
tags:
 - R
 - outlier detection
 - robustness
 - survey
authors:
 - name: Tobias Schoch
   orcid: 0000-0002-1640-3395
   affiliation: 1
affiliations:
 - name: University of Applied Sciences and Arts Northwestern Switzerland, School of Business
index: 1
date: 26 March 2021
bibliography: paper.bib
---

# Summary

A summary describing the high-level functionality and purpose of the software for a diverse, non-specialist audience.

# The BACON algorithms

Outlier nomination (detection) and robust regression are computationally hard problems. This is all the more true when the number of variables and observations grow rapidly. Among all candidate methods, the BACON (blocked adaptive computationally efficient outlier nominators) algorithm of @billoretal2000 has favorable computational characteristics as it requires only a few model evaluation irrespective of the sample size. This makes it a superior algorithm for big data applications. It is also a quite popular method (google scholar reports more than 500 citations in February 2021). The BACON algorithms for multivariate outliers nomination and robust linear regression are available  from the R package `robustX` [@robustx] for the R statistical software [@r2021].

The BACON algorithm for multivariate outlier nomination has been by @beguinhulliger2008 extended to weighted problems (in the context of survey sampling) and incomplete data. The extended method is available from the R package `modi` [@modi].

# Statement of need

The existing implementations (packages `robustx` and `modi`) are not explicitly targeted at big data applications. There is a clear need for this.

# What the package offers

Technically, the BACON algorithms consist of the application of series of simple statistical estimation methods such as coordinate-wise means, covariance matrix, Mahalanobis distances, or least squares regression on subsets of the data. A naive implementation would call the estimation methods iteratively on a sequence of growing subsets of the data without bothering too much with reusing or updating existing blocks of data. This leads to an excessively large number of copy/ modify operations and (unnecessary) recomputations. Overall, we would end up with a computationally inefficient implementation.

The implementation of the wbacon package is targeted at medium to large datasets. In order to ensure reasonable compute time for such data, the code is mainly implemented in the C language with an API for the R statistical software. BLAS  [@blas2002] and LAPACK [@lapack1999] and parallelization by OpenMP [@openmp2018]. OpenBLAS

we use a partial sorting device (based on Quicksort), reuse of computed estimates, and employ an up-/downdating scheme for the Cholesky and QR factorizations. The computational costs of the up-/downdating schemes are far less than recomputing the entire decomposition repeatedly.

# Outlier: no automatic

[@willemsetal2009]

[@rousseeuwvanzomeren1990]

[@qiujoe2006]

hexbin [@hexbin]

# Acknowledgements

# References
