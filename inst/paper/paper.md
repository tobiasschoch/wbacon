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
date: "16 May 2021"
bibliography: paper.bib
---

# Summary
Outlier nomination (detection) and robust regression are computationally hard problems. This is all the more true when the number of variables and observations grow rapidly. Among all candidate methods, the two BACON (blocked adaptive computationally efficient outlier nominators) algorithms of @billoretal2000 have favorable computational characteristics as they require only a few model evaluations irrespective of the sample size. This makes them popular algorithms for multivariate outlier nomination/ detection and robust linear regression (at the time of writing [Google Scholar](https://scholar.google.com) reports more than 500 citations of the @billoretal2000 paper).

`wbacon` is a package for the `R` statistical software [@r2021]. It is aimed at medium to large datasets that can possibly have (sampling) weights (e.g., data from complex survey samples). The package has a user-friendly `R` interface (with plotting methods, etc.) and is written mainly in the `C` language (with OpenMP support for parallelization; see @openmp2018) for performance reasons.

# The BACON algorithms
Technically, the BACON algorithms consist of the application of series of simple statistical estimation methods such as coordinate-wise means/medians, covariance matrix, Mahalanobis distances, or least squares regression on subsets of the data. The algorithms start from an initial small subset of "good" data and keep adding those observations to the subset whose distances/ discrepancies are smaller than a predefined threshold value. The algorithms terminate if the subset cannot be increased further. The observations not in the final subset are nominated as outliers. We follow @billoretal2000 and use the term "nomination" of outliers instead of "detection" to emphasize that the algorithms should not go beyond nominating observations as potential outliers. It is left to the analyst to finally label outlying observations as such.

The BACON algorithm for multivariate outlier nomination can be initialized in two ways: version "V1" or "V2" [see @billoretal2000]. In version "V1", the algorithm is started from the coordinatewise median. As a consequence, the resulting estimators of location and scatter are robust (the breakdown point is approximately 40\%, see @billoretal2000) but not affine equivariant estimators of the population location and scatter. However, @billoretal2000 show that the estimators are *nearly* affine equivariant. The initialization by version "V1" yields estimators that are affine equivariant by design because the algorithm is started from the coordinatewise mean but the estimators have a very low breakdown point.

A naive implementation of the BACON algorithms would call the (simple) estimation methods iteratively on a sequence of growing subsets of the data without bothering too much with reusing or updating existing blocks of data. This leads to an excessively large number of copy/ modify operations and (unnecessary) recomputations. Overall, we would end up with a computationally inefficient implementation. With small datasets, the inefficiencies would probably not be noticed. With large amounts of data, however, the situation is quite different.

# Statement of need
The two BACON algorithms are available from the package `robustX` [@robustx] for the `R` statistical software. The BACON algorithm for multivariate outlier nomination has been extended  to weighted problems (in the context of survey sampling) and incomplete data by @beguinhulliger2008. The extended method is available from the `R` package `modi` [@modi]. Both implementations are not explicitly targeted at large data sets. Our package fills this gap.

# What the package offers
The implementation of the `wbacon` package is targeted at medium to large datasets and is mainly implemented in the `C` language. The code depends heavily on the subroutines in the libraries `BLAS` [@blas2002] and `LAPACK` [@lapack1999]. If computation time is of great importance to the user, we recommend replacing the reference implementation of the `BLAS` library that ships with `R` by a version that has been adapted to the user's hardware (see e.g., OpenBLAS). The non-`BLAS` components of `wbacon` use multiple threads (if the compiler supports `OpenMP`; see @openmp2018) for the computations over the $p$ variables/ columns because the computational time complexity is dominated by $p$. For instance, the complexity of the BACON algorithm for multivariate outlier nomination is determined by the complexity of the covariance matrix computation, which is of order $\mathcal{O}(np^2)$, where $n$ denotes the number of observations. The major improvements of `wbacon` (in terms of computation time) over the naive implementation, however, are achieved by using partial sorting (in place of a full sort), reusing computed estimates, and employing an up-/downdating scheme for the Cholesky and QR factorizations. The computational costs of the up-/downdating schemes are far less than recomputing the entire decomposition repeatedly.

# Diagnostic tools
The BACON algorithms assume that the underlying model is an appropriate description of the non-outlying observations. More precisely [@billoretal2000],

* the outlier nomination method assumes that the "good" data have (roughly) an elliptically contoured distribution (this includes the Gaussian distribution as a special case);
* the regression method assumes that the non-outlying ("good") data are described by a linear (homoscedastic) regression model and that the independent variables (having removed the regression intercept/constant, if there is a constant) follow (roughly) an elliptically contoured distribution.

It is strongly recommended that the structure of the data be examined and whether the assumptions made about the "good" observations are reasonable.

> "Although the algorithms will often do something reasonable even when these assumptions are violated, it is hard to say what the results mean." [@billoretal2000, p. 289]

The `wbacon` library provides the analyst with tools to study potentially outlying observations. For multivariate outlier nomination, the package implements several diagnostic plots. Worth mentioning is the graph which plots the robust (Mahalanobis) distances against the univariate projection of the data that maximizes the separation criterion of @qiujoe2006. This kind of diagnostic graph attempts to separate outlying from non-outlying observations as much as possible; see @willemsetal2009. It is particularly helpful when the outliers are clustered or show patterns. For robust linear regression, the package offers the standard plotting methods that are available for objects of the class `lm`. In addition, it implements the plot of the robust distances of the (non-constant) design variables against the standardized residuals. This diagnostic plot been proposed by @rousseeuwvanzomeren1990. All plotting methods can be displayed as hexagonally binned scatter plots, using the functionality of the `hexbin` [@hexbin] package. This option is recommended for large datasets.

# Illustration

In this section, we study the BACON algorithm for robust linear regression. Our data are on education expenditures in 50 US states in 1975 [@chatterjeehadi2006, Chap. 5.7]. The data can be loaded from the `robustbase` package [@robustbase] package.

```{.r}
library(wbacon)
data(education, package = "robustbase")

names(education)[3:6] <- c("RES", "INC", "YOUNG", "EXP")
head(education)
```

The variables are:

 * `State`: State
 * `Region`: group variable with outcomes: 1=Northeastern, 2=North central, 3=Southern, and 4=Western
 * `RES`: Number of residents per thousand residing in urban areas in 1970
 * `INC`: Per capita personal income in 1973 (\$US)
 * `YOUNG`: Number of residents per thousand under 18 years of age in 1974
 * `EXP`: Per capita expenditure on public education in a state (\$US), projected for 1975

Our goal is to regress education expenditures (`EXP`) on the variables `RES`, `INC`, and `YOUNG`. For the BACON robust linear regression algorithm, we have

```{.r}
reg <- wBACON_reg(EXP ~ RES + INC + YOUNG, data = education)

reg

#> Call:
#> wBACON_reg(formula = EXP ~ RES + INC + YOUNG, data = education)
#>
#> Regression on the subset of 49 out of 50 observations (98%)
#>
#> Coefficients:
#> (Intercept)          RES          INC        YOUNG
#>  -277.57731      0.06679      0.04829      0.88693
```

By default, `wBACON_reg()` uses the parametrization `alpha = 0.05`, `collect = 4`, and `version = "V2"`. These parameters are used to call the `wBACON()` multivariate outlier nomination/ detection algorithm on the design matrix. Then, the same parameters are used to compute the robust linear regression.

To ensure a high breakdown point, `version = "V2"` should not be changed to "V1" unless you have good reasons to do so. The main "turning knob" to tune the algorithm is `alpha`, which defines the $(1 - \alpha)$ quantile of the Student $t$-distribution. All observations whose distances/discrepancies are smaller (in absolute value) than the quantile (times a constant) are selected into the subset of "good" data [see document `methods.pdf` in the folder `doc` of the package]. By choosing larger values for `alpha` (e.g., 0.2), more observations are selected (ceteris paribus) into the subset of "good" data (and vice versa).

The parameter `collect` specifies the initial subset size, which is defined as $m = p \cdot collect$. It should be chosen such that $m$ is considerably smaller than the number of observations $n$. Otherwise we are at risk of selecting too many "bad" observations into the initial subset, which will eventually bias the regression estimates.

The instance `reg` is an object of the class `wbaconlm`. The printed output of `wBACON_reg()` is identical with the one of the `stats::lm()` function. In addition, we are told the size of the subset on which the regression has been computed. The observations not in the subset are considered potential outliers (here 1 out of 50 observations). The `summary()` method can be used to summarize the estimated model.

```{.r}
summary(reg)

#> Call:
#> wBACON_reg(formula = EXP ~ RES + INC + YOUNG, data = education)
#>
#> Residuals:
#>     Min      1Q  Median      3Q     Max
#> -81.128 -22.154  -7.542  22.542  80.890
#>
#> Coefficients:
#>               Estimate Std. Error t value Pr(>|t|)
#> (Intercept) -277.57731  132.42286  -2.096 0.041724 *
#> RES            0.06679    0.04934   1.354 0.182591
#> INC            0.04829    0.01215   3.976 0.000252 ***
#> YOUNG          0.88693    0.33114   2.678 0.010291 *
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#>
#> Residual standard error: 35.81 on 45 degrees of freedom
#> Multiple R-squared:  0.4967,    Adjusted R-squared:  0.4631
#> F-statistic:  14.8 on 3 and 45 DF,  p-value: 7.653e-07
```

The methods `coef()`, `vcov()`, and `predict()` work exactly the same as their `lm()` counterparts. This is also true for the first three `plot()` types, that is

 * `which = 1`: Residuals vs Fitted,
 * `which = 2`: Normal Q-Q,
 * `which = 3`: Scale-Location

The plot types `4:6` of `plot.lm()` are not implemented for objects of the class `wbaconlm` because it is not sensible to study the standard regression influence diagnostics in the presence of outliers in the model's design space (leverage observations). Instead, type four (`which = 4`) plots the robust Mahalanobis distances with respect to the non-constant design variables against the standardized residual. This plot has been proposed by [Rousseeuw and van Zomeren](#biblio) (1990).

See vignette to learn more about the package.

# Benchmarking

We benchmark our implementation against `robustX::BACON()`. Our implementation is intentionally *limited to single-threaded* computations (`n_threads = 1`; no OpenMP parallelization support). It is clear that when the number of variables is large, the parallelized computations are (usually) much faster.

We consider estimating a robust linear regression for a Gaussian mixture distribution, where a proportion of $1- \epsilon$ of the observations on the $p$ independent variables is generated by the Gaussian model, while a proportion of $\epsilon$ (the outliers) is generated by a shifted Gaussian distribution. For the outlying observations (i.e., $\epsilon$ proportion of the data), the response variable is generated by a regression coefficient which is 10 times larger than the coefficient of the non-outlying observations. We chose $\epsilon = 0.05$; see Appendix for more details on the setup. The number of variables ($p$) and the number of observations ($n$) are varied. Table 1 shows the ratio of average compute time of the two implementation for some configurations of $(n,p)$. A ratio $> 1.0$ ($< 1.0$) implies that `wBACON_reg()` is faster (slower) than `robustX:BACON()`. The average ratio refers to compute time averages over repeated benchmarks using the R package `microbenchmark` [@microbenchmark]. It is evident from the results in Table 1 that `wBACON_reg()` is considerably faster than its competitor, e.g., `wBACON_reg()` is on average 4.4 times faster for the setup $p=5$ and $n=100$. More importantly, the differences become larger as we increase $n$ or $p$. The differences in compute speed are mainly due to the fact that `wBACON_reg()` updates the regression estimates as the subset of non-outlying observations grows, while `robustX:BACON()` recomputes the estimates at each iteration.

| No. of variables $p$ | No. of observations $n$ | Ratio |
| :------------------- | :---------------------- | ----: |
|  5 |  100 |  4.4 |
| 10 |  100 |  5.3 |
| 20 |  100 |  7.1 |
|  5 | 1000 | 43.9 |
| 10 | 1000 | 52.5 |
| 20 | 1000 | 55.0 |

Table 1: Benchmarks: Robust linear regression


| No. of variables | No. of observations | Ratio |
| :--------------- | :------------------ | ----: |
|  5 |    1,000 | 5.9 |
| 10 |    1,000 | 3.9 |
| 20 |    1,000 | 2.5 |
|  5 |   10,000 | 6.5 |
| 10 |   10,000 | 3.7 |
| 20 |   10,000 | 2.8 |
|  5 |  100,000 | 7.3 |
| 10 |  100,000 | 4.2 |
| 20 |  100,000 | 3.3 |
|  5 | 1,000,000 | 6.2 |
| 10 | 1,000,000 | 4.1 |
| 20 | 1,000,000 | 3.2 |

Table 2: Benchmarks: Multivariate outlier nomination/ detection



| No. of variables | No. of observations | Ratio |
| :--------------- | :------------------ | ----: |
|  5 |    1,000 | 5.9 |
| 10 |    1,000 | 3.9 |
| 20 |    1,000 | 2.5 |

|  5 |   10,000 | 6.5 |
| 10 |   10,000 | 3.7 |
| 20 |   10,000 | 2.8 |

|  5 |  100,000 | 7.3 |
| 10 |  100,000 | 4.2 |
| 20 |  100,000 | 3.3 |

|  5 | 1,000,000 | 6.2 |
| 10 | 1,000,000 | 4.1 |
| 20 | 1,000,000 | 3.2 |


| n         |  p = 5 | p = 10 | p = 20 |
| :-------- | -----: | -----: | -----: |
|     1,000 | 5.9 | 3.9 | 2.5 |
|    10,000 | 6.5 | 3.7 | 2.8 |
|   100,000 | 7.3 | 4.2 | 3.3 |
| 1,000,000 | 6.2 | 5.1 | 3.2 |


# Community guidelines

## Submitting an issue

If you have any suggestions for feature additions or any problems with the software that you would like addressed with the development community, please
submit an issue on the Issues tab of the project GitHub repository. You may want to search the existing issues before submitting, to avoid asking a
question or requesting a feature that has already been discussed.

## How to contribute

If you are interested in modifying the code, you may fork the project for your own use, as detailed in the GNU GPL License we have adopted for the
project. In order to contribute, please contact Tobias Schoch after making the desired changes.

# Acknowledgements

I would like to acknowledge many fruitful discussions with Beat Hulliger. This research did not receive any special grant from funding agencies in the public, commercial, or not-for-profit sectors.

# Appendix

Consider the Gaussian mixture distribution $G = (1 - \epsilon) \cdot N(0 \cdot 1_p, I_p) + \epsilon \cdot N(4 \cdot 1_p, I_p)$, where $\epsilon = 0.05$ (amount of contamination), $N$ is the cumulative distribution function of the $p$-variate standard Gaussian distribution, $I_p$ and $1_p$ are, respectively, the $(p \times p)$ identity matrix and the $p$-vector of ones. We generate the $(\lfloor \epsilon n\rfloor \times p)$ design matrix $X_{good}$ of "good" observations from the $N(0 \cdot 1_p, I_p)$ distribution; $n$ denotes the sample size. The design matrix $X_{bad}$ consisting of $\lceil (1 - \epsilon) n \rceil$ "bad" observations is generated from the $N(4 \cdot 1_p, I_p)$ distribution. 


For the regression analysis, we generate the vectors of the response variable $y_{good} = X_{good} 1_p + e$ and $y_{bad} = X_{bad} (10 \cdot 1_p) + e$, where $e$ is a random error with standard Gaussian distribution. In the simulation, $y$ is regressed on $X$, where $y = (y_{good}^T, y_{bad}^T)^T$ and $X = (X_{good}^T, X_{bad}^T)^T$.

Compute environment: R version 4.0.3 (2020-10-10), x86_64-w64-mingw32, Windows 10 x64 (build 19041), Intel Core i7-8550U CPU (8th generation mobile processor, released in 2017), 1.80 GHz base clock.

# References
