\name{philips}
\alias{philips}
\docType{data}
\title{Philips data}
\description{
    The data set consists of 677 observations on 9 variables/characteristics of
    diaphragm parts for television sets.
}
\usage{data(philips)}
\format{
    A \code{data.frame} with 677 observations on the following variables:
    \describe{
        \item{\code{X1}}{\code{[double]}, characteristic 1.}
        \item{\code{X2}}{\code{[double]}, characteristic 2.}
        \item{\code{X3}}{\code{[double]}, characteristic 3.}
        \item{\code{X4}}{\code{[double]}, characteristic 4.}
        \item{\code{X5}}{\code{[double]}, characteristic 5.}
        \item{\code{X6}}{\code{[double]}, characteristic 6.}
        \item{\code{X7}}{\code{[double]}, characteristic 7.}
        \item{\code{X8}}{\code{[double]}, characteristic 8.}
        \item{\code{X9}}{\code{[double]}, characteristic 9.}
    }
}
\details{
    The data have been studied in Rousseeuw and van Driessen (1999) and
    Billor et al. (2000). They have been published in Raymaekers and Rousseeuw
    (2023).
}
\source{
Billor, N., A. S. Hadi, and P. F. Vellemann (2000).
    BACON: Blocked Adaptive Computationally-efficient Outlier Nominators.
    \emph{Computational Statistics and Data Analysis} \bold{34}, 279--298.
    \doi{10.1016/S0167-9473(99)00101-2}

Raymaekers, J. and P. Rousseeuw (2023).
    cellWise: Analyzing Data with Cellwise Outliers.
    R package version 2.5.3,
    \url{https://CRAN.R-project.org/package=cellWise}

Rousseeuw, P. J. and K. van Driessen (1999).
    A fast algorithm for the Minimum Covariance Determinant estimator.
    \emph{Technometrics} \bold{41}, 212--223.
    \doi{10.2307/1270566}
}
\examples{
head(philips)
}
\keyword{datasets}
