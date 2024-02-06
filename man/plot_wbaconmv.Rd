\name{plot.wbaconmv}
\alias{plot.wbaconmv}
\alias{SeparationIndex}
\title{Plot Diagnostics for an Object of Class \code{wbaconmv}}
\usage{
\method{plot}{wbaconmv}(x, which = 1:2,
    caption = c("Robust distance vs. Index",
    "Robust distance vs. Univariate projection"), hex = FALSE, col = 2,
    pch = 19, ask = prod(par("mfcol")) < length(which) && dev.interactive(),
    alpha = 0.05, maxiter = 20, tol = 1e-5, ...)
SeparationIndex(object, alpha = 0.05, tol = 1e-5, maxiter = 20)
}
\arguments{
	\item{x}{object of class \code{wbaconmv}}
	\item{which}{if a subset of the plots is required, specify a subset of
		the numbers \code{1:2}, \code{[integer]}.}
	\item{caption}{captions to appear above the plots;
		\code{[character]} vector of valid graphics annotations. It
		can be set to \code{""} or \code{NA} to suppress all captions.}
	\item{hex}{toogle the hexagonal bin plot on/off \code{[logical]} (default:
        \code{hex = FALSE})}
    \item{col}{color of outliers, \code{[integer]} (default: \code{col = 2})}
    \item{pch}{plot character of outliers, \code{[integer]} (default:
        \code{pch = 19})}
	\item{ask}{\code{[logical]}; if \code{TRUE}, the user is \emph{ask}ed before
		each plot, see \code{\link{par}(ask=.)}.}
	\item{alpha}{\code{[numeric]} tuning constant, level of significance,
		\eqn{0 < \alpha < 1}{0 < alpha < 1}; (default: \code{alpha = 0.05}).}
	\item{maxiter}{\code{[integer]} maximal number of iterations (default:
		\code{maxiter = 20}).}
    \item{tol}{numerical termination criterion, \code{[numeric]} (default:
        \code{tol = 1e-5})}
	\item{object}{object of class \code{wbaconmv}}
	\item{...}{additional arguments passed to the method.}
}
\description{
Two plots (selectable by \code{which}) are available for an object of class
\code{wbaconmv}: (1) Robust distance vs. Index and (2) Robust distance
vs. Univariate projection.
}
\details{
The first plot (\code{which = 1}) is a standard diagnostic tool which plots
the observations' index (\code{1:n}) against.the robust (Mahalanobis)
distances; see. e.g., Rousseeuw and van Driessen (1999).

The second plot (\code{which = 2}) plots the univariate projection of
the data which maximizes the separation criterion for clusters of
Qui and Joe (2006) against.the robust (Mahalanobis) distances. This plot
is due to Willems et al. (2009).

For large data sets, it is recommended to specify the argument
\code{hex = TRUE}. This option shows a hexagonally binned scatterplot
in place of the classical scatterplot.
}
\value{
[no return value]
}
\references{
Rousseeuw, P.J. and K. van Driessen (1999).
    A Fast Algorithm for the Minimum Covariance Determinant,
    \emph{Technometrics} \bold{41},  212--223. \doi{10.2307/1270566}

Qiu, W. and H. Joe (2006).
    Separation index and partial membership for clustering,
    \emph{Computational Statistics and Data Analysis} \bold{50}, 585--603.
    \doi{10.1016/j.csda.2004.09.009}

Willems, G., H. Joe, and R. Zamar (2009).
    Diagnosing Multivariate Outliers Detected by Robust Estimators,
    \emph{Journal of Computational and Graphical Statistics} \bold{18}, 73--91.
    \doi{10.1198/jcgs.2009.0005}
}
\seealso{
\code{\link{wBACON}}
}
