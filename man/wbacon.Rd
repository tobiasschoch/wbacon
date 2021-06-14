\name{wBACON}
\alias{wBACON}
\alias{print.wbaconmv}
\alias{summary.wbaconmv}
\alias{vcov.wbaconmv}
\alias{distance}
\alias{center}
\title{Weighted BACON Algorithm for Multivariate Outlier Detection}
\usage{
wBACON(x, weights = NULL, alpha = 0.05, collect = 4, version = c("V2", "V1"),
    na.rm = FALSE, maxiter = 50, verbose = FALSE, n_threads = 2)
distance(x)
\method{print}{wbaconmv}(x, digits = max(3L, getOption("digits") - 3L), ...)
\method{summary}{wbaconmv}(object, ...)
center(object)
\method{vcov}{wbaconmv}(object, ...)
}
\arguments{
	\item{x}{\code{[matrix]} or \code{[data.frame]}.}
	\item{weights}{\code{[numeric]} sampling weight (default
        \code{weights = NULL}).}
	\item{alpha}{\code{[numeric]} tuning constant, level of significance,
		\eqn{0 < \alpha < 1}{0 < alpha < 1}; (default: \code{alpha = 0.05}).}
	\item{collect}{determines the size \eqn{m} of the initial subset to be
		\eqn{m = collect \cdot p}{m = collect * p}, where \eqn{p} is the
		number of variables, \code{[integer]}.}
	\item{version}{\code{[character]} method of initialization; \code{"V1"}:
		weighted Mahalanobis distances (not robust but affine equivariant);
		\code{"V2"} (\code{default}): Euclidean norm of the data centered by
		the coordinate-wise weighted median.}
	\item{na.rm}{\code{[logical]} indicating whether \code{NA} values should
		be removed before the computation proceeds (default: \code{FALSE}).}
	\item{maxiter}{\code{[integer]} maximal number of iterations (default:
		\code{maxiter = 50}).}
	\item{verbose}{\code{[logical]} indicating whether additional information
		is printed to the console (default: \code{TRUE}).}
    \item{n_threads}{\code{[integer]} number of threads used for OpenMP
        (\code{default: 2}).}
	\item{digits}{\code{[integer]} minimal number of significant digits.}
	\item{...}{additional arguments passed to the method.}
	\item{object}{object of class \code{wbaconmv}.}
}
\description{
\code{wBACON} is an iterative method for the computation of multivariate
location and scatter (under the assumption of a Gaussian distribution).
}
\details{
The algorithm is initialized from a set of uncontaminated data. Then the
subset is iteratively refined; i.e., additional observations are included
into the subset if their Mahalanobis distance is below some threshold
(likewise, observations are removed from the subset if their distance
larger than the threshold). This process iterates until the set of good
data remain stable. Observations not among the good data are outliers;
see Billor et al. (2000). The weighted Bacon algorithm is due to Béguin
and Hulliger (2008).

The threshold for the (squared) Mahalanobis distances is defined as
the standardized chi-square \eqn{1 - \alpha}{1 - alpha} quantile. All
observations whose squared Mahalanobis distances is larger than
the threshold are regarded as outliers.

If the sampling weights \code{weights} are not explicitly specified (i.e.,
\code{weights = NULL}), they are taken to be 1.0.

\subsection{Incomplete/missing data}{
The \code{wBACON} \emph{cannot} deal with missing values. In contrast,
function \code{\link[modi]{BEM}} in package \pkg{modi} implements
the BACON-EEM algorithm of Béguin and Hulliger (2008), which
is tailored to work with outlying and missing values.

If the argument \code{na.rm} is set to \code{TRUE} the method behaves
like \code{\link{na.omit}}.
}

\subsection{Assumptions}{
The BACON algorithm \emph{assumes} that the non-outlying data have (roughly)
an elliptically contoured distribution (this includes the Gaussian
distribution as a special case). "Although the algorithms will often
do something reasonable even when these assumptions are violated, it
is hard to say what the results mean." (Billor et al., 2000, p. 289)

In line with Billor et al. (2000, p. 290), we use the term
outlier "nomination" rather than "detection" to highlight that algorithms
should not go beyond nominating observations as \emph{potential} outliers;
see also Béguin and Hulliger (2008). It is left to the analyst
to finally label outlying observations as such.
}

\subsection{Utility functions and tools}{
Diagnostic plots are available by the \code{\link[=plot.wbaconmv]{plot}}
method.

The method \code{center} and \code{vcov} return, respectively, the
estimated center/location and covariance matrix.

The \code{distance} method returns the robust Mahalanobis distances.

The function \link[=is_outlier.wbaconmv]{is_outlier} returns a vector of
logicals that flags the nominated outliers.
}
}
\value{
An object of class \code{wbaconmv} with slots

	\item{x}{see function arguments}
	\item{weights}{see function arguments}
	\item{center}{estimated center of the data}
	\item{dist}{Mahalanobis distances}
	\item{n}{number of observations}
	\item{p}{number of variables}
	\item{alpha}{see function arguments}
	\item{subset}{final subset of outlier-free data}
	\item{cutoff}{see function arguments}
	\item{maxiter}{number of iterations until convergence}
	\item{version}{see functions arguments}
	\item{collect}{see functions arguments}
	\item{cov}{covariance matrix}
	\item{converged}{logical that indicates whether the algorithm converged}
	\item{call}{the matched call}
}
\references{
Billor N., Hadi A.S., Vellemann P.F. (2000). BACON: Blocked Adaptive
	Computationally efficient Outlier Nominators. \emph{Computational
	Statistics and Data Analysis} 34, pp. 279-298.

Béguin C., Hulliger B. (2008). The BACON-EEM Algorithm for Multivariate
	Outlier Detection in Incomplete Survey Data. \emph{Survey Methodology}
	34, pp. 91-103.
}
\seealso{
\code{\link[=plot.wbaconmv]{plot}} and
\code{\link[=is_outlier.wbaconmv]{is_outlier}}
}
\examples{
data(swiss)
dt <- swiss[, c("Fertility", "Agriculture", "Examination", "Education",
    "Infant.Mortality")]
m <- wBACON(dt)
m
which(is_outlier(m))
\dontshow{stopifnot(identical(which(is_outlier(m)), 45L))}
}
