\name{wBACON_reg}
\alias{wBACON_reg}
\alias{print.wbaconlm}
\alias{summary.wbaconlm}
\title{Robust Fitting Linear Regression Models by the BACON Algorithm}
\usage{
wBACON_reg(formula, weights = NULL, data, collect = 4, na.rm = FALSE,
    alpha = 0.05, version = c("V2", "V1"), maxiter = 50, verbose = FALSE)

\method{print}{wbaconlm}(x, digits = max(3L, getOption("digits") - 3L), ...)
\method{summary}{wbaconlm}(object, ...)
}
\arguments{
	\item{formula}{an object of class \code{formula}: a symbolic description
		of the model to be fitted.}
	\item{weights}{\code{[numeric]} sampling weight (default
		\code{weights = NULL}).}
	\item{data}{a \code{data.frame} object.}
	\item{collect}{determines the size \eqn{m} of the initial subset to be
		\eqn{m = collect \cdot p}{m = collect * p}, where \eqn{p} is the
		number of variables, \code{[integer]}.}
	\item{na.rm}{\code{[logical]} indicating whether \code{NA} values should
		be removed before the computation proceeds (default: \code{FALSE}).}
	\item{alpha}{\code{[numeric]} tuning constant, level of significance,
		\eqn{0 < \alpha < 1}{0 < alpha < 1}; (default: \code{alpha = 0.05}).}
	\item{version}{method to initialize the basic subset, \code{[character]}:
		Version \code{"V1"} of Billor et al. (2000) yields affine
		equivariant but \emph{not} robust estimators; Version \code{"V1"}
		yields estimators that are robust but not affine equivariant;
		(default: \code{V2}).}
	\item{maxiter}{\code{[integer]} maximal number of iterations (default:
		\code{maxiter = 50}).}
	\item{verbose}{\code{[logical]} indicating whether additional information
		is printed to the console (default: \code{TRUE}).}
	\item{digits}{\code{[integer]} minimal number of significant digits.}
	\item{object}{object of class \code{wbaconlm}.}
	\item{x}{object of class \code{wbaconlm}.}
	\item{...}{additional arguments passed to the method.}
}
\description{
The weighted BACON algorithm is a robust method to fit weighted linear
regression models. The method is robust against outlier in the response
variable and the design matrix (leverage observation).
}
\details{
First, the \code{\link{wBACON}} method is applied to the model's design
matrix (having removed the regression intercept/constant, if there is
a constant) to establish a subset of observations which is supposed to
be free of outliers. Second, the so generated subset is regressed onto
the corresponding subset of response variables. The subset is iteratively
enlarged to include as many ``good'' observations as possible.

Models for \code{wBACON_reg} are specified symbolically. A typical model
has the form \code{response ~ terms}, where \code{response} is the
(numeric) response vector and \code{terms} is a series of terms
which specifies a linear predictor for response.

A \code{formula} has an implied intercept term. To remove this use
either \code{y ~ x - 1} or \code{y ~ 0 + x}. See \code{\link{formula}}
or \code{\link{lm}} for for more details.

The \code{weights} argument can be used to specify sampling weights or
case weights.

It is \emph{not} possible to fit multiple response variables (on the r.h.s.
of the formula, i.e. multivariate models) in one call.

The method \emph{cannot} deal with missing values. If the argument
\code{na.rm} is set to \code{TRUE} the method behaves like
\code{\link{na.omit}}.

\subsection{Assumptions}{
The algorithm \emph{assumes} that the non-outlying data follow
a \emph{linear} (homoscedastic) regression model and that the independent
variables have (roughly) an elliptically contoured distribution. ``Although
the algorithms will often do something reasonable even when these
assumptions are violated, it is hard to say what the results mean.''
(Billor et al., 2000, p. 289)

In line with Billor et al. (2000, p. 290), we use the term
outlier ``nomination'' rather than ``detection'' to highlight that algorithms
should not go beyond nominating observations as \emph{potential} outliers.
It is left to the analyst to finally label outlying observations as such.
}

\subsection{Utility functions and tools}{
The generic function \code{summary} summarizes the estimated model.

Additional utility functions are listed in the \emph{See also} section below.
}
}
\value{
An object of class \code{wbaconlm} with slots

	\item{coefficients}{a named vector of coefficients}
	\item{residuals}{the residuals (for all observations in the data.frame
		not only the ones in the final subset}
	\item{rank}{the numeric rank of the fitted linear model (i.e.. number of
		variables in the design matrix}
	\item{fitted.values}{fitted values}
	\item{df.residual}{the residual degrees of freedom (computed for the
		observations in the final subset)}
	\item{call}{the matched call}
	\item{terms}{the \code{\link{terms}} object}
	\item{model}{the \code{\link{model.frame}} used}
	\item{weights}{weights}
	\item{qr}{the \code{\link{qr}} object of the linear model fit for
		the final subset}
	\item{subset}{the subset}
	\item{reg}{a list with additional details on \code{wBACON_reg}}
	\item{mv}{a list with details on the results of \code{\link{wBACON}}
		that have been used to initialize \code{wBACON_reg}}
}
\references{
Billor N., Hadi A.S., Vellemann P.F. (2000). BACON: Blocked Adaptive
	Computationally efficient Outlier Nominators. \emph{Computational
	Statistics and Data Analysis} 34, pp. 279-298.
}
\seealso{
The generic functions \code{\link[=coef.wbaconlm]{coef}},
\code{\link[=fitted.wbaconlm]{fitted}},
\code{\link[=residuals.wbaconlm]{residuals}},
and \code{\link[=vcov.wbaconlm]{vcov}} extract the estimate coefficients,
fitted values, residuals, and the covariance matrix of the estimated
coefficients.

\code{\link[=plot.wbaconlm]{plot}} gives diagnostic plots for an
\code{wbaconlm} object.

\code{\link[=predict.wbaconlm]{predict}} is used for prediction (incl.
confidence and prediction intervals).
}
