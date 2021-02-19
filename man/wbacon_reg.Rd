\name{wBACON_reg}
\alias{wBACON_reg}
\alias{print.roblm}
\alias{summary.roblm}
\title{Weighted BACON algorithm for robust linear regression}
\usage{
wBACON_reg(formula, weights = NULL, data, na.rm = FALSE, collect = 4,
    alpha = 0.95, version = c("V2", "V1"), maxiter = 50, verbose = FALSE) 

\method{print}{roblm}(x, digits = max(3L, getOption("digits") - 3L), ...)
\method{summary}{roblm}(object, ...)
}
\arguments{
	\item{formula}{an object of class \code{formula}: a symbolic description of
		the model to be fitted.}
	\item{weights}{\code{[numeric]} sampling weight (default 
		\code{weights = NULL}).}
	\item{data}{a \code{data.frame} object.}
	\item{na.rm}{\code{[logical]} indicating whether \code{NA} values should be
		removed before the computation proceeds (default: \code{FALSE}).}
	\item{collect}{determines the size \eqn{m} of the initial subset to be
		\eqn{m = collect \cdot p}{m = collect * p}, where \eqn{p} is the number
		of variables, \code{[integer]}.}
	\item{alpha}{\code{[double]} tuning constant, 
		\eqn{0 < \alpha < 1}{0 < alpha < 1}; (default: \code{alpha = 0.95}).}
	\item{version}{method to initialize the basic subset: Version \code{"V1"} of
		Billor et al. (2000) yields affine equivariant but \emph{not} robust
		estimators; Version \code{"V1"} yields estimators that are robust but
		not affine equivariant; (default: \code{V2}).}
	\item{maxiter}{\code{[integer]} maximal number of iterations (default:
		\code{maxiter = 50}).}
	\item{verbose}{\code{[logical]} indicating whether additional information is
		printed to the console (default: \code{TRUE}).}
	\item{digits}{\code{[integer]} minimal number of significant digits.}
	\item{object}{object of class \code{roblm}.}
	\item{...}{additional arguments passed to the method.}
}
\description{
The weighted BACON regression
}
\details{
A
}
\value{
An object of class \code{roblm} with slots

	\item{coefficients}{a named vector of coefficients}
	\item{residuals}{the residuals (for all observations in the data.fram not
		only the ones in the final subset}
	\item{rank}{the numeric rank of the fitted linear model (i.e.. number of
		variables in the design matrix}
	\item{fitted.values}{fitted values}
	\item{df.residual}{the residual degrees of freedom (computed for the
		observations in the final subset)}
	\item{call}{the matched call}
	\item{terms}{the \code{\link{terms}} object}
	\item{model}{the \code{\link{model.frame}} used}
	\item{weights}{weights}
	\item{qr}{the \code{\link{qr}} object of the linear model fit for the final
		subset}
	\item{subset}{the subset}
	\item{reg}{a list with additional details on \code{wBACON_reg}}
	\item{mv}{a list with details on the results of \code{\link{wBACON}} that
		have been used to initialize \code{wBACON_reg}}
}
\references{
Billor N., Hadi A.S., Vellemann P.F. (2000). BACON: Blocked Adaptative
	Computationally efficient Outlier Nominators. \emph{Computational
	Statistics and Data Analysis} 34, pp. 279-298.
} 
\seealso{
\code{\link[=fitted.roblm]{fitted}}, \code{\link[=residuals.roblm]{residuals}}, 
\code{\link[=coef.roblm]{coef}}, and \code{\link[=vcov.roblm]{vcov}}.
}
