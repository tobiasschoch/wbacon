\name{wbaconlm}
\alias{fitted.wbaconlm}
\alias{residuals.wbaconlm}
\alias{coef.wbaconlm}
\alias{vcov.wbaconlm}
\title{Accessing Weighted BACON Linear Model Fit}
\usage{
\method{fitted}{wbaconlm}(object, ...)
\method{residuals}{wbaconlm}(object, ...)
\method{coef}{wbaconlm}(object, ...)
\method{vcov}{wbaconlm}(object, ...)
}
\arguments{
	\item{object}{object of class \code{wbaconlm}.}
	\item{...}{additional arguments passed to the method.}
}
\description{
The generic accessor functions \code{coef}, \code{fitted}, \code{residuals},
and \code{vcov} can be used to extract various useful features of the value
returned by \code{\link{wBACON_reg}}.
}
\seealso{
\code{\link{wBACON_reg}}, \code{\link[=summary.wbaconlm]{summary}},
\code{\link[=predict.wbaconlm]{predict}}, and \code{\link[=plot.wbaconlm]{plot}}
}
