\name{roblm}
\alias{fitted.roblm}
\alias{residuals.roblm}
\alias{coef.roblm}
\alias{vcov.roblm}
\title{Accessing Weighted BACON Linear Model Fit}
\usage{
\method{fitted}{roblm}(object, ...)
\method{residuals}{roblm}(object, ...)
\method{coef}{roblm}(object, ...)
\method{vcov}{roblm}(object, ...)
}
\arguments{
	\item{object}{object of class \code{roblm}.}
	\item{x}{object of class \code{roblm}.}
	\item{...}{additional arguments passed to the method.}
}
\description{
The generic accessor functions \code{coef}, \code{fitted}, \code{residuals},
and \code{vcov} can be used to extract various useful features of the value
returned by \code{\link{wBACON_reg}}.
}
\seealso{
\code{\link{wBACON_reg}}, \code{\link[=summary.roblm]{summary}},
\code{\link[=predict.roblm]{predict}}, and \code{\link[=plot.roblm]{plot}}
}
