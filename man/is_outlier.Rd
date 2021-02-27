\name{is_outlier}
\alias{is_outlier.roblm}
\alias{is_outlier.robmv}
\title{Flag Outliers}
\usage{
is_outlier(object, ...)
\method{is_outlier}{roblm}(object, ...)
\method{is_outlier}{robmv}(object, ...)
}
\arguments{
	\item{object}{object of class \code{roblm} or \code{robmv}.}
	\item{...}{additional arguments passed to the method.}
}
\description{
Returns a logical vector that indicates which observations were declared
outlier by the method.
}
\value{
A logical vector.
}
\seealso{
\code{\link{wBACON_reg}} and \code{\link{wBACON}}
}
