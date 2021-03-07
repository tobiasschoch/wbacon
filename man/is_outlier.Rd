\name{is_outlier}
\alias{is_outlier.wbaconlm}
\alias{is_outlier.wbaconmv}
\title{Flag Outliers}
\usage{
is_outlier(object, ...)
\method{is_outlier}{wbaconlm}(object, ...)
\method{is_outlier}{wbaconmv}(object, ...)
}
\arguments{
	\item{object}{object of class \code{wbaconmv} or \code{wbaconlm}.}
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
