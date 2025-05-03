\name{is_outlier}
\alias{is_outlier}
\alias{is_outlier.wbaconlm}
\alias{is_outlier.wbaconmv}
\title{Flag Outliers}
\usage{
is_outlier(object, ...)
\method{is_outlier}{wbaconmv}(object, names = FALSE, ...)
\method{is_outlier}{wbaconlm}(object, names = FALSE, ...)
}
\arguments{
	\item{object}{object of class \code{wbaconmv} or \code{wbaconlm}.}
    \item{names}{\code{[logical]}, if \code{TRUE} (default), a vector is
        returned that identifies the (potential) outliers; otherwise, the
        row names of the outliers are returned.}
	\item{...}{additional arguments passed to the method.}
}
\description{
By default the function returns a logical vector that indicates which
observations were identified or declared as (potential) outliers by the method;
if \code{names = TRUE} is set in the function call, the row names of the
(potential) outliers are returned.
}
\value{
A logical vector or vector with row names.
}
\seealso{
\code{\link{wBACON_reg}} and \code{\link{wBACON}}
}
\examples{
data(swiss)
m <- wBACON(swiss)

# indicator vector of potential outliers
is_outlier(m)

# names of the potential outliers
is_outlier(m, names = TRUE)

\dontshow{stopifnot(identical(unname(which(is_outlier(m))),
    as.integer(c(2, 3, 6:11, 31:38, 45:47))))}
}
