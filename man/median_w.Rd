\name{median_w}
\alias{median_w}
\title{Weighted Median}
\usage{
median_w(x, w, na.rm = FALSE)
}
\arguments{
\item{x}{\code{[numeric vector]} observations.}

\item{w}{\code{[numeric vector]} weights (same length as vector \code{x}).}

\item{na.rm}{\code{[logical]} indicating whether \code{NA} values should be
	removed before the computation proceeds (default: \code{FALSE}).}
}
\value{
Weighted estimate of the population median.
}
\description{
\code{median_w} computes the weighted population median.
}
\details{
Weighted sample median; see \code{\link{quantile_w}} for more
information.
}
\seealso{
\code{\link{quantile_w}}
}
