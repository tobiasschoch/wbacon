\name{predict.wbaconlm}
\alias{predict.wbaconlm}
\title{Predicted Values Based on the Weighted BACON Linear Regression}
\usage{
\method{predict}{wbaconlm}(object, newdata, se.fit = FALSE, scale = NULL,
    df = Inf, interval = c("none", "confidence", "prediction"), level = 0.95,
    type = c("response", "terms"), terms = NULL, na.action = na.pass, ...)
}
\arguments{
	\item{object}{Object of class inheriting from \code{"lm"}}
	\item{newdata}{An optional data frame in which to look for variables with
		which to predict.  If omitted, the fitted values are used.}
	\item{se.fit}{A switch \code{[logical]} indicating if standard errors
		are required.}
	\item{scale}{Scale parameter for std.err. calculation, \code{[numeric]}.}
	\item{df}{Degrees of freedom for scale, \code{[integer]}.}
	\item{interval}{Type of interval calculation, \code{[character]}. Can be
		abbreviated.}
	\item{level}{Tolerance/confidence level, \code{[numeric]}.}
	\item{type}{Type of prediction (response or model term),
		\code{[character]}. Can be abbreviated.}
	\item{terms}{If \code{type = "terms"}, which terms (default is all
		terms), a \code{[character]} vector.}
	\item{na.action}{function determining what should be done with missing
		values in \code{newdata}. The default is to predict \code{NA}.}
	\item{\dots}{further arguments passed to
        \code{\link[=predict.lm]{predict.lm}}}
}
\description{
This function does exactly what \code{\link[=predict.lm]{predict}} does for
the linear model \code{lm}; see \code{\link[=predict.lm]{predict.lm}} for
more details.
}
\value{
	\code{predict.wbaconlm} produces a vector of predictions or a matrix of
	predictions and bounds with column names \code{fit}, \code{lwr}, and
	\code{upr} if \code{interval} is set.  For \code{type = "terms"} this
	is a matrix with a column per term and may have an attribute
	\code{"constant"}.

	If \code{se.fit} is
	\code{TRUE}, a list with the following components is returned:
	\item{fit}{vector or matrix as above}
	\item{se.fit}{standard error of predicted means}
	\item{residual.scale}{residual standard deviations}
	\item{df}{degrees of freedom for residual}
}
\seealso{
\code{\link{wBACON_reg}}
}
\examples{
data(iris)
m <- wBACON_reg(Sepal.Length ~ Sepal.Width + Petal.Length + Petal.Width,
    data = iris)
predict(m, newdata = data.frame(Sepal.Width = 1, Petal.Length = 1,
    Petal.Width = 1))
\dontshow{
yhat <- predict(m, newdata = data.frame(Sepal.Width = 1, Petal.Length = 1,
    Petal.Width = 1))
stopifnot(abs(yhat - 2.6594839512004644) < .Machine$double.eps)
}
}
