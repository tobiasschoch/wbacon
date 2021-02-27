\name{predict.roblm}
\alias{predict.roblm}
\title{Predicted Values Based on the Weighted BACON Linear Regression}
\usage{
predict.roblm(object, newdata, se.fit = FALSE, scale = NULL, df = Inf,
    interval = c("none", "confidence", "prediction"), level = 0.95,
    type = c("response", "terms"), terms = NULL, na.action = na.pass,
    pred.var = res.var/weights, weights = 1, ...)
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
	\item{pred.var}{the variance(s) for future observations to be assumed
		for prediction intervals.  See \sQuote{Details}.}
	\item{weights}{variance weights for prediction \code{[numeric vector]}.
		This can be a numeric vector or a one-sided model formula. In the
		latter case, it is interpreted as an expression evaluated in
		\code{newdata}.}
	\item{\dots}{further arguments passed to or from other methods.}
}
\description{
This function does exactly what \code{\link[=predict.lm]{predict}} does for
the linear model \code{lm}; see \code{\link[=predict.lm]{predict.lm}} for
more details.
}
\value{
	\code{predict.roblm} produces a vector of predictions or a matrix of
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
