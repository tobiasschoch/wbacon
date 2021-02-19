\name{predict.roblm}
\alias{predict.roblm}
\title{Predicted values based on the weighted BACON linear regression}
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
	\item{se.fit}{A switch indicating if standard errors are required.}
	\item{scale}{Scale parameter for std.err. calculation.}
	\item{df}{Degrees of freedom for scale.}
	\item{interval}{Type of interval calculation.  Can be abbreviated.}
	\item{level}{Tolerance/confidence level.}
	\item{type}{Type of prediction (response or model term).  Can be abbreviated.}
	\item{terms}{If \code{type = "terms"}, which terms (default is all
		terms), a \code{\link{character}} vector.}
	\item{na.action}{function determining what should be done with missing
		values in \code{newdata}.  The default is to predict \code{NA}.}
	\item{pred.var}{the variance(s) for future observations to be assumed
		for prediction intervals.  See \sQuote{Details}.}
	\item{weights}{variance weights for prediction.  This can be a numeric
		vector or a one-sided model formula.  In the latter case, it is
		interpreted as an expression evaluated in \code{newdata}.}
	\item{\dots}{further arguments passed to or from other methods.}
}
\description{
This function does exactly what \code{\link[=predict.lm]{predict}} does for
the linear model \code{{lm}}. See \code{\link[=predict.lm]{predict.lm}} for
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
\code{\link{wBACON_reg}}.
}
