predict.wbaconlm <- function(object, newdata, se.fit = FALSE, scale = NULL,
	df = Inf, interval = c("none", "confidence", "prediction"),
    level = 0.95, type = c("response", "terms"), terms = NULL,
	na.action = na.pass, ...)
{
	# on the subset, the weighted BACON regression works like a lm model
	in_subset <- object$subset == 1
	# cast 'object' to an object of class 'lm'
	ans <- object
	ans$residuals <- ans$residuals[in_subset]
	ans$fitted.values <- ans$fitted.values[in_subset]
	ans$weights <- ans$weights[in_subset]
	ans$qr$qr = ans$qr$qr[in_subset, ]
	class(ans) <- "lm"
	stats::predict.lm(ans, newdata, se.fit, scale, df, interval, level, type,
		terms, na.action, ...)
}
