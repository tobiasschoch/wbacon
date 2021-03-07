fitted.wbaconlm <- function(object, ...)
{
	object$fitted.values
}

residuals.wbaconlm <- function(object, ...)
{
	object$residuals
}

coef.wbaconlm <- function(object, ...)
{
	object$coefficients
}

vcov.wbaconlm <- function(object, ...)
{
	tmp <- summary.wbaconlm(object, ...)
	tmp$sigma^2 * tmp$cov.unscaled
}
