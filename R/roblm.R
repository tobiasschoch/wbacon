fitted.roblm <- function(object, ...)
{
	object$fitted.values
}

residuals.roblm <- function(object, ...)
{
	object$residuals
}

coef.roblm <- function(object, ...)
{
	object$coefficients
}

vcov.roblm <- function(object, ...)
{
	tmp <- summary.roblm(object, ...)
	tmp$sigma^2 * tmp$cov.unscaled
}
