is_outlier <- function(object, ...)
{
	UseMethod("is_outlier", object)
}

is_outlier.robmv <- function(object, ...)
{
	object$subset == 0
}

is_outlier.roblm <- function(object, ...)
{
	object$subset == FALSE
}
