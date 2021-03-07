is_outlier <- function(object, ...)
{
	UseMethod("is_outlier", object)
}

is_outlier.wbaconmv <- function(object, ...)
{
	object$subset == 0
}

is_outlier.wbaconlm <- function(object, ...)
{
	object$subset == FALSE
}
