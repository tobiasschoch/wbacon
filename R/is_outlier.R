is_outlier <- function(object, ...)
{
	UseMethod("is_outlier", object)
}

is_outlier.wbaconmv <- function(object, names = FALSE, ...)
{
    out <- object$subset == 0
    if (names)
        names(which(out))
    else
        out
}

is_outlier.wbaconlm <- function(object, names = FALSE, ...)
{
    if (names)
        names(which(!object$subset))
    else
        !object$subset
}
