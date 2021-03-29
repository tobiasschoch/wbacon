median_w <- function(x, w, na.rm = FALSE)
{
	quantile_w(x, w, 0.5, na.rm)
}
