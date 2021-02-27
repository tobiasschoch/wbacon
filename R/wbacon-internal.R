# some sanity checks (univariate)
.check <- function(x, w, na.rm)
{
	if (is.factor(x) || is.factor(w) || is.data.frame(x))
		stop("Arguments data and weights must be numeric vectors\n")

	n <- length(x); nw <- length(w)
	if (nw != n)
		stop("Data vector and weights are not of the same dimension\n",
			call. = FALSE)
	if (n == 0)
		return(NA)

	# check for missing values
	cc <- stats::complete.cases(x, w)
	if (sum(cc) != n) {
		if (na.rm) {
			x <- x[cc]
			w <- w[cc]
		} else {
			return(NULL)
		}
	}
	n <- length(x)

	# check if data vector and weights are finite
	if (sum(is.finite(c(x, w))) != 2 * n) {
		warning("Some observations are not finite\n", call. = FALSE,
			immediate. = TRUE)
		return(NULL)
	}

	list(x = x, w = w, n = n)
}
