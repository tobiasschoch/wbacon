wBACON <- function(x, weights = NULL, alpha = 0.05, collect = 4,
	version = c("V2", "V1"), na.rm = FALSE, maxiter = 50, verbose = FALSE,
    n_threads = 2)
{
	n <- NROW(x); p <- NCOL(x)
	stopifnot(n > p, p > 0, 0 < alpha, alpha < 1, maxiter > 0, collect > 1,
              n_threads > 0)

	if (version[1] == "V2")
		vers <- 1
	else if (version[1] == "V1")
		vers <- 0
	else
		stop(paste0("Argument '", version, "' is not defined\n"))

	if (!is.matrix(x))
		x <- as.matrix(x)

	if (is.null(weights))
		weights <- rep(1, n)

	stopifnot(n == length(weights))

	# NA treatment
	cc <- stats::complete.cases(x, weights)
	if (sum(cc) != n) {
		if (na.rm) {
			x <- x[cc, ]
			weights <- weights[cc]
		} else
			stop("Data must not contain missing values; see argument 'na.rm'\n",
				 call. = FALSE)
	}
	n <- nrow(x)

	# check if any element is not finite
	chk <- sum(is.finite(c(x, weights))) != (1 + p) * n
	if (chk)
		stop("Some observations are not finite\n", call. = FALSE)

	# check if collect is corretly specified
	if (collect >= n / p)
		stop("Argument 'collect' must be an integer smaller than ",
             floor(n / p), "\n")
	if (collect * p / n > 0.6 && verbose)
		cat("Note: initial subset > 60% (use a smaller value for 'collect')\n")

	# compute weighted BACON algorithm
	tmp <- .C(C_wbacon, x = as.double(x), w = as.double(weights),
              center = as.double(numeric(p)),
              scatter = as.double(numeric(p * p)),
              dist = as.double(numeric(n)), n = as.integer(n),
              p = as.integer(p), alpha = as.double(alpha),
              subset = as.integer(rep(0, n)), cutoff = as.double(numeric(1)),
              maxiter = as.integer(abs(maxiter)),
              verbose = as.integer(verbose), version = as.integer(vers),
              collect = as.integer(collect), success = as.integer(1),
              n_threads = as.integer(n_threads))

    tmp$cutoff <- sqrt(tmp$cutoff)
 	tmp$verbose <- NULL
	tmp$converged <- tmp$success == 1
	tmp$success <- NULL
    tmp$x <- matrix(tmp$x, ncol = p)

    if (!tmp$converged) {
        tmp$center <- rep(NA, p)
        tmp$cov <- matrix(rep(NA, p * p), ncol = p)
        tmp$dist <- rep(NA, n)
        tmp$subset <- rep(NA, n)
        tmp$cutoff <- NA
    } else {
        tmp$cov <- matrix(tmp$scatter, ncol = p)
	    tmp$cov <- tmp$cov + t(tmp$cov * lower.tri(tmp$cov))
    }
	names(tmp$center) <- colnames(x)
	colnames(tmp$cov) <- colnames(x)
	rownames(tmp$cov) <- colnames(x)
    tmp$scatter <- NULL

	tmp$call <- match.call()
	class(tmp) <- "wbaconmv"
	tmp
}

print.wbaconmv <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
	if (x$converged) {
        cat("\nWeighted BACON: Robust location, covariance, and distances\n")
		cat(paste0("Converged in ", x$maxiter, " iterations (alpha = ",
                   x$alpha, ")\n"))
        n_outlier <- x$n - sum(x$subset)
        cat(paste0("Number of potential outliers: ", n_outlier, " (",
                   round(100 * n_outlier / x$n, 2), "%)\n\n"))
	} else
		cat(paste0("Weighted BACON did not converge in ", x$maxiter,
                   " iterations!\n\n"))
}

summary.wbaconmv <- function(object, ...)
{
	digits <- max(3L, getOption("digits") - 3L)
	cat("\nWeighted BACON: Robust location, covariance, and distances\n")
    cat("Initialized by method:", ifelse(object$version == 1, "V2", "V1"), "\n")
    if (object$converged)
        cat(paste0("Converged in ", object$maxiter, " iterations (alpha = ",
                   object$alpha, ")\n"))
    else
        cat(paste0("\nDID NOT CONVERGE in ", object$maxiter,
                   " iterations (alpha = ", object$alpha, ")\n"))
    n <- length(object$subset)
    n_outlier <- object$n - sum(object$subset)
    cat(paste0("\nNumber of potential outliers: ", n_outlier, " (",
               round(100 * n_outlier / n, 2), "%)\n"))
    cat("\nRobust estimate of location:\n")
    print(object$center, digits = digits)
    cat("\nRobust estimate of covariance:\n")
    print(object$cov, digits = digits)
    cat(paste0("\nDistances (cutoff: ",
               format(object$cutoff, digits = digits), "):\n"))
    if (object$converged)
        print(summary(object$dist), digits = digits)
    else
        print(NA)
	cat("\n")
}

distance <- function(x)
{
	if (!inherits(x, "wbaconmv"))
		cat("not defined for this type of argument\n")
	else
		x$dist
}

vcov.wbaconmv<- function(object, ...)
{
	object$cov
}

center <- function(object)
{
	object$center
}
