wBACON <- function(x, w = NULL, alpha = 0.95, version = c("V2", "V1"),
	na.rm = FALSE, maxiter = 50, verbose = FALSE) 
{ 
	n <- nrow(x); p <- ncol(x)
	stopifnot(n > p, p > 1,  0 < alpha, alpha < 1, maxiter > 0)

	if (version[1] == "V2")
		vers <- 1
	else if (version[1] == "V1")
		vers <- 0
	else
		stop(paste0("Argument '", version, "' is not defined\n"))	

	if (!is.matrix(x))
		x <- as.matrix(x)

	if (is.null(w)) 
		w <- rep(1, n)

	stopifnot(n == length(w))

	# NA treatment
	cc <- stats::complete.cases(x, w) 
	if (sum(cc) != n) {
		if (na.rm) { 
			x <- x[cc, ]
			w <- w[cc] 
		} else 
			stop("Data must not contain missing values; see argument 'na.rm'\n", 
				call. = FALSE)
	} 
	n <- nrow(x)

	# check if any element is not finite 
	chk <- sum(is.finite(c(x, w))) != (1 + p) * n 
	if (chk) 
		stop("Some observations are not finite\n", call. = FALSE)

	# compute weighted BACON algorithm
	tmp <- .C("wbacon", x = as.double(x), w = as.double(w), 
		center = as.double(numeric(p)), scatter = as.double(numeric(p * p)),
		dist = as.double(numeric(n)), n = as.integer(n), p = as.integer(p),
		alpha = as.double(alpha), subset = as.integer(rep(0, n)), 
		cutoff = as.double(numeric(1)), maxiter = as.integer(abs(maxiter)),
		verbose = as.integer(verbose), version = as.integer(vers), 
		success = as.integer(1),  
		PACKAGE = "wbacon")

	tmp$cov <- matrix(tmp$scatter, ncol = p)
	tmp$verbose <- NULL
	tmp$converged <- tmp$success == 1
	tmp$call <- match.call()
	names(tmp$center) <- colnames(x)
	colnames(tmp$cov) <- colnames(x)
	rownames(tmp$cov) <- colnames(x)
	class(tmp) <- "robmv"
	tmp
}

print.robmv <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
	cat("\nWeighted BACON: Robust location, covariance, and distances\n")
	if (x$converged) {
		cat("Initialized by method:", ifelse(x$version == 1, "V2", "V1"), "\n")
		cat(paste0("Converged in ", x$maxiter, " iterations (alpha = ", x$alpha, 
			")\n\n"))
	} else
		cat(paste0("Algorithm did not converge in ", x$maxiter, 
			" iterations!\n\n"))
}

summary.robmv <- function(object, ...)
{
	digits <- max(3L, getOption("digits") - 3L)
	cat("\nWeighted BACON: Robust location, covariance, and distances\n")
	if (object$converged) {
		cat("Initialized by method:", ifelse(object$version == 1, "V2", "V1"),
			"\n")
		cat(paste0("Converged in ", object$maxiter, " iterations (alpha = ",
			object$alpha, ")\n"))
		n <- length(object$subset)
		n_outlier <- sum(object$subset)
		cat(paste0("\nNumber of detected outliers: ", n_outlier, " (",
			round(100 * n_outlier / n, 2), "%)\n"))
		cat("\nRobust estimate of location:\n")
		print(object$center, digits = digits)
		cat("\nRobust estimate of covariance:\n")
		print(object$cov, digits = digits)
		cat(paste0("\nDistances (cutoff: ", format(object$cutoff, 
			digits = digits), "):\n"))
		print(summary(object$dist), digits = digits) 
	} else
		cat(paste0("Algorithm did not converge in ", object$maxiter, 
			" iterations!\n"))

	cat("\n")
}

distance <- function(x)
{
	if (!inherits(x, "robmv"))
		cat("not defined for this type of argument\n")
	else
		x$dist      
}

vcov.robmv <- function(object, ...)
{
	object$cov
}

center <- function(object)
{
	object$center
}


