wBACON_reg <- function(formula, weights = NULL, data, collect = 4,
	na.rm = FALSE, alpha = 0.05, version = c("V2", "V1"), maxiter = 50,
    verbose = FALSE, original = FALSE, n_threads = 2)
{
	stopifnot(alpha < 1, alpha > 0, collect > 0, maxiter > 0, collect > 0,
              n_threads > 0)
	if (!inherits(formula, "formula"))
		stop("Argument '", formula, "' must be a formula\n", call. = FALSE)

	# data preparation
	mf <- stats::model.frame(formula, data, na.action = stats::na.pass)
	mt <- stats::terms(mf)
    if (any(attr(mt, "dataClasses") == "factor"))
        stop("Factor variables are not allowed\n")

	response <- attr(mt, "response")
	y <- as.numeric(stats::model.response(mf))
	n <- length(y)
	yname <- names(mf)[response]
	x <- stats::model.matrix(mt, mf)
	if (is.null(weights))
		weights <- rep(1, n)

	# NA treatment
	cc <- stats::complete.cases(y, x, weights)
	if (sum(cc) != n) {
		if (na.rm) {
	        mf <- mf[cc, ]
			x <- x[cc, ]
			y <- y[cc]
			weights <- weights[cc]
		} else {
			stop("Data must not contain missing values; see 'na.rm'\n",
				call. = FALSE)
		}
	}
	n <- NROW(x); p <- NCOL(x)

	# check if any element is not finite
	if (sum(is.finite(c(x, y, weights))) != (2 + p) * n)
		stop("Some observations are not finite\n", call. = FALSE)

	# Algorithm 3
	if (verbose)
		cat("\nOutlier detection (Algorithm 3)\n---\n")
	wb <- wBACON(if (attr(mt, "intercept")) x[, -1] else x, weights, alpha,
                 collect, version, na.rm, maxiter, verbose)

	if (isFALSE(wb$converged))
		stop("wBACON on the design matrix failed\n")

	# Algorithms 4 and 5
	if (verbose)
		cat("\nRegression\n---\n")
	collect <- min(collect, floor(n / p))
	tmp <- .C(C_wbacon_reg, x = as.double(x), y = as.double(y),
              w = as.double(weights), resid = as.double(numeric(n)),
              beta = as.double(numeric(p)), subset = as.integer(wb$subset),
              dist = as.double(wb$dist), n = as.integer(n), p = as.integer(p),
              m = as.integer(sum(wb$subset)), verbose = as.integer(verbose),
              sucess = as.integer(1), collect = as.integer(collect),
              alpha = as.double(alpha), maxiter = as.integer(maxiter),
              original = as.integer(original),
              n_threads = as.integer(n_threads))

	# cast the QR factorization as returned by LAPACK:dgeqrf to a 'qr' object
	QR <- structure(
		list(qr = matrix(tmp$x, ncol = p), qraux = rep(NA, p), pivot = 1L:p,
             tol = NA, rank = p), class = "qr")

	# return value
	res <- list(coefficients = tmp$beta, residuals = tmp$resid, rank = p,
                fitted.values = y - tmp$resid,
                df.residual = sum(weights[tmp$subset == 1]) - p,
                call = match.call(), terms = mt, model = mf, weights = weights,
                qr = QR, subset = (tmp$subset == 1),
                reg = list(converged = as.logical(tmp$sucess),
                           collect = collect, version = version, alpha = alpha,
                           maxiter = tmp$maxiter, dist = tmp$dist,
                           cutoff = qt(alpha / (2 * (tmp$m + 1)),
                                       tmp$m - p, lower.tail = FALSE)),
                mv = list(center = wb$center, cov = wb$cov, dist = wb$dist,
                          cutoff = wb$cutoff))
	names(res$coefficients) <- colnames(x)
	class(res) <- "wbaconlm"
	res
}

print.wbaconlm <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
	if (x$reg$converged){
		n <- length(x$residuals)
		cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
            "\n", sep = "")
		n_subset <- sum(x$subset)
		cat(paste0("\nRegression on the subset of ", n_subset, " out of ", n,
                   " observations (", round(100 * n_subset / n, 1), "%)\n"))
		cat("\nCoefficients:\n")
		print.default(format(x$coefficients, digits = digits), print.gap = 2L,
                      quote = FALSE)
	} else {
		cat(paste0("Algorithm did not converge in ", x$reg$maxiter,
                   " iterations!\n\n"))
	}
	invisible(x)
}

summary.wbaconlm <- function (object, ...)
{
	in_subset <- object$subset
    r <- object$residuals[in_subset]
    w <- object$weights[in_subset]
    f <- object$fitted.values[in_subset]
    qr <- object$qr
    p <- object$rank
    rdf <- sum(w) - p
    n <- sum(w)
    mss <- if (attr(object$terms, "intercept")) {
        m <- sum(w * f / sum(w))
        sum(w * (f - m)^2)
    } else {
        sum(w * f^2)
    }
    rss <- sum(w * r^2)
    resvar <- rss / rdf
    if (is.finite(resvar) && resvar < (mean(f)^2 + stats::var(c(f))) * 1e-30)
        warning("essentially perfect fit: summary may be unreliable")
    p1 <- 1L:p
    R <- chol2inv(qr$qr[p1, p1, drop = FALSE])
    se <- sqrt(diag(R) * resvar)
    est <- object$coefficients[qr$pivot[p1]]
    tval <- est / se
    ans <- object[c("call", "terms", if (!is.null(object$weights)) "weights")]
    ans$residuals <- r
    ans$coefficients <- cbind(Estimate = est, `Std. Error` = se,
                              `t value` = tval,
                              `Pr(>|t|)` = 2 * stats::pt(abs(tval), rdf,
                                                         lower.tail = FALSE))
    ans$aliased <- is.na(object$coefficients)
    ans$sigma <- sqrt(resvar)
    ans$df <- c(p, rdf, NCOL(qr$qr))
    if (p != attr(object$terms, "intercept")) {
        df.int <- if (attr(object$terms, "intercept")) 1L else 0L
        ans$r.squared <- mss / (mss + rss)
        ans$adj.r.squared <- 1 - (1 - ans$r.squared) * ((n - df.int) / rdf)
        ans$fstatistic <- c(value = (mss / (p - df.int)) / resvar,
                            numdf = p - df.int, dendf = rdf)
    } else {
        ans$r.squared <- ans$adj.r.squared <- 0
    }
    ans$cov.unscaled <- R
    dimnames(ans$cov.unscaled) <- dimnames(ans$coefficients)[c(1, 1)]
    if (!is.null(object$na.action))
        ans$na.action <- object$na.action
    class(ans) <- c("summary.lm", "summary.wbaconlm")
    ans
}

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
