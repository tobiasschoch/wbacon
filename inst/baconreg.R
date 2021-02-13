library(robustX)

data(iris)
x <- as.matrix(iris[, 1:3])
y <- iris[,4]
w <- rep(1,150)

#set.seed(1)
#w <- w + abs(rnorm(150, 2, 2))

tmp <- mvBACON(x)

a <- robustX:::.lmBACON(x, y, init.dis = tmp$dis)
b <- wlmBACON(x, y, w, init.dis = tmp$dis)

wlmBACON <- function (x, y, w, intercept = TRUE, init.dis, init.fraction = 0, 
    collect = 4, alpha = 0.95, maxsteps = 100, verbose = TRUE, alg5 = TRUE) 
{
    trace1 <- function(i, r, n, init.steps = FALSE, skip.init = FALSE) {
        cat("Reg-BACON (", if (init.steps) 
            "init ", "subset no. ", i, if (skip.init) 
            " after skipping init", "): ", r, " of ", 
            n, " (", round(r/n * 100, digits = 2), " %)", 
            sep = "", fill = TRUE)
    }
    ordered.indices <- order(init.dis)
    if (intercept) 
        x <- cbind(1, x)
    n <- nrow(x)
    p <- ncol(x)
    steps <- 0L
    ordered.x <- x[ordered.indices, , drop = FALSE]
    ordered.y <- y[ordered.indices]
	ordered.w <- w[ordered.indices]
    skip.init <- (init.fraction > sqrt(.Machine$double.eps))
    m <- as.integer(if (skip.init) round(init.fraction * n) else collect * p)
    if (m <= 0) {
        message(gettextf(".lmBACON(): m = %d replaced by m = 1", m), domain = NA)
        m <- 1L
    }
    subset <- is.element(1:n, ordered.indices[1:m])
    tis <- GiveTis(x, y, w, subset, ordered.x, ordered.y, ordered.w, m, TRUE)
    m <- tis$m
    tis <- abs(tis$tis)
    if (verbose) 
        trace1(steps, m, n, skip.init = skip.init, init.steps = !skip.init)
    if (skip.init) {
        r <- m
    }
    else {
		# Algorithm 4
        r <- p + 1L
        if (verbose) 
            trace1(steps, r, n, init.steps = TRUE)
        while (r < n && r < m) {
            ordered.indices <- order(tis)
            ordered.x <- x[ordered.indices, , drop = FALSE]
            ordered.y <- y[ordered.indices]
			ordered.w <- w[ordered.indices]
            subset <- is.element(1:n, ordered.indices[1:r])
            tis <- GiveTis(x, y, w, subset, ordered.x, ordered.y, ordered.w, r,
				TRUE)
			beta <- tis$beta
            r <- tis$m + 1L
            tis <- abs(tis$tis)
            steps <- steps + 1L
            if (verbose) 
                trace1(steps, r, n, init.steps = TRUE)
        }
    }

	ordered.indices <- order(tis)
	subset <- is.element(1:n, ordered.indices[1:r])

	# Algorithm 5
	if (!alg5) {
		beta <- lm.wfit(x[subset, ], y[subset], w[subset])$coefficients
	} else {
		presubset <- FALSE
		prepre.r <- pre.r <- 0L
		steps <- 0L
		while (steps <= maxsteps && !(pre.r == r && (!any(xor(presubset, 
				subset)) || prepre.r == r))) {
			presubset <- subset
			prepre.r <- pre.r
			pre.r <- r
			tmp <- GiveTis(x, y, w, subset, x[subset, , drop = FALSE],
				y[subset], w[subset])
			tis <- abs(tmp$tis)
			beta <- tmp$beta
			limit <- qt(alpha/(2 * (r + 1)), r - p, lower.tail = FALSE)
			subset <- tis < limit
			r <- sum(subset)
			steps <- steps + 1L
			if (verbose) 
				trace1(steps, r, n)
		}
	}
    list(subset = subset, tis = tis, steps = steps, beta = beta)
}

GiveTis <- function(x, y, w, subset, ordered.x, ordered.y, ordered.w, m = 1L, 
	check.rank = FALSE) {
	n <- nrow(x)
	p <- ncol(x)
	if (check.rank) {
		while (m < n && p > qr(ordered.x[1:m, , drop = FALSE])$rank) m <- m + 1L
		xm <- ordered.x[1:m, , drop = FALSE]
		ym <- ordered.y[1:m]
		wm <- ordered.w[1:m]
	}
	else {
		xm <- ordered.x
		ym <- ordered.y
		wm <- ordered.w
	}
	nm <- nrow(xm)
	stopifnot(is.matrix(xm))
	fit.m <- lm.wfit(xm, ym, wm)
	betahatm <- fit.m$coefficients
	x <- as.matrix(x)
	resid <- y - as.vector(x %*% betahatm)
	em <- ym - xm %*% betahatm 
	sigmahatm <- sqrt(sum(wm * em^2) / (sum(wm) - p))
	Rinv <- backsolve(qr.R(fit.m$qr), diag(fit.m$rank))
	Hii <- rowSums((x %*% Rinv)^2) * w
	sqroot <- 1 + (1 - 2 * subset) * Hii
	tis <- resid / (sigmahatm * sqrt(sqroot))
	list(m = m, tis = tis, beta = betahatm)
} 
