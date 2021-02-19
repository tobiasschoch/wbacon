library(robustX)

wlmBACON <- function (x, y, w, bacon, intercept = TRUE, 
    collect = 4, alpha = 0.95, maxsteps = 100, verbose = TRUE, alg5 = TRUE) 
{
	ordered.indices <- order(bacon$dis)
    if (intercept) 
        x <- cbind(1, x)
    n <- nrow(x)
    p <- ncol(x)
    steps <- 0L
    ordered.x <- x[ordered.indices, , drop = FALSE]
    ordered.y <- y[ordered.indices]

	ordered.w <- w[ordered.indices]					 			# TS
	subset <- bacon$subset[ordered.indices]						# TS
	m <- as.integer(sum(subset))								# TS

    tmp <- GiveTis(x, y, w, subset, ordered.x, ordered.y, ordered.w, m, TRUE)
    m <- tmp$m
    tis <- abs(tmp$tis)
    if (verbose)
        trace1(steps, m, n, skip.init = FALSE, init.steps = FALSE)
    
	# Algorithm 4
	m <- collect * p											# TS
	r <- p + 1L
	if (verbose) 
		trace1(steps, r, n, init.steps = TRUE)
	while (r < n && r < m) {
		ordered.indices <- order(tis)
		ordered.x <- x[ordered.indices, , drop = FALSE]
		ordered.y <- y[ordered.indices]
		ordered.w <- w[ordered.indices]							# TS
		subset <- is.element(1:n, ordered.indices[1:r])
		tmp <- GiveTis(x, y, w, subset, ordered.x, ordered.y, ordered.w, r,
			TRUE)
		sigma <- tmp$sigma
		beta <- tmp$beta
		r <- tmp$m + 1L
		tis <- abs(tmp$tis)
		steps <- steps + 1L
		if (verbose) 
			trace1(steps, r, n, init.steps = TRUE)
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
	fit.m <- lm.wfit(xm, ym, wm)
	betahatm <- fit.m$coefficients
	resid <- y - as.vector(x %*% betahatm)
	sigmahatm <- sqrt(sum(wm * resid[subset]^2) / (sum(wm) - p))

	Rinv <- backsolve(qr.R(fit.m$qr), diag(fit.m$rank))
	Hii <- rowSums((x %*% Rinv)^2) * w

#FIXME:
#	Hii <- diag(x %*% solve(crossprod(sqrt(wm) * xm)) %*% t(x)) * w

	sqroot <- 1 + (1 - 2 * subset) * Hii
	tis <- resid / (sigmahatm * sqrt(sqroot))
	list(m = m, tis = tis, beta = betahatm, sigma = sigmahatm)
}

trace1 <- function(i, r, n, init.steps = FALSE, skip.init = FALSE) {
	cat("Reg-BACON (", if (init.steps) 
		"init ", "subset no. ", i, if (skip.init) 
		" after skipping init", "): ", r, " of ", 
		n, " (", round(r/n * 100, digits = 2), " %)", 
		sep = "", fill = TRUE)
}

