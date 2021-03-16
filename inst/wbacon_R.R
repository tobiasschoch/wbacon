bc <- function(x, init = "V2", collect = 4, m = min(collect * p, n * 0.5), 
		alpha = 0.95, maxsteps = 100, verbose = TRUE){
	x <- as.matrix(x); n <- nrow(x); p <- ncol(x)

	# initialize the subsets 
	subset0 <- rep(FALSE, n)
	subset1 <- switch(init,
		"dUniMedian" = init_dUniMedian(x, m, verbose),
		"V1" = init_v1(x, m, verbose),
		"V2" = init_v2(x, m, verbose),
		"median" = init_median(x, m, verbose))
	if (verbose)
		cat("subset:", sum(subset1), "of", n, "\n")

	# iterate
	converged <- FALSE; steps <- 1L
	repeat {
		dis <- mahalanobis_distance(x, subset1)
#FIXME:
print(which(subset1 == TRUE))		
		subset0 <- subset1
		subset1 <- dis < chi2_cutoff(alpha, n, p, sum(subset1))
		if (verbose)
			cat("subset:", sum(subset1), "of", n, "\n")
		converged <- !any(xor(subset0, subset1))
		if (converged | steps > maxsteps)
			break
		steps <- steps + 1L
	}

	list(dis = dis, subset = subset1, center = colMeans(x[subset1, ]), 
		scatter = cov(x[subset1, ]), steps = steps, converged = converged)
} 

init_median <- function(x, m, verbose)
{
	n <- NROW(x)
	center <- apply(x, 2L, quantile, probs = 0.5, type = 1)
	x_centered <- sweep(x, 2L, center)
	scatter <- crossprod(x_centered) / (n - 1)
	dist <- mahalanobis(x_centered, 0, scatter)
	# sort by dist
	order_index <- order(dist)
	x_ordered <- x[order_index, , drop = FALSE]
	# check for rank deficiency of scatter matrix
	m <- check_rank(x_ordered, n, m, verbose)
	1:n %in% order_index[1:m]
}

init_dUniMedian <- function(x, m, verbose)
{
	n <- NROW(x)
	x_centered <- sweep(x, 2, colMedians(x))
	# scatter: covariance of the centered observations!
	scatter <- cov(x_centered) 
	dist <- mahalanobis(x_centered, 0, scatter)
	# sort by dist
	order_index <- order(dist)
	x_ordered <- x[order_index, , drop = FALSE]
	# check for rank deficiency of scatter matrix
	m <- check_rank(x_ordered, n, m, verbose)
	1:n %in% order_index[1:m]
}

init_v1 <- function(x, m, verbose)
{
	n <- NROW(x)
	# Mahalanobis distances about the mean
	dist <- mahalanobis(x, colMeans(x), cov(x))
 	# sort by dist
	order_index <- order(dist)
	x_ordered <- x[order_index, , drop = FALSE]
	# check for rank deficiency of scatter matrix
	m <- check_rank(x_ordered, n, m, verbose)
	1:n %in% order_index[1:m]
}

init_v2 <- function(x, m, verbose)
{
	n <- NROW(x)
	center <- apply(x, 2L, quantile, probs = 0.5, type = 1)
	x_centered <- sweep(x, 2L, center)
	# Euclidean norm
	dist <- apply(x_centered, 1, function(u) sqrt(sum(u^2)))
 	# sort by dist
	order_index <- order(dist)
	x_ordered <- x[order_index, , drop = FALSE]
	# check for rank deficiency of scatter matrix
	m <- check_rank(x_ordered, n, m, verbose)
	1:n %in% order_index[1:m]
}

check_rank <- function(x_ordered, n, m, verbose)
{
	p <- NCOL(x_ordered) 
	sigma <- cov(x_ordered[1:m, ])	
	rnk <- qr(sigma)$rank
	while (m < n & rnk < p) {
		if (verbose)
			cat("rank deficient\n")
		m <- m + 1
		rnk <- qr(cov(x_ordered[1:m, ]))$rank
	}
	return(m)	
}

mahalanobis_distance <- function(x, subset)
{
	x_subset <- x[subset, , drop = FALSE]
	center <- colMeans(x_subset)
	x_centerd <- sweep(x_subset, 2L, center)
	scatter <- crossprod(x_centerd) / (sum(subset) - 1)
	sqrt(mahalanobis(x, center, scatter))
}

chi2_cutoff <- function(alpha, n, p, r)
{
	h <- (n + p + 1) / 2
	chr <- max(0, (h - r) / (h + r))
	cnp <- 1 + (p + 1) / (n - p) + 2 / (n - 1 - 3*p)
	cnpr <- cnp + chr
	cnpr * sqrt(qchisq(alpha / n, p, lower.tail = FALSE))
}


#----------------------------------------------------------------------------
chi2Crit <- function(n, p, r, alpha) {
    h <- (n + p + 1)/2
    chr <- max(0, (h - r)/(h + r))
    cnp <- 1 + (p + 1)/(n - p) + 2/(n - 1 - 3 * p)
    cnpr <- cnp + chr
    cnpr^2 * qchisq(alpha/n, p, lower.tail = FALSE)
}


mvBACON2 <- function (x, collect = 4, m = min(collect * p, n * 0.5),
    alpha = 0.95, init.sel = c("Mahalanobis", "dUniMedian", "random",
        "manual", "V2"), man.sel, maxsteps = 100,
    allowSingular = FALSE)
{
    init.sel <- match.arg(init.sel)
    n <- nrow(x)
    p <- ncol(x)
    stopifnot(n > p, p > 0, 0 < alpha, alpha < 1)
    ordered.indices <- switch(init.sel,
        Mahalanobis = {
            order(mahalanobis(x, center = colMeans(x), cov = var(x)))
        },
        random = sample(n),
        manual = {
            m <- length(man.sel)
            stopifnot(is.numeric(man.sel), 1 <= man.sel, man.sel <= n,
                man.sel == as.integer(man.sel))
            c(man.sel, c(1:n)[-man.sel])
        },
        dUniMedian = {
            x.centr <- sweep(x, 2, colMedians(x))
            order(mahalanobis(x.centr, 0, cov(x.centr)))
        },
        V2 = {
            x.centr <- sweep(x, 2, colMedians(x))
            order(apply(x.centr, 1, crossprod))
        },
        stop("invalid 'init.sel' -- should not happen; please report!")
    )
    m <- as.integer(m)
    stopifnot(n >= m, m > 0)
    ordered.x <- x[ordered.indices, , drop = FALSE]

    while (m < n) {
        rnk <- qr(var(ordered.x[1:m, , drop = FALSE]))$rank
        if (rnk == p)
            break
        else
            m <- m + 1
    }

    if (rnk < p && !allowSingular)
        stop("matrix-rank ( x[1:m,] ) < p  for all m <= n")
    subset <- 1:n %in% ordered.indices[1:m]
    presubset <- rep(FALSE, n)
    converged <- FALSE
    steps <- 1L
    repeat {
        r <- sum(subset)
        x. <- x[subset, , drop = FALSE]
        center <- colMeans(x.)
        cov <- var(x.)
        dis <- mahalanobis(x, center, cov)
        converged <- !any(xor(presubset, subset))
        if (converged)
            break
        presubset <- subset
        limit <- chi2Crit(n, p, r, alpha)
        subset <- dis < limit
        steps <- steps + 1L
        if (steps > maxsteps)
            break
    }
    if (steps > maxsteps)
        warning("basic subset has not converged in ", maxsteps, " steps")
    list(dis = dis, subset = subset, center = center, cov = cov,
        steps = steps, limit = limit, converged = converged)
}
