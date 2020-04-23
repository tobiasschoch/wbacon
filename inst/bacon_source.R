


bc <- function(x, w = NULL, collect = 4, m = min(collect * p, n * 0.5), alpha = 0.95,
		    maxsteps = 100){
   x <- as.matrix(x)
   n <- nrow(x); p <- ncol(x)
   stopifnot(n > p, p > 0,  0 < alpha, alpha < 1)
   m <- as.integer(m)
   stopifnot(n >= m, m > 0)
   if(is.null(w)) w <- rep(1, n)

   # initial center and scatter
   x_centered <- sweep(x, 2, apply(x, 2, quantile, probs = 0.5, type = 1), check.margin = FALSE)

   # mahalanobis distance 
   # ordered.indices <- order(rowSums(x_centered %*% cov(x_centered) * x_centered))
   # ordered.indices <- order(mahalanobis(x_centered, 0, cov(x_centered)))
   ordered.indices <- order(mahalanobis(x_centered, 0, t(x_centered) %*% x_centered / (n-1)))


   ordered.x <- x[ordered.indices, , drop = FALSE]
   while (m < n && p > (rnk <- qr(var(ordered.x[1:m, , drop = FALSE]))$rank)){
      print(rnk) 
      m <- m + 1L
   }
   if(rnk < p ) stop("matrix not of full rank\n")



   subset <- 1:n %in% ordered.indices[1:m]
   presubset <- rep(FALSE, n)
   converged <- FALSE; steps <- 1L
   repeat {
      r <- sum(subset)
	 x_subset <- x[subset, , drop = FALSE]
	 w_subset <- w[subset]
	 sum_w <- sum(w_subset)
	 
	 center <- colSums(x_subset * w_subset) / sum_w 
	 x_sweeped <- sweep(x, 2L, center, check.margin = FALSE)
	 tmp <- x_sweeped[subset, , drop = FALSE]
	 scatter <- crossprod(tmp * w_subset, tmp) / sum_w
	 # mahalanobis distance 
	 dis <- sqrt(rowSums(tcrossprod(x_sweeped, solve(scatter)) * x_sweeped))

	 # convergence
	 converged <- !any(xor(presubset, subset))
	 if(converged)
	    break

	 presubset <- subset

	 # cut-off point (chi-squared)
	 h <- (n + p + 1) / 2
	 chr <- max(0, (h - r) / (h + r))
	 cnp <- 1 + (p + 1) / (n - p) + 2 / (n - 1 - 3*p)
	 cnpr <- cnp + chr
	 limit <- cnpr * sqrt(qchisq(alpha / n, p, lower.tail = FALSE))

	 subset <- dis < limit
	 steps <- steps + 1L
	 if (steps > maxsteps)
	    break
   }
   if(steps > maxsteps)
      warning("basic subset has not converged in ", maxsteps, " steps")

   list(dis = dis, subset = subset, center = center, scatter = scatter,
         steps = steps, limit = limit, converged = converged)
} 


collect = 4; m = min(collect * p, n * 0.5); alpha = 0.99; maxsteps = 100

