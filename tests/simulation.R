#===============================================================================
# SUBJECT  Replication of the simulation of Billor et al. 2000
# AUTHORS  Tobias Schoch, tobias.schoch@gmail.com, March 5, 2021
# LICENSE  GPL >= 2
# COMMENT
#===============================================================================
library(wbacon)
# simulation using the contamination model of Billor et al. (2000, p. 290):
# F(x) = (1-epsilon) * N(0, I_p) + epsilon * N(4, I_p),
# where I_p is the p-dimensional identity matrix, and N is the c.d.f. of the
# Gaussian distr.

criteria <- function(object, epsilon)
{
	n <- object$n
	# relative number of nominated outliers (perfect performance: 1.0)
	Out <- sum(object$subset == 0)

	# relative number of correctly identified outliers (perfect performance
	# 1.0)
	TrueOut <- sum(object$subset[1:floor(n * epsilon)] == 0)

	if (epsilon > 0)
		reference <- floor(epsilon * n)
	else
		reference <- n

	Out <- Out / reference
	TrueOut <- TrueOut / reference
	c(Out = Out, TrueOut = TrueOut)
}

p <- 5
n <- 500
replicates <- 500
epsilon <- 0.1

set.seed(1)
res <- NULL
for (i in 1:replicates){
	good_obs <- matrix(rnorm(p * floor(n * (1 - epsilon)), 0, 1), ncol = p)
	bad_obs <- matrix(rnorm(p * floor(n * epsilon), 4, 1), ncol = p)
	if (NROW(bad_obs) != 0) {
		all <- rbind(bad_obs, good_obs)
	} else {
		all <- good_obs
	}

	m <- wBACON(all, alpha = 0.95)
	res <- rbind(res, criteria(m, epsilon))
}
result <- cbind(apply(res, 2, mean), apply(res, 2, sd), apply(res, 2, max))
colnames(result) <- c("mean", "sd", "max")
result


