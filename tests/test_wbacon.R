#===============================================================================
# SUBJECT  Test the implementation of 'wbacon' against the function BACON in the
#          package 'robustX'
# AUTHORS  Tobias Schoch, tobias.schoch@gmail.com, January 23, 2021
# LICENSE  GPL >= 2
# COMMENT  pkg 'robustbase' must be installed
#===============================================================================
library(robustX)
library(wbacon)

#===============================================================================
# Comparison function
#===============================================================================
# this function compares our implementation with the reference implementation
# in pkg 'robustX'
compare <- function(data, name, init = "V2"){
	# we test 11 different values for the parameter 'alpha'
	set_of_alphas <- c(seq(0.3, 0.9, 0.1), 0.95, 0.99, 0.999, 0.9999)
	data <- as.matrix(data)

	res <- NULL
	acc <- sqrt(.Machine$double.eps)
	for (i in 1:length(set_of_alphas)) {

		# reference implementation robustX::BACON with the "V2" (version 2)
		# intialization method of Billor et al. (2000)
		ref.init <- ifelse(init == "V2", "V2", "Mahalanobis")
		reference <- BACON(data, alpha = set_of_alphas[i], init.sel = ref.init,
			verbose = FALSE)

		# our implementation
		new_imp <- wBACON(data, alpha = set_of_alphas[i], version = init)
		deviations <- NULL

		# absolute deviations
		dev_location <- max(abs(reference$center - new_imp$center))
		dev_scatter <- max(abs(reference$cov - new_imp$cov))
		dev_distance <- max(abs(reference$dis - new_imp$dist))

		# deviations in terms of the selected subset
		dev_subset <- sum(reference$subset - new_imp$subset)

		at <- matrix(rep(0, 4), ncol = 4)
		rownames(at) <- paste(set_of_alphas[i])

		# throw an error if the deviations are larger than 'acc'
		if (dev_location > acc)
			at[1] <- dev_location
		if (dev_scatter > acc)
			at[2] <- dev_scatter
		if (dev_distance> acc)
			at[3] <- dev_distance
		if (dev_subset> acc)
			at[4] <- dev_subset

		if (any(at > 0))
		res <- rbind(res, at)
	}

	if (!is.null(res)) {
		cat(name, ": differences detected\n")
		colnames(res) <- c("center", "cov", "distance", "subset")
	}

	res
}

#===============================================================================
# Tests I
#===============================================================================
# We test the implementations on 7 well known data sets. The function 'compare'
# is expected to return 'NULL'; otherwise the two implementations differ

data(hbk, package = "robustbase")
d_hbk <- data.matrix(hbk[, 1:3])
compare(d_hbk, "hbk", "V1")
compare(d_hbk, "hbk", "V2")

data(bushfire, package = "robustbase")
compare(bushfire, "bushfire", "V1")
compare(bushfire, "bushfire", "V2")

data(aircraft, package = "robustbase")
d_aircraft <- data.matrix(aircraft[, 1:4])
compare(d_aircraft, "aircraft", "V1")
compare(d_aircraft, "aircraft", "V2")

data(education, package = "robustbase")
d_education<- data.matrix(education[, 2:4])
compare(d_education, "education", "V1")
compare(d_education, "education", "V2")

data(heart, package = "robustbase")
d_heart <- data.matrix(heart[, 1:2])
compare(d_heart, "heart", "V1")
compare(d_heart, "heart", "V2")

data(milk, package = "robustbase")
compare(milk, "milk", "V1")
compare(milk, "milk", "V2")

data(pulpfiber, package = "robustbase")
pulp <- as.matrix(pulpfiber)
compare(pulp, "pulpfiber", "V1")
compare(pulp, "pulpfiber", "V2")

#===============================================================================
# Tests II
#===============================================================================
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

#plot(m$dist, col = 2 - m$subset, pch = 19)
#abline(v = floor(epsilon * n))





library(microbenchmark)
#-------------------------------------------------------------------------------
p <- 10
n <- 10000
epsilon <- 0.1
set.seed(1)
good_obs <- matrix(rnorm(p * floor(n * (1 - epsilon)), 0, 1), ncol = p)
bad_obs <- matrix(rnorm(p * floor(n * epsilon), 4, 1), ncol = p)
all <- rbind(bad_obs, good_obs)
microbenchmark(wBACON(all, alpha = 0.95), times = 100)
#Unit: milliseconds
#                      expr     min      lq     mean   median      uq     max neval
# wBACON(all, alpha = 0.95) 12.7512 13.6846 13.94739 13.94055 14.2525 15.2571   100

#-------------------------------------------------------------------------------
p <- 20
n <- 100000
epsilon <- 0.1
set.seed(1)
good_obs <- matrix(rnorm(p * floor(n * (1 - epsilon)), 0, 1), ncol = p)
bad_obs <- matrix(rnorm(p * floor(n * epsilon), 4, 1), ncol = p)
all <- rbind(bad_obs, good_obs)
microbenchmark(wBACON(all, alpha = 0.95), times = 10)
#Unit: milliseconds
#                      expr      min       lq     mean   median       uq      max neval
# wBACON(all, alpha = 0.95) 414.3206 415.9524 419.2129 417.3597 422.7726 430.0953    10

#-------------------------------------------------------------------------------
p <- 200
n <- 10000
epsilon <- 0.1
set.seed(1)
good_obs <- matrix(rnorm(p * floor(n * (1 - epsilon)), 0, 1), ncol = p)
bad_obs <- matrix(rnorm(p * floor(n * epsilon), 4, 1), ncol = p)
all <- rbind(bad_obs, good_obs)
microbenchmark(wBACON(all, alpha = 0.95), times = 10)
#Unit: seconds
#                      expr      min       lq     mean   median       uq      max neval
# wBACON(all, alpha = 0.95) 2.224406 2.232228 2.238915 2.236654 2.243672 2.262959    10

#-------------------------------------------------------------------------------
p <- 20
n <- 1000000
epsilon <- 0.1
set.seed(1)
good_obs <- matrix(rnorm(p * floor(n * (1 - epsilon)), 0, 1), ncol = p)
bad_obs <- matrix(rnorm(p * floor(n * epsilon), 4, 1), ncol = p)
all <- rbind(bad_obs, good_obs)
microbenchmark(wBACON(all, alpha = 0.95), times = 10)
#Unit: seconds
#                      expr      min       lq     mean  median       uq      max neval
# wBACON(all, alpha = 0.95) 5.216988 5.249868 5.304857 5.27994 5.352519 5.462912    10

