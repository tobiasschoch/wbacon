#===============================================================================
# SUBJECT  Benchmarking of the implementation
# AUTHORS  Tobias Schoch, tobias.schoch@gmail.com, March 5, 2021
# LICENSE  GPL >= 2
# COMMENT
#===============================================================================
library(wbacon)
library(microbenchmark)
VERS <- "V2"
#-------------------------------------------------------------------------------
cat("----------------------------------------------------------\n")
cat("negative: faster; positive: slower\n")
cat("----------------------------------------------------------\n")

p <- 10
n <- 10000
epsilon <- 0.1
set.seed(1)
good_obs <- matrix(rnorm(p * floor(n * (1 - epsilon)), 0, 1), ncol = p)
bad_obs <- matrix(rnorm(p * floor(n * epsilon), 4, 1), ncol = p)
all <- rbind(bad_obs, good_obs)
res <- microbenchmark(wBACON(all, alpha = 0.95, version = VERS), times = 100)
avg <- mean(res$time / 1e6)
med <- median(res$time / 1e6)
ref_avg <- 12.446
ref_med <- 12.290
cat("----------------------------------------------------------\n")
cat("p = 10, n = 10'000\n")
cat("average:", 100 * (avg - ref_avg) / ref_avg, "%\n")
cat("median:", 100 * (med - ref_med) / ref_med, "%\n\n")

#-------------------------------------------------------------------------------
p <- 20
n <- 100000
epsilon <- 0.1
set.seed(1)
good_obs <- matrix(rnorm(p * floor(n * (1 - epsilon)), 0, 1), ncol = p)
bad_obs <- matrix(rnorm(p * floor(n * epsilon), 4, 1), ncol = p)
all <- rbind(bad_obs, good_obs)
res <- microbenchmark(wBACON(all, alpha = 0.95, version = VERS), times = 10)
avg <- mean(res$time / 1e6)
med <- median(res$time / 1e6)
ref_avg <- 324.580
ref_med <- 323.581
cat("----------------------------------------------------------\n")
cat("p = 20, n = 100'000\n")
cat("average:", 100 * (avg - ref_avg) / ref_avg, "%\n")
cat("median:", 100 * (med - ref_med) / ref_med, "%\n\n")

#-------------------------------------------------------------------------------
p <- 200
n <- 10000
epsilon <- 0.1
set.seed(1)
good_obs <- matrix(rnorm(p * floor(n * (1 - epsilon)), 0, 1), ncol = p)
bad_obs <- matrix(rnorm(p * floor(n * epsilon), 4, 1), ncol = p)
all <- rbind(bad_obs, good_obs)
res <- microbenchmark(wBACON(all, alpha = 0.95, version = VERS), times = 10)
avg <- mean(res$time / 1e6)
med <- median(res$time / 1e6)
ref_avg <- 1393.797
ref_med <- 1383.028
cat("----------------------------------------------------------\n")
cat("p = 200, n = 10'000\n")
cat("average:", 100 * (avg - ref_avg) / ref_avg, "%\n")
cat("median:", 100 * (med - ref_med) / ref_med, "%\n\n")

#-------------------------------------------------------------------------------
p <- 20
n <- 1000000
epsilon <- 0.1
set.seed(1)
good_obs <- matrix(rnorm(p * floor(n * (1 - epsilon)), 0, 1), ncol = p)
bad_obs <- matrix(rnorm(p * floor(n * epsilon), 4, 1), ncol = p)
all <- rbind(bad_obs, good_obs)
res <- microbenchmark(wBACON(all, alpha = 0.95, version = VERS), times = 10)
avg <- mean(res$time / 1e6)
med <- median(res$time / 1e6)
ref_avg <- 4269.193
ref_med <- 4229.532
cat("----------------------------------------------------------\n")
cat("p = 20, n = 1'000'000\n")
cat("average:", 100 * (avg - ref_avg) / ref_avg, "%\n")
cat("median:", 100 * (med - ref_med) / ref_med, "%\n\n")

#-------------------------------------------------------------------------------
cat("\nEOF\n")
