#===============================================================================
# SUBJECT  Benchmarking of the implementation against the R package robustX
# AUTHORS  Tobias Schoch, tobias.schoch@gmail.com, May 15, 2021
# LICENSE  GPL >= 2
# COMMENT
#===============================================================================
library(wbacon)
library(robustX)
# check that version 1.25 (or newer) of robustX is installed
robustX_version <- as.numeric(gsub("-", "", getNamespaceVersion("robustX")))
if (robustX_version < 1.25) {
    stop(paste0("Version >= 1.25 of package 'robustX' is required, you have ",
        robustX_version, "\nGrab the new version from R-Forge\n"))
}
library(microbenchmark)
#===============================================================================
# 0 SETUP
#===============================================================================
VERSION   <- "V2"    # Initialization (version "V1" or "V2")
ALPHA     <- 0.05    # Cutoff value: (1-ALPHA) quantile of the chi square distr.
N_THREADS <- 1       # OpenMP number of threads

#===============================================================================
# 1 BACON robust linear regression
#===============================================================================
# scenarios
#-------------------------------------------------------------------------------
scenario <- data.frame(matrix(c(
#---------------------------------
#       p         n   eps   times
#---------------------------------
       5,       100,  0.05,    10,
       10,      100,  0.05,    10,
       20,      100,  0.05,    10,
       5,      1000,  0.05,    10,
       10,     1000,  0.05,    10,
       20,     1000,  0.05,    10),
    byrow = TRUE, ncol = 4))
colnames(scenario) <- c("p", "n", "eps", "times")
# n: number observations; p: number of variables; eps: amount of contamination
# times: number of microbenchmark evaluations

#-------------------------------------------------------------------------------
# microbenchmark
#-------------------------------------------------------------------------------
scenario$ratio <- NA
for (i in 1:NROW(scenario)) {
    p <- scenario[i, "p"]
    n <- scenario[i, "n"]
    epsilon <- scenario[i, "eps"]
    times <- scenario[i, "times"]
    set.seed(1)

    # outlier-free ('good') and contaminated ('bad') data
    good_obs <- matrix(rnorm(p * floor(n * (1 - epsilon)), 0, 1), ncol = p)
    bad_obs <- matrix(rnorm(p * floor(n * epsilon), 4, 1), ncol = p)
    x <- rbind(bad_obs, good_obs)
    beta <- rep(1, p)
    y_good <- good_obs %*% beta
    y_bad <- bad_obs %*% (beta + 10)
    y <- c(y_bad, y_good) + rnorm(n)
    dat <- data.frame(y, x)

    # benchmarking
    m0 <- microbenchmark(BACON(x, y, alpha = ALPHA, init.sel = VERSION,
        verbose = FALSE), times = scenario$times[i])
    m1 <- microbenchmark(wBACON_reg(y ~ ., data = dat, alpha = ALPHA,
        version = VERSION, n_threads = N_THREADS), times = scenario$times[i])
    scenario$ratio[i] <- mean(m0$time) / mean(m1$time)

    cat(i, "of", NROW(scenario), "\n")
    flush.console()
}

# The variable 'ratio' measures by how many times wBACON_reg is faster than
# robustX:BACON
scenario

#===============================================================================
# 2 BACON multivariate outlier nomination (detection)
#===============================================================================
# scenarios
#-------------------------------------------------------------------------------
scenario <- data.frame(matrix(c(
#------------------------------------
#         p        n     eps  times
#------------------------------------
          5,    1000,   0.05,   50,
         10,    1000,   0.05,   50,
         20,    1000,   0.05,   50,
         30,    1000,   0.05,   50,
         40,    1000,   0.05,   50,
         50,    1000,   0.05,   50,
        100,    1000,   0.05,   50,
        200,    1000,   0.05,   50,
          5,   10000,   0.05,   30,
         10,   10000,   0.05,   30,
         20,   10000,   0.05,   30,
         30,   10000,   0.05,   30,
         40,   10000,   0.05,   30,
         50,   10000,   0.05,   30,
        100,   10000,   0.05,   30,
        200,   10000,   0.05,   30,
          5,  100000,   0.05,   20,
         10,  100000,   0.05,   20,
         20,  100000,   0.05,   20,
         30,  100000,   0.05,   20,
         40,  100000,   0.05,   20,
         50,  100000,   0.05,   20,
        100,  100000,   0.05,   20,
        200,  100000,   0.05,   20,
          5, 1000000,   0.05,   10,
         10, 1000000,   0.05,   10,
         20, 1000000,   0.05,   10,
         30, 1000000,   0.05,   10,
         40, 1000000,   0.05,   10,
         50, 1000000,   0.05,   10,
        100, 1000000,   0.05,   10,
        200, 1000000,   0.05,   10),
    byrow = TRUE, ncol = 4))
colnames(scenario) <- c("p", "n", "eps", "times")
# n: number observations; p: number of variables; eps: amount of contamination
# times: number of microbenchmark evaluations

#-------------------------------------------------------------------------------
# microbenchmark
#-------------------------------------------------------------------------------
scenario$ratio <- NA
for (i in 1:NROW(scenario)) {
    p <- scenario[i, "p"]
    n <- scenario[i, "n"]
    epsilon <- scenario[i, "eps"]
    times <- scenario[i, "times"]
    set.seed(1)
    good_obs <- matrix(rnorm(p * floor(n * (1 - epsilon)), 0, 1), ncol = p)
    bad_obs <- matrix(rnorm(p * floor(n * epsilon), 4, 1), ncol = p)
    x <- rbind(bad_obs, good_obs)

    # benchmarking
    m0 <- microbenchmark(BACON(x, alpha = ALPHA, init.sel = VERSION,
        verbose = FALSE), times = scenario$times[i])
    m1 <- microbenchmark(wBACON(x, alpha = ALPHA, version = VERSION,
        n_threads = N_THREADS), times = scenario$times[i])
    scenario$ratio[i] <- mean(m0$time) / mean(m1$time)

    cat(i, "of", NROW(scenario), "\n")
    flush.console()
}

# The variable 'ratio' measures by how many times wBACON is faster than
# robustX:BACON
scenario
