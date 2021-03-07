#===============================================================================
# SUBJECT  Benchmarking of the implementation
# AUTHORS  Tobias Schoch, tobias.schoch@gmail.com, March 5, 2021
# LICENSE  GPL >= 2
# COMMENT
#===============================================================================
library(wbacon)
library(microbenchmark)

ref <- read.csv("data_benchmark.csv")
#-------------------------------------------------------------------------------
# TUNING OF THE METHOD
VERS <- "V2"
ALPHA <- 0.95

#-------------------------------------------------------------------------------
# SCENARIOS
scenario <- data.frame(matrix(c(
#------------------------------
#     p         n   eps   times
#------------------------------
     10,    10000,  0.1,  100,
     20,   100000,  0.1,  100,
    200,    10000,  0.1,   10,
     20,  1000000,  0.1,   10), byrow = TRUE, ncol = 4))
colnames(scenario) <- c("p", "n", "eps", "times")

#-------------------------------------------------------------------------------
# RUN BENCHMARK
cat("----------------------------------------------------------\n")
cat("negative: faster; positive: slower\n")
cat("----------------------------------------------------------\n")

res <- matrix(nrow = NROW(scenario), ncol = 2)
for (i in 1:NROW(scenario)) {
    p <- scenario[i, "p"]
    n <- scenario[i, "n"]
    epsilon <- scenario[i, "eps"]
    times <- scenario[i, "times"]
    set.seed(1)
    good_obs <- matrix(rnorm(p * floor(n * (1 - epsilon)), 0, 1), ncol = p)
    bad_obs <- matrix(rnorm(p * floor(n * epsilon), 4, 1), ncol = p)
    all <- rbind(bad_obs, good_obs)

    this <- microbenchmark(wBACON(all, alpha = ALPHA, version = VERS),
        times = times)

    res[i, 1] <- mean(this$time / 1e6)
    res[i, 2] <- median(this$time / 1e6)
    cat("----------------------------------------------------------\n")
    cat(paste0("p = ", p, ", n = ", n, "\n"))
    cat("average:", round(100 * (res[i, 1] - ref[i, 1]) / ref[i, 1], 2), "%\n")
    cat("median:", round(100 * (res[i, 2] - ref[i, 2]) / ref[i, 2], 2), "%\n\n")
}

write.csv(res, file = "new.csv", row.names = FALSE)

#-------------------------------------------------------------------------------
cat("\nEOF\n")
