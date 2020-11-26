setwd("C:/My/code/cbacon/src")
dyn.load("cbacon.dll")

# set.seed(1)
# x <- matrix(c(rnorm(20, 10, 2)), ncol = 2)
# n <- nrow(x); p <- ncol(x)   
# w <- rep(1, n)
#
# x[10,1] <- 100

# l <- .C("weightedmean", x = as.double(x), 
#    w = as.double(w), 
#    center = as.double(numeric(p)), 
#    n = as.integer(n), 
#    p = as.integer(p))
#
# s <- .C("weightedscatter", x = as.double(x), 
#    w = as.double(w), 
#    center = as.double(l$center), 
#    scatter = as.double(numeric(p * p)), 
#    n = as.integer(n), 
#    p = as.integer(p))
#
# m <- .C("mahalanobis", x = as.double(x), 
#    center = as.double(l$center), 
#    scatter = as.double(s$scatter),
#    dist = as.double(numeric(n)), 
#    n = as.integer(n), 
#    p = as.integer(p))
#
#
# mahalanobis(x, colMeans(x), var(x))
#
#
# xcenter <- sweep(x, 2, colMeans(x))
# L <- t(chol(var(x)))
# colSums((solve(L) %*% t(xcenter))^2)
#


# s <- .C("weightedmeancov", x = as.double(x), 
#    w = as.double(w), 
#    center = as.double(numeric(p)), 
#    scatter = as.double(numeric(p * p)), 
#    n = as.integer(n), 
#    p = as.integer(p))


WB <- function(x, w = NULL, alpha = 0.95, maxiter = 50, verbose = FALSE){ 
   x <- as.matrix(x)
   n <- nrow(x); p <- ncol(x)
   if (is.null(w)) w <- rep(1, n)
   stopifnot(n > p, p > 0,  0 < alpha, alpha < 1, n == length(w))
   tmp <- .C("wbacon", x = as.double(x), w = as.double(w), 
      center = as.double(numeric(p)), scatter = as.double(numeric(p * p)),
      dist = as.double(numeric(n)), n = as.integer(n), p = as.integer(p),
      alpha = as.double(alpha), subset = as.integer(numeric(n)), 
      cutoff = as.double(numeric(1)), maxiter = as.integer(abs(maxiter)),
      verbose = as.integer(verbose))
   tmp$scatter <- matrix(tmp$scatter, ncol = p)
   tmp$verbose <- NULL
   tmp$converged <- ifelse(tmp$maxiter < maxiter, TRUE, FALSE)
   names(tmp$center) <- colnames(x)
   colnames(tmp$scatter) <- colnames(x)
   rownames(tmp$scatter) <- colnames(x)
   class(tmp) <- "robmv"
   tmp
}

# print.robmv <- function(x, ...){
#    cat("Robust estimates of location and scatter\n\n")   
#    cat("Location\n")
#    print(x$center)
#    cat("\nScatter\n")
#    print(x$scatter)
#    cat(paste0("\nDeclared outliers: ", x$n - sum(x$subset), " (sample size n = ", 
#       x$n,")\n"))
# }
library(modi)
data(bushfire)
w <- wBACON(bushfire, rep(1,38), alpha = 0.99, verbose = T)

#library(robustX)

library(cellWise)
data(philips)
wBACON(philips, verbose = T)

o <- mvBACON(as.matrix(philips), alpha = 0.99, init.sel = "dUniMedian")




