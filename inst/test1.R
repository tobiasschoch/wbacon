setwd("C:/My/code/cbacon/src")
dyn.load("cbacon.dll")

library(robustbase)
library(robustX)

test_wBACON <- function(x, w = NULL, alpha = 0.95, intercept = FALSE, 
   na.rm = FALSE, maxiter = 50, verbose = FALSE) 
{ 
   if (!is.matrix(x))
      x <- as.matrix(x)

   if (intercept) {
      x <- x[, - which(apply(x, 2, var) == 0)]
   }

   n <- nrow(x); p <- ncol(x)

   if (is.null(w)) w <- rep(1, n)

   stopifnot(n > p, p > 1,  0.1 <= alpha, alpha < 1, n == length(w))

   # NA treatment
   cc <- stats::complete.cases(x, w) 
   if (sum(cc) != n) {
      if (na.rm) { 
	 x <- x[cc, ]
	 w <- w[cc] 
      } else 
	 stop("Data must not contain missing values; see argument 'na.rm'\n", 
	    call. = FALSE)
   } 
   n <- nrow(x)

   # check if any element is not finite 
   chk <- sum(is.finite(c(x, w))) != (1 + p) * n 
   if (chk) 
      stop("Some observations are not finite\n", call. = FALSE)
   
   # compute weighted BACON algorithm
   tmp <- .C("wbacon", x = as.double(x), w = as.double(w), 
      center = as.double(numeric(p)), scatter = as.double(numeric(p * p)),
      dist = as.double(numeric(n)), n = as.integer(n), p = as.integer(p),
      alpha = as.double(alpha), subset = as.integer(rep(0, n)), 
      cutoff = as.double(numeric(1)), maxiter = as.integer(abs(maxiter)),
      verbose = as.integer(verbose))

   tmp$scatter <- matrix(tmp$scatter, ncol = p)
   tmp$verbose <- NULL
   tmp$converged <- ifelse(tmp$maxiter < maxiter, TRUE, FALSE)
   tmp$call <- match.call()
   names(tmp$center) <- colnames(x)
   colnames(tmp$scatter) <- colnames(x)
   rownames(tmp$scatter) <- colnames(x)
   tmp
}

# compare our implementation with the one of Ueli Oetliker (robustX)
compare <- function(data, name){
   acc <- sqrt(.Machine$double.eps)
   set_of_alphas <- c(seq(0.1, 0.9, 0.1), 0.95, 0.99, 
      0.999, 0.9999)
   res <- NULL 

   if (!is.matrix(data))
      data <- as.matrix(data)

   for (i in 1:length(set_of_alphas)) {
      oetliker <- BACON(data, alpha = set_of_alphas[i], verbose = FALSE,
	 init.sel = "dUniMedian") 
      my <- test_wBACON(data, alpha = set_of_alphas[i])
      deviations <- NULL

      dev_location <- max(abs(oetliker$center - my$center))
      dev_scatter <- max(abs(oetliker$cov - my$scatter))
      dev_distance <- max(abs(oetliker$dis - my$dist))
      dev_subset <- sum(oetliker$subset - my$subset)

      at <- matrix(rep(0, 4), ncol = 4)
      rownames(at) <- paste(set_of_alphas[i])

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

   if (!is.null(res))
      colnames(res) <- c("center", "scatter", "distance", "subset")

   res
}

#===============================================================================
data(hbk)
d_hbk <- data.matrix(hbk[, 1:3])
compare(d_hbk, "hbk")

#===============================================================================
data(bushfire)
compare(bushfire, "bushfire")

#===============================================================================
data(aircraft)
d_aircraft <- data.matrix(aircraft[, 1:4])
compare(d_aircraft, "aircraft")

#===============================================================================
data(education)
d_education<- data.matrix(education[, 2:4])
compare(d_education, "education")

#===============================================================================
data(heart)
d_heart <- data.matrix(heart[, 1:2])
compare(d_heart, "heart")

#===============================================================================
data(milk)
compare(milk, "milk")

#===============================================================================
data(pulpfiber)
d_pulp <- as.matrix(pulpfiber)
compare(d_pulp, "pulpfiber")



