library(robsurvey)
setwd("C:/My/code/cbacon/src")
dyn.load("wbacon.dll")
#dyn.unload("wbacon.dll")

data(iris)
at <- c(27, 139, 96, 55, 128, 116, 90, 13, 36, 23, 34, 91, 147, 50, 17)
#at <- 1:150
x <- as.matrix(iris[at, 1:3])
n <- nrow(x)
x <- cbind(x, rep(1, n)) 
y <- iris[at, 4]
p <- ncol(x)

#set.seed(1)
#w <- abs(rnorm(n))

w <- rep(1, n)

bc <- wBACON(x, w, intercept = TRUE)
dist <- bc$dist
subset0 <- bc$subset
m <- sum(subset0)

#-----------------------
#subset0 <- x[, 3] == 1.4
#m <- sum(subset0)

#-------------------------------------------------------------------------------
library(robustX)
setwd("C:/My/code/cbacon/src")
dyn.load("wbacon.dll")
data(iris)
x <- as.matrix(iris[, 1:3])
y <- iris[,4]
w <- rep(1,150)

set.seed(1)
w <- w + abs(rnorm(length(y), 3, 2))

tmp <- mvBACON(x)
source("C:/My/code/cbacon/inst/baconreg.R")

#a <- robustX:::.lmBACON(x, y, init.dis = tmp$dis)
b <- wlmBACON(x, y, w, tmp, alg5 = FALSE)
z <- foo(x, y, w, tmp)



foo <- function(x, y, w, bc, intercept = TRUE, alpha = 0.95, collect = 4,
	verbose = TRUE)
{
	n <- length(y)
	if (intercept)
		x <- cbind(rep(1, n), x)
	p <- ncol(x)
	m <- sum(bc$subset)
	collect <- min(collect, floor(n / p))
	tmp <- .C("wbacon_reg", x = as.double(x), y = as.double(y), w = as.double(w),
		resid = as.double(numeric(n)), beta = as.double(numeric(p)), 
		subset = as.integer(bc$subset), dist = as.double(bc$dis), n = as.integer(n), 
		p = as.integer(p), m = as.integer(m), sucess = as.integer(1), 
		verbose = as.integer(verbose), collect = as.integer(collect), 
		alpha = as.double(alpha), maxiter = as.integer(50))
	tmp
}

lm.fit(x[tmp$subset == 1, ], y[tmp$subset == 1])$coefficients


lm.wfit(x[subset0 == 1, ], y[subset0 == 1], w[subset0 == 1])$coefficients


w <- 1:15
wx <- sqrt(w[1:10]) * x[1:10, ]
L <- t(chol(crossprod(wx)))
chol_update(L, x[11,])


t(chol(crossprod(x[1:11,])))


L <- t(chol(crossprod(x)))
xty <- crossprod(x, y)
beta <- lm.fit(x,y)$coefficients

x <- cbind(x, rep(1, 150), rep(2, 150))

at <- 1:11
L <- t(chol(crossprod(x[at, ])))
xty <- crossprod(x[at, ], y[at])
subset0 <- rep(0, n); subset0[at] <- 1
m <- sum(subset0)
subset1 <- rep(0, n); subset1[c(at, 12)] <- 1


dist <- rep(1, n)



# tmp <- .C("hat_matrix", L = as.double(L), x = as.double(x), hat = as.double(numeric(n)),
#    work_pp = as.double(numeric(p*p)), work_np = as.double(numeric(n*p)),
#    n = as.integer(n), p = as.integer(p))
#
# tmp <- .C("update_chol_xty", x = as.double(x), y = as.double(y), 
#    xty = as.double(xty), L = as.double(L), subset0 = as.integer(subset0), 
#    subset1 = as.integer(subset1), work_p = as.double(numeric(p)),
#    n = as.integer(n), p = as.integer(p))
#
# tmp <- .C("full_rank_subset", L = as.double(L), x = as.double(x), dist = as.double(dist), 
#    subset = as.integer(subset0), iarray = as.integer(numeric(n)), 
#    work = as.double(numeric(p)), n = as.integer(n), p = as.integer(p), 
#    m = as.integer(m))
#


Algorithm 4
   + check if x has full rank
      -> add obs. with smallest dist (from Alg 3) until x has full rank
   + solve for b in Xb=y
   + compute t[i]
   + select m <- (p+1) obs. with smallest t[i] 
   while(m < n)
      + check if x has full rank
	 -> add obs. with smallest t[i] until x has full rank (update m)
      + solve for b in Xb = y
      + compute t[i]
      + select m + 1 obs. with smallest t[i]
   
Implementation
   dgels
      -> if error (i.e. not full rank): add obs. until full rank
   L <- chol(crossprod(x))
   check rank



