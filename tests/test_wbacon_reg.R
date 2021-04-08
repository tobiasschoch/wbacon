#===============================================================================
# SUBJECT  Test the implementation of 'wBACON_reg'
# AUTHORS  Tobias Schoch, tobias.schoch@gmail.com, February 24, 2021
# LICENSE  GPL >= 2
# COMMENT  pkg 'robustbase' must be installed
#===============================================================================
library(wbacon)

#===============================================================================
# Tests
#===============================================================================
# We test the implementation on 5 well known data sets.
ERR <- 0
NEWLINE <- FALSE

#-------------------------------------------------------------------------------
data(hbk, package = "robustbase")
m <- wBACON_reg(Y ~ ., data = hbk, alpha = 0.95)

ref_beta <- c(-0.1804616286, 0.0813787106, 0.0399018125, -0.0516655770)
tmp <- all.equal(coef(m), ref_beta, check.attributes = FALSE)
if (is.character(tmp)) {
	ERR <- ERR + 1
	cat("hbk: coefficients differ\n")
    print(tmp)
}

ref_is_outlier <- 1L:10L
tmp <- all.equal(which(is_outlier(m)), ref_is_outlier)
if (is.character(tmp)) {
	ERR <- ERR + 1
	cat("hbk: identified outliers differ\n")
}
if (NEWLINE)
    cat("\n")

#-------------------------------------------------------------------------------
data(aircraft, package = "robustbase")
m <- wBACON_reg(Y ~ ., collect = 3, data = aircraft, alpha = 0.95)

ref_beta <- c(14.25805560, -5.02245065, 1.87011709, 0.00197066, -0.00116736)
tmp <- all.equal(coef(m), ref_beta, check.attributes = FALSE)
if (is.character(tmp)) {
	ERR <- ERR + 1
	cat("aircraft: coefficients differ\n")
    print(tmp)
}

ref_is_outlier <- c(1L, 2L, 3L, 4L, 12L, 16L, 19L, 22L)
tmp <- all.equal(which(is_outlier(m)), ref_is_outlier)
if (is.character(tmp)) {
	ERR <- ERR + 1
	cat("aircraft: identified outliers differ\n")
    print(tmp)
}
if (NEWLINE)
    cat("\n")

#-------------------------------------------------------------------------------
data(education, package = "robustbase")
m <- wBACON_reg(Y ~ Region + X1 + X2 + X2, data = education, alpha = 0.95)

ref_beta <- c(7.7854265771, 18.0872441054, 0.0437458151, 0.038913849)
tmp <- all.equal(coef(m), ref_beta, check.attributes = FALSE)
if (is.character(tmp)) {
	ERR <- ERR + 1
	cat("education: coefficients differ\n")
    print(tmp)
}

ref_is_outlier <- c(1L, 3L, 5L, 6L, 7L, 9L, 10L, 13L, 14L, 15L, 16L, 17L, 21L,
    22L, 23L, 24L, 25L, 29L, 30L, 31L, 32L, 36L, 38L, 40L, 42L, 43L, 44L, 45L,
    47L, 49L, 50L)
tmp <- all.equal(which(is_outlier(m)), ref_is_outlier)
if (is.character(tmp)) {
	ERR <- ERR + 1
	cat("education: identified outliers differ\n")
    print(tmp)
}
if (NEWLINE)
    cat("\n")

#-------------------------------------------------------------------------------
data(heart, package = "robustbase")
m <- wBACON_reg(clength ~ ., collect = 3, data = heart, alpha = 0.95)

ref_beta <- c(31.0708746739, -0.1896775299, 0.3425817244)
tmp <- all.equal(coef(m), ref_beta, check.attributes = FALSE)
if (is.character(tmp)) {
	ERR <- ERR + 1
	cat("heart: coefficients differ\n")
    print(tmp)
}

ref_is_outlier <- c(8L, 11L)
tmp <- all.equal(which(is_outlier(m)), ref_is_outlier)
if (is.character(tmp)) {
	ERR <- ERR + 1
	cat("heart: identified outliers differ\n")
    print(tmp)
}
if (NEWLINE)
    cat("\n")

#-------------------------------------------------------------------------------
data(pulpfiber, package = "robustbase")
m <- wBACON_reg(Y1 ~ X1 + X2 + X3, data = pulpfiber, alpha = 0.95,
    collect = 12)
ref_beta <- c(6.2632015426, 2.7453578623, 0.3046841523, 0.1283534872)
tmp <- all.equal(coef(m), ref_beta, check.attributes = FALSE)
if (is.character(tmp)) {
	ERR <- ERR + 1
	cat("pulpfiber: coefficients differ\n")
    print(tmp)
}

ref_is_outlier <- c(22L, 51L, 52L, 56L, 60L, 61L, 62L)
tmp <- all.equal(which(is_outlier(m)), ref_is_outlier)
if (is.character(tmp)) {
	ERR <- ERR + 1
	cat("pulpfiber: identified outliers differ\n")
    print(tmp)
}
if (NEWLINE)
    cat("\n")

#-------------------------------------------------------------------------------
if (ERR == 0)
	cat("\nno errors\n\n")
