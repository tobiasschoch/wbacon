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

#-------------------------------------------------------------------------------
data(hbk, package = "robustbase")
m <- wBACON_reg(Y ~ ., data = hbk, alpha = 0.95)

ref_beta <- c(-0.18046162865083981, 0.08137871068818925, 0.03990181252317023,
	-0.05166557707657377)
if (!all.equal(coef(m), ref_beta, check.attributes = FALSE)) {
	ERR <- ERR + 1
	cat("hbk: coefficients differ\n")
}

ref_is_outlier <- 1:10
if (!all(which(is_outlier(m)) == ref_is_outlier)) {
	ERR <- ERR + 1
	cat("hbk: identified outliers differ\n")
}

#-------------------------------------------------------------------------------
data(aircraft, package = "robustbase")
m <- wBACON_reg(Y ~ ., collect = 3, data = aircraft, alpha = 0.95)

ref_beta <- c(9.50074034620055308, -3.04879688663568649, 1.21003298781689805,
	0.00138096419402889, -0.00055485757094962)
if (!all.equal(coef(m), ref_beta, check.attributes = FALSE)) {
	ERR <- ERR + 1
	cat("aircraft: coefficients differ\n")
}

ref_is_outlier <- c(16, 22)
if (!all(which(is_outlier(m)) == ref_is_outlier)) {
	ERR <- ERR + 1
	cat("aircraft: identified outliers differ\n")
}

#-------------------------------------------------------------------------------
data(education, package = "robustbase")
m <- wBACON_reg(Y ~ Region + X1 + X2 + X2, data = education, alpha = 0.95)

ref_beta <- c(20.8925307788478278, 9.6910657876150736, 0.0614905565468016,
	0.0408586963725780)
if (!all.equal(coef(m), ref_beta, check.attributes = FALSE)) {
	ERR <- ERR + 1
	cat("education: coefficients differ\n")
}

ref_is_outlier <- c(15, 50)
if (!all(which(is_outlier(m)) == ref_is_outlier)) {
	ERR <- ERR + 1
	cat("education: identified outliers differ\n")
}

#-------------------------------------------------------------------------------
data(heart, package = "robustbase")
m <- wBACON_reg(clength ~ ., collect = 3, data = heart, alpha = 0.95)

ref_beta <- c(31.0708746739284010, -0.1896775299044007, 0.3425817244946800)
if (!all.equal(coef(m), ref_beta, check.attributes = FALSE)) {
	ERR <- ERR + 1
	cat("heart: coefficients differ\n")
}

ref_is_outlier <- c(8, 11)
if (!all(which(is_outlier(m)) == ref_is_outlier)) {
	ERR <- ERR + 1
	cat("heart: identified outliers differ\n")
}

#-------------------------------------------------------------------------------
data(pulpfiber, package = "robustbase")
m <- wBACON_reg(Y1 ~ X1 + X2 + X3, data = pulpfiber, alpha = 0.95)

ref_beta <- c(5.1075264612314788, 2.1740952964309233, 0.3177965854148064,
	0.1433700922951207)
if (!all.equal(coef(m), ref_beta, check.attributes = FALSE)) {
	ERR <- ERR + 1
	cat("pulpfiber: coefficients differ\n")
}

ref_is_outlier <- c(56, 60, 61, 62)
if (!all(which(is_outlier(m)) == ref_is_outlier)) {
	ERR <- ERR + 1
	cat("pulpfiber: identified outliers differ\n")
}
#-------------------------------------------------------------------------------
if (ERR == 0)
	cat("\nno errors\n")

