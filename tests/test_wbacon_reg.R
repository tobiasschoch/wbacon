#===============================================================================
# SUBJECT  Test the implementation of 'wBACON_reg'
# AUTHORS  Tobias Schoch, tobias.schoch@gmail.com
# LICENSE  GPL >= 2
# COMMENT  pkg 'robustbase' must be installed
#===============================================================================
library(wbacon)
library(robustX)

#===============================================================================
# Tests
#===============================================================================
# We test our implementation against the method robustX::BACON for 5 well known
# data sets.

compare <- function(formula, data, alpha = 0.05, original = TRUE,
    verbose = FALSE)
{
    m <- wBACON_reg(formula, data = data, alpha = alpha, original = original,
        verbose = verbose)
	y <- as.numeric(model.response(m$model))
	x <- model.matrix(m$terms, m$model)
    ref <- suppressWarnings(BACON(x[, -1], y, init.sel = "V2", alpha = alpha,
        verbose = verbose))
    res <- sum(xor(m$subset, unname(ref$subset)))
    call <- match.call()
    if (res > 0)
        cat(paste0("\n The subsets differ in ", res, " places for '",
            call$data, "'\n"))
    res
}

errors <- 0

# check that version 1.25 (or newer) of robustX is installed
robustX_version <- as.numeric(gsub("-", "", getNamespaceVersion("robustX")))
if (robustX_version >= 1.25) {
    data(hbk, package = "robustbase")
    errors <- errors + compare(Y ~ ., hbk)

    data(aircraft, package = "robustbase")
    errors <- errors + compare(Y ~ ., data = aircraft)

    data(education, package = "robustbase")
    errors <- errors + compare(Y ~ Region + X1 + X2 + X2, data = education)

    data(heart, package = "robustbase")
    errors <- errors + compare(clength ~ ., data = heart)

    data(pulpfiber, package = "robustbase")
    errors <- errors + compare(Y1 ~ X1 + X2 + X3, data = pulpfiber)

    if (errors == 0) {
        cat("\nno errors\n\n")
    } else {
        cat("\n")
    }
} else {
    cat("Version >= 1.25 of package 'robustX' is required, you have",
        robustX_version, "\n")
}
