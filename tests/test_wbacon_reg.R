#===============================================================================
# SUBJECT  Test the implementation of 'wBACON_reg'
# AUTHORS  Tobias Schoch, tobias.schoch@gmail.com
# LICENSE  GPL >= 2
# COMMENT  pkg 'robustbase' and 'robustX' must be installed
#===============================================================================
library(wbacon)

if (requireNamespace("robustX", quietly = TRUE) &
    requireNamespace("robustbase", quietly = TRUE)) {

    library(robustX)
    library(robustbase)

    #---------------------------------------------------------------------------
    # function to compare wBACON_reg against BACON
    #---------------------------------------------------------------------------
    all_equal <- function(target, current, label,
        tolerance = sqrt(.Machine$double.eps), scale = NULL,
        check.attributes = FALSE)
    {
        if (missing(label))
            stop("Argument 'label' is missing\n")
        res <- all.equal(target, current, tolerance, scale,
            check.attributes = check.attributes)
        if (is.character(res))
            cat(paste0(label, ": ", res, "\n"))
    }
    compare <- function(formula, data, name, alpha = 0.05, original = TRUE,
        verbose = FALSE)
    {
        # our implementation
        m <- wBACON_reg(formula, data = data, alpha = alpha,
                        original = original, verbose = verbose)

        # we extract the response variable and the design matrix
        y <- as.numeric(model.response(m$model))
        x <- model.matrix(m$terms, m$model)

        # reference implementation (robustX)
        ref <- suppressWarnings({
            BACON(x[, -1], y, init.sel = "V2", alpha = alpha,
                  verbose = verbose)
        })

        all_equal(m$subset, ref$subset, name)
    }

    # check that version 1.25 (or newer) of robustX is installed
    robustX_version <-
        as.numeric(gsub("-", "", getNamespaceVersion("robustX")))
    #---------------------------------------------------------------------------
    # Tests
    #---------------------------------------------------------------------------
    # We test our implementation against the method robustX::BACON for 5 well
    # known data sets.
    if (robustX_version >= 1.25) {
        data(hbk, package = "robustbase")
        compare(Y ~ ., hbk, "hbk")

        data(aircraft, package = "robustbase")
        compare(Y ~ ., aircraft, "aircraft")

        data(education, package = "robustbase")
        compare(Y ~ Region + X1 + X2 + X2, education, "education")

        data(heart, package = "robustbase")
        compare(clength ~ ., heart, "heart")

        data(pulpfiber, package = "robustbase")
        compare(Y1 ~ X1 + X2 + X3, pulpfiber, "pulpfiber")
    }
}
