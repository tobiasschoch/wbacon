#===============================================================================
# SUBJECT  Test the implementation of 'wbacon' against the function BACON in the
#          package 'robustX'
# AUTHORS  Tobias Schoch, tobias.schoch@gmail.com, Feburary 7, 2024
# LICENSE  GPL >= 2
# COMMENT  pkg 'robustX' is required
#===============================================================================
library(wbacon)

# pkg robustX is required
if (requireNamespace("robustX", quietly = TRUE)) {
    library(robustX)

    # utility function
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

    #---------------------------------------------------------------------------
    # Comparison function
    #---------------------------------------------------------------------------
    # this function compares our implementation with the reference
    # implementation in pkg 'robustX'
    compare <- function(data, name, init = "V2"){
        # we test 11 different values for the parameter 'alpha'
        set_of_alphas <- c(seq(0.3, 0.9, 0.1), 0.95, 0.99, 0.999, 0.9999)
        data <- as.matrix(data)

        for (i in 1:length(set_of_alphas)) {
            # reference implementation robustX::BACON with the "V2" (version 2)
            # intialization method of Billor et al. (2000)
            ref.init <- ifelse(init == "V2", "V2", "Mahalanobis")
            reference <- suppressWarnings({
                mvBACON(data, alpha = set_of_alphas[i], init.sel = ref.init,
                        verbose = FALSE)
            })

            # our implementation
            new_imp <- wBACON(data, alpha = set_of_alphas[i], version = init)

            # comparison
            meth <- paste0("data = ", name, "init" = init, ", ",
                           "alpha = ", set_of_alphas[i])
            all_equal(reference$center, new_imp$center, paste0("center: ", meth))
            all_equal(reference$cov, new_imp$cov, paste0("cov: ", meth))
            all_equal(reference$dis, new_imp$dist, paste0("dist: ", meth))
            all_equal(reference$subset, new_imp$subset == 1,
                      paste0("subset: ", meth))
        }
    }

    #---------------------------------------------------------------------------
    # Tests
    #---------------------------------------------------------------------------
    # We test the implementations on 7 well known data sets.
    setup <- matrix(c(
    #------------------------------
    #    DATASET        VARIABLES
    #------------------------------
    "hbk",          "1:3",
    "bushfire",     "1:5",
    "aircraft",     "1:4",
    "education",    "2:4",
    "heart",        "1:2",
    "milk",         "1:8",
    "pulpfiber",    "1:8"), byrow = TRUE, ncol = 2)

    robustX_version <-
        as.numeric(gsub("-", "", getNamespaceVersion("robustX")))

    # check that version 1.25 (or newer) of robustX is installed
    if (robustX_version >= 1.25) {
        for (i in 1:nrow(setup)) {
            data_name <- setup[i, 1]
            data(list = data_name, package = "robustbase")
            dt <- data.matrix(get(data_name))
            eval(parse(text = paste0("dt <- dt[,", setup[i, 2], "]")))
            compare(dt, data_name, "V1")
            compare(dt, data_name, "V2")
        }
    }
}
