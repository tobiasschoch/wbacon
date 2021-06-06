#===============================================================================
# SUBJECT  Test the implementation of 'wbacon' against the function BACON in the
#          package 'robustX'
# AUTHORS  Tobias Schoch, tobias.schoch@gmail.com, January 23, 2021
# LICENSE  GPL >= 2
# COMMENT  pkg 'robustbase' must be installed
#===============================================================================
library(wbacon)
library(robustX)

#===============================================================================
# Comparison function
#===============================================================================
# this function compares our implementation with the reference implementation
# in pkg 'robustX'
compare <- function(data, name, init = "V2"){
	# we test 11 different values for the parameter 'alpha'
	set_of_alphas <- c(seq(0.3, 0.9, 0.1), 0.95, 0.99, 0.999, 0.9999)
	data <- as.matrix(data)

	res <- NULL
	acc <- sqrt(.Machine$double.eps)
	for (i in 1:length(set_of_alphas)) {

		# reference implementation robustX::BACON with the "V2" (version 2)
		# intialization method of Billor et al. (2000)
		ref.init <- ifelse(init == "V2", "V2", "Mahalanobis")
		reference <- BACON(data, alpha = set_of_alphas[i], init.sel = ref.init,
			verbose = FALSE)

		# our implementation
		new_imp <- wBACON(data, alpha = set_of_alphas[i], version = init)
		deviations <- NULL

		# absolute deviations
		dev_location <- max(abs(reference$center - new_imp$center))
		dev_scatter <- max(abs(reference$cov - new_imp$cov))
		dev_distance <- max(abs(reference$dis - new_imp$dist))

		# deviations in terms of the selected subset
		dev_subset <- sum(reference$subset - new_imp$subset)

		at <- matrix(rep(0, 4), ncol = 4)
		rownames(at) <- paste(set_of_alphas[i])

		# throw an error if the deviations are larger than 'acc'
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

	if (!is.null(res)) {
		cat(name, ": differences detected\n")
		colnames(res) <- c("center", "cov", "distance", "subset")
        print(res)
	}

    return(NROW(res))
}

#===============================================================================
# Tests I
#===============================================================================
# We test the implementations on 7 well known data sets.
errors <- 0

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

# check that version 1.25 (or newer) of robustX is installed
robustX_version <- as.numeric(gsub("-", "", getNamespaceVersion("robustX")))
if (robustX_version >= 1.25) {
    for (i in 1:nrow(setup)) {
        data_name <- setup[i, 1]
        data(list = data_name, package = "robustbase")
        dt <- data.matrix(get(data_name))
        eval(parse(text = paste0("dt <- dt[,", setup[i, 2], "]")))
        errors <- errors + compare(dt, data_name, "V1")
        errors <- errors + compare(dt, data_name, "V2")
    }
    if (errors == 0)
        cat("\nno errors\n\n")
} else {
    cat("Version >= 1.25 of package 'robustX' is required, you have",
        robustX_version, "\n")
}
