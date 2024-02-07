plot.wbaconmv <- function(x, which = 1:2,
    caption = c("Robust distance vs. Index",
    "Robust distance vs. Univariate projection"), hex = FALSE, col = 2,
    pch = 19, ask = prod(par("mfcol")) < length(which) && dev.interactive(),
    alpha = 0.05, maxiter = 20, tol = 1e-5, ...)
{
    if (!inherits(x, "wbaconmv"))
        stop("use only with 'wbaconmv' objects")

    if (!is.numeric(which) || any(which < 1) || any(which > 2))
        stop("'which' must be in 1:2")

    show <- rep(FALSE, 2)
    show[which] <- TRUE

    if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }

    if (show[1]) {
        plot(x$dist, xlab = "Index", ylab = "Robust distance",
             main = caption[1], type = "n", ...)
        at <- x$subset == 1
        points(which(at), x$dist[at])
        points(which(!at), x$dist[!at], pch = pch, col = col)
        abline(h = x$cutoff, lty = 2, col = col)
        dev.flush()
    }

    if (show[2]) {
        tmp <- SeparationIndex(x, alpha, tol, maxiter)
        if (tmp$failed) {
            plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "",
                 main = caption[2], ...)
            text(0, 1, labels = "[too few observations]")
        } else {
            if (!tmp$converged)
                warning("Optimal projection: not converged; see arguments 'maxiter'")
            if (hex) {
                requireNamespace("hexbin")
                hb <- hexbin(tmp$proj, x$dist)
                hvp <- hexbin::plot(hb, xlab = "Univariate projection",
                                    ylab = "Robust distance", main = caption[2],
                                    ...)
                hexVP.abline(hvp$plot, h = x$cutoff, lty = 2, col = 2)
            } else {
                if (length(tmp$proj) > 10000)
                    message("Tool tip: A hexbin scatterplot is available (hex = TRUE)\n")
                plot(tmp$proj, x$dist, xlab = "Univariate projection",
                    ylab = "Robust distance", main = caption[2], type = "n",
                    ...)
                at <- x$subset == 1
                points(tmp$proj[at], x$dist[at])
                points(tmp$proj[!at], x$dist[!at], pch = pch, col = col)
                abline(h = x$cutoff, lty = 2, col = col)
            }
        }
        grDevices::dev.flush()
    }
    invisible()
}

# Separation index of Qiu and Joe (2006): Separation index and partial
# membership for clustering, Computational Statistics and Data Analysis 50,
# pp. 585-603
SeparationIndex <- function(object, alpha = 0.05, tol = 1e-5, maxiter = 20)
{
    n_outlier <- object$n - sum(object$subset)
    if (n_outlier < object$p)
        return(list(failed = TRUE))

    stopifnot(alpha > 0, alpha < 1, maxiter > 0, tol > 0)
    C1 <- object$cov
    C2 <- cov(object$x[object$subset == 0, ])
    L1 <- chol(C1)
    L2 <- chol(C2)
    delta <- colMeans(object$x[!object$subset, ]) - object$center

    a0 <- delta / norm(as.matrix(delta), type = "F")
    iter <- 0
    repeat {
        iter <- iter + 1
        p1 <- sqrt(crossprod(L1 %*% a0))[1, 1]
        p2 <- sqrt(crossprod(L2 %*% a0))[1, 1]
        a1 <- (p1 + p2) * solve(C1 / p1 + C2 / p2) %*% delta
        a1 <- a1 / norm(as.matrix(a1), type = "F")
        if (norm(as.matrix(a1 - a0), type = "F") < tol || iter >= maxiter)
            break
        a0 <- a1
    }

    if (crossprod(a1, delta) < 0)
        a1 <- -a1

    d <- crossprod(a1, delta)
    p1 <- sqrt(crossprod(L1 %*% a1))[1, 1]
    p2 <- sqrt(crossprod(L2 %*% a1))[1, 1]
    q_alpha <- qnorm(1 - alpha / 2)
    p12 <- q_alpha * (p1 + p2)
    J <- (d - p12) / (d + p12)

    list(J = J, a = a1, proj = as.vector(object$x %*% a1),
         converged = iter < maxiter, failed = FALSE)
}
