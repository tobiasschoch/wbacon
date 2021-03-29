plot.wbaconlm <- function(x, which = c(1, 2, 3, 4),
    hex = FALSE, caption = c("Residuals vs Fitted", "Normal Q-Q",
    "Scale-Location", "Standardized Residuals vs Robust Mahalanobis Distance"),
	panel = if (add.smooth) function(x, y, ...) panel.smooth(x, y,
    iter = iter.smooth, ...) else points, sub.caption = NULL, main = "",
	ask = prod(par("mfcol")) < length(which) && dev.interactive(), ...,
	id.n = 3, labels.id = names(residuals(x)), cex.id = 0.75, qqline = TRUE,
	add.smooth = getOption("add.smooth"), iter.smooth = 3,
	label.pos = c(4, 2), cex.caption = 1, cex.oma.main = 1.25)
{
    dropInf <- function(x, h)
	{
		if (any(isInf <- h >= 1)) {
			warning(gettextf("not plotting observations with leverage one:\n %s",
                paste(which(isInf), collapse = ", ")), call. = FALSE,
				domain = NA)
            x[isInf] <- NaN
        }
        x
    }

    if (!is.numeric(which) || any(which < 1) || any(which > 6))
        stop("'which' must be in 1:4")

    show <- rep(FALSE, 6)
    show[which] <- TRUE

    subset0 <- x$subset
    r <- x$residuals[subset0]
    yh <- x$fitted.values[subset0]
    w <- x$weights[subset0]
    n <- length(r)

	# regression scale
	s <- sqrt(sum(w * r^2) / x$df.residual)

    if (any(show[c(2L:4L)])) {
        ylab5 <- ylab23 <- "Standardized residuals"
		xmat <- stats::model.matrix(x$terms, x$model)[subset0, ]
		Q <- qr.Q(qr(sqrt(w) * xmat))
		hii <- rowSums(Q^2) * w
        rs <- dropInf(r / (s * sqrt(1 - hii)), hii)
	}
    if (any(show[c(1L, 3L)]))
        l.fit <- "Fitted values"
    if (is.null(id.n)) {
        id.n <- 0
    } else {
        id.n <- as.integer(id.n)
        if (id.n < 0L || id.n > n)
            stop(gettextf("'id.n' must be in {1,..,%d}", n), domain = NA)
    }
    if (id.n > 0L) {
        if (is.null(labels.id))
            labels.id <- paste(1L:n)
        iid <- 1L:id.n
        show.r <- sort.list(abs(r), decreasing = TRUE)[iid]
        if (any(show[2L:3L]))
            show.rs <- sort.list(abs(rs), decreasing = TRUE)[iid]
        text.id <- function(x, y, ind, adj.x = TRUE)
		{
            labpos <- if (adj.x)
                label.pos[1 + as.numeric(x > mean(range(x)))]
            else
				3
            text(x, y, labels.id[ind], cex = cex.id, xpd = TRUE,
                pos = labpos, offset = 0.25)
        }
    }
    getCaption <- function(k)
	{
		if (length(caption) < k)
        	NA_character_
	    else
			as.graphicsAnnot(caption[[k]])
	}
    if (is.null(sub.caption)) {
        cal <- x$call
        if (!is.na(m.f <- match("formula", names(cal)))) {
            cal <- cal[c(1, m.f)]
            names(cal)[2L] <- ""
        }
        cc <- deparse(cal, 80)
        nc <- nchar(cc[1L], "c")
        abbr <- length(cc) > 1 || nc > 75
        sub.caption <- if (abbr)
            paste(substr(cc[1L], 1L, min(75L, nc)), "...")
        else
			cc[1L]
    }
    one.fig <- prod(par("mfcol")) == 1
    if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }
	# Residuals vs Fitted
    if (show[1L]) {
        ylim <- range(r, na.rm = TRUE)
        if (id.n > 0)
            ylim <- extendrange(r = ylim, f = 0.08)
        dev.hold()
        if (hex) {
            requireNamespace("hexbin")
            hb <- hexbin(yh, r, ybnds = ylim)
            hvp <- hexbin::plot(hb, xlab = l.fit, ylab = "Residuals",
                main = main)
            hexVP.abline(hvp$plot, h = 0, lty = 3, col = "gray")
            hexVP.loess(hb, hvp = hvp$plot, span = 2 / 3, ...)
        } else {
            plot(yh, r, xlab = l.fit, ylab = "Residuals", main = main,
                ylim = ylim, type = "n", ...)
            panel(yh, r, ...)
            if (one.fig)
                title(sub = sub.caption, ...)
            mtext(getCaption(1), 3, 0.25, cex = cex.caption)
            if (id.n > 0) {
                y.id <- r[show.r]
                y.id[y.id < 0] <- y.id[y.id < 0] - strheight(" ") / 3
                text.id(yh[show.r], y.id, show.r)
            }
            abline(h = 0, lty = 3, col = "gray")
        }
        dev.flush()
    }
	# Normal Q-Q
    if (show[2L]) {
        ylim <- range(rs, na.rm = TRUE)
        ylim[2L] <- ylim[2L] + diff(ylim) * 0.075
        dev.hold()
        qq <- qqnorm(rs, main = main, ylab = ylab23, ylim = ylim, ...)
        if (qqline)
            qqline(rs, lty = 3, col = "gray50")
        if (one.fig)
            title(sub = sub.caption, ...)
        mtext(getCaption(2), 3, 0.25, cex = cex.caption)
        if (id.n > 0)
            text.id(qq$x[show.rs], qq$y[show.rs], show.rs)
        dev.flush()
    }
	# Scale-Location
    if (show[3L]) {
        sqrtabsr <- sqrt(abs(rs))
        ylim <- c(0, max(sqrtabsr, na.rm = TRUE))
        yl <- as.expression(substitute(sqrt(abs(YL)),
            list(YL = as.name(ylab23))))
        yhn0 <- yh
        dev.hold()
        if (hex) {
            hb <- hexbin(yhn0, sqrtabsr, ybnds = ylim)
            hvp <- plot(hb, xlab = l.fit, ylab = yl, main = main)
            hexVP.loess(hb, hvp = hvp$plot, span = 2 / 3, ...)
        } else {
            plot(yhn0, sqrtabsr, xlab = l.fit, ylab = yl, main = main,
                ylim = ylim, type = "n", ...)
            panel(yhn0, sqrtabsr, ...)
            if (one.fig)
                title(sub = sub.caption, ...)
            mtext(getCaption(3), 3, 0.25, cex = cex.caption)
            if (id.n > 0)
                text.id(yhn0[show.rs], sqrtabsr[show.rs], show.rs)
        }
        dev.flush()
    }
	# Standardized residuals vs. robust Mahalanobis distances
    if (show[4L]) {
		xlim <- c(0, max(x$mv$dist, na.rm = TRUE))
		dev.hold()
        if (hex) {
            hb <- hexbin(x$mv$dist, x$residuals/s, xbnds = xlim)
            hvp <- plot(hb, xlab = "Robust distance",
                ylab = "Standardized residuals", main = main)
            hexVP.abline(hvp$plot, h = 0, lty = 3, col = "gray")
            hexVP.abline(hvp$plot, v = x$mv$cutoff, h = c(-x$reg$cutoff,
                x$reg$cutoff), lty = 2, col = 2)
        } else {
            plot(x$mv$dist, x$residuals / s, xlab = "Robust distance",
                ylab = "Standardized residuals", main = main, xlim = xlim,
                type = "n", ...)
            points(x$mv$dist[subset0], x$residuals[subset0] / s, ...)
            points(x$mv$dist[!subset0], x$residuals[!subset0] / s, pch = 19,
                col = 2, ...)
            abline(h = 0, lty = 3, col = "gray")
            abline(v = x$mv$cutoff, h = c(-x$reg$cutoff, x$reg$cutoff),
                lty = 2, col = 2)
           if (one.fig)
                title(sub = sub.caption, ...)
            mtext(getCaption(4), 3, 0.25, cex = cex.caption)
        }
		dev.flush()
    }
    if (!one.fig && par("oma")[3L] >= 1)
        mtext(sub.caption, outer = TRUE, cex = cex.oma.main)
    invisible()
}
