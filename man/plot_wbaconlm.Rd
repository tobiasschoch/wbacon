\name{plot.wbaconlm}
\alias{plot.wbaconlm}
\title{Plot Diagnostics for an Object of Class \code{wbaconlm}}
\usage{
\method{plot}{wbaconlm}(x, which = c(1, 2, 3, 4), hex = FALSE,
	caption = c("Residuals vs Fitted", "Normal Q-Q", "Scale-Location",
		"Standardized Residuals vs Robust Mahalanobis Distance"),
	panel = if (add.smooth) function(x, y, ...)
		panel.smooth(x, y, iter = iter.smooth, ...) else points,
    sub.caption = NULL, main = "",
	ask = prod(par("mfcol")) < length(which) && dev.interactive(),
	...,
	id.n = 3, labels.id = names(residuals(x)), cex.id = 0.75,
	qqline = TRUE,
	add.smooth = getOption("add.smooth"), iter.smooth = 3,
	label.pos = c(4, 2), cex.caption = 1, cex.oma.main = 1.25)
}
\arguments{
	\item{x}{object of class \code{wbaconlm}.}
	\item{which}{if a subset of the plots is required, specify a subset of
		the numbers \code{1:4}, \code{[integer]}.}
	\item{hex}{toogle a hexagonally binned plot, \code{[logical]}, default
        \code{hex = FALSE}.}
	\item{caption}{captions to appear above the plots;
		\code{[character]} vector of valid graphics annotations. It
		can be set to \code{""} or \code{NA} to suppress all captions.}
	\item{panel}{panel function. The useful alternative to
		\code{\link{points}}, \code{\link{panel.smooth}} can be chosen
		by \code{add.smooth = TRUE}.}
	\item{sub.caption}{common title \code{[character]}---above the figures
		if there are more than one; used as \code{sub} (s.\code{\link{title}})
		otherwise. If \code{NULL}, as by default, a possible abbreviated
		version of \code{deparse(x$call)} is used.}
	\item{main}{title to each plot \code{[character]}---in addition
		to \code{caption}.}
	\item{ask}{\code{[logical]}; if \code{TRUE}, the user is \emph{ask}ed before
		each plot, see \code{\link{par}(ask=.)}.}
	\item{\dots}{other parameters to be passed through to plotting
		functions.}
	\item{id.n}{number of points to be labelled in each plot, starting
		with the most extreme, \code{[integer]}.}
	\item{labels.id}{vector of labels \code{[character]}, from which the
		labels for extreme points will be chosen. \code{NULL} uses
		observation numbers.}
	\item{cex.id}{magnification of point labels, \code{[numeric]}.}
	\item{qqline}{\code{[logical]} indicating if a \code{\link{qqline}()}
		should be added to the normal Q-Q plot.}
	\item{add.smooth}{\code{[logical]} indicating if a smoother should be
		added to most plots; see also \code{panel} above.}
	\item{iter.smooth}{the number of robustness iterations \code{[integer]},
		the argument \code{iter} in \code{\link{panel.smooth}()}.}
	\item{label.pos}{positioning of labels \code{[numeric]}, for the left
		half and right half of the graph respectively, for plots 1-3.}
	\item{cex.caption}{controls the size of \code{caption}, \code{[numeric]}.}
	\item{cex.oma.main}{controls the size of the \code{sub.caption} only if
		that is \emph{above} the figures when there is more than one,
		\code{[numeric]}.}
}
\description{
Four plots (selectable by \code{which}) are available for an object of
class \code{wbaconlm} (see \code{\link{wBACON_reg}}): A plot
of residuals against fitted values, a scale-location plot of
\eqn{\sqrt{| residuals |}}{sqrt(| residuals |)} against fitted values,
a Normal Q-Q plot, and a plot of the standardized residuals versus the
robust Mahalanobis distances.
}
\details{
The plots for \code{which \%in\% 1:3} are identical with the
plot method for linear models (see \code{\link[=plot.lm]{plot.lm}}).
There you can find details on the implementation and references.

The standardized residuals vs. robust Mahalanobis distance plot
(\code{which = 4}) has been proposed by Rousseeuw and van Zomeren (1990).
}
\value{
[no return value]
}
\references{
Rousseeuw, P.J. and B.C. van Zomeren (1990).  Unmasking Multivariate Outliers
and Leverage Points,
\emph{Journal of the American Statistical Association} 411, pp. 633-639.
}
\seealso{
\code{\link{wBACON_reg}}
}
