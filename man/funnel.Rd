\name{funnel}
\alias{funnel}
\alias{funnel.meta.nsue}
\concept{plot}
\title{Funnel Plots for \dQuote{meta.nsue} Objects}
\description{Draws a funnel plot.}
\usage{funnel(x, ...)}
\arguments{
  \item{x}{an object of class \code{"meta.nsue"}.}
  \item{...}{other arguments (currently ignored).}
}
\details{The plot shows the residual effect size of the studies on the x-axis against their standard errors on the y-axis. A light gray shadow shows the ellipse approximately containing 95\% of the imputations for studies with Non Statistically-Significant Unreported Effects (NSUEs). A pseudo confidence interval region is drawn with a dashed line.

Asymmetry in the funnel plot may be an indication of publication bias, but this should be taken with caution, especially if sample sizes are too similar, if there are subgroups, or if dichotomous outcomes have been used (some effect estimates are naturally correlated with their standard errors).

Also, publication bias cannot be excluded in case of negative results, whereas it is only one of the possible explanations in case of positive results - others being selective outcome reporting, poor methodological quality in the smaller studies, true heterogeneity (i.e. effect size truly depends on study size), and etcetera.
}
\references{
  Egger, M., Smith, G.D., Schneider, M., Minder, C. (1997)
  Bias in meta-analysis detected by a simple, graphical test.
  \emph{BMJ}, \bold{315}, 629--634, doi:10.1136/bmj.315.7109.629.
}
\author{Joaquim Radua}
\seealso{
  \code{\link{metabias}} for testing for funnel plot asymmetry.

  \code{\link{forest}} for plotting forest plots.

  \code{\link{meta}} for conducting a meta-analysis.

  \code{\link{pdf}} for saving the plot in a PDF file.
}
\examples{
t <- c(3.4, NA, NA, NA, NA, 2.8, 2.1, 3.1, 2.0, 3.4)
n <- c(40, 20, 22, 24, 18, 30, 25, 30, 16, 22)
m <- meta(smc_from_t(t, n))
funnel(m)
}
\keyword{bias}
\keyword{meta-analysis}
\keyword{plot}
