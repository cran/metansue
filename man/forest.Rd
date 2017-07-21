\name{forest}
\alias{forest}
\alias{forest.meta.nsue}
\alias{plot.meta.nsue}
\concept{plot}
\title{Forest Plots for \dQuote{meta.nsue} Objects}
\description{Draws a forest plot.}
\usage{forest(x, ...)}
\arguments{
  \item{x}{an object of class \code{"meta.nsue"}.}
  \item{...}{other arguments (currently ignored).}
}
\details{The plot shows the observed or imputed effect sizes of the studies with the corresponding 95\% confidence intervals. A light gray shadow shows the interval containing 95\% of the imputations for studies with Non Statistically-Significant Unreported Effects (NSUEs).

The effect size of the hypothesis of the model with the corresponding 95\% confidence intervals is added to the bottom of the forest plot.
}
\value{The function invisibly returns the optimal width and height of the plot, useful when calling the function \code{\link{pdf}}.}
\author{Joaquim Radua}
\seealso{
  \code{\link{funnel}} for plotting funnel plots.

  \code{\link{meta}} for conducting a meta-analysis.

  \code{\link{pdf}} for saving the plot in a PDF file.
}
\examples{
t <- c(3.4, NA, NA, NA, NA, 2.8, 2.1, 3.1, 2.0, 3.4)
n <- c(40, 20, 22, 24, 18, 30, 25, 30, 16, 22)
m <- meta(smc_from_t(t, n))
forest(m)
}
\keyword{meta-analysis}
\keyword{plot}
