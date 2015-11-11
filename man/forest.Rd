\name{forest}
\alias{forest}
\concept{plot}
\title{Forest Plots}
\description{Draws a forest plot.}
\usage{forest(x, ...)}
\arguments{
  \item{x}{an object of class \code{"meta.nsue"}, or of other classes in other meta-analytic packages (e.g. \code{"meta"} or \code{"rma"}).}
  \item{...}{other arguments.}
}
\details{The function \code{forest} is generic. See \code{\link{forest.meta.nsue}} for plotting forest plots for \code{"meta.nsue"} objects, or other methods from other meta-analytic packages (e.g. \code{\link[meta]{forest.meta}} or \code{\link[metafor]{forest.rma}}).}
\seealso{
  \code{\link{forest.meta.nsue}} for plotting forest plots for \code{"meta.nsue"} objects.

  \code{\link[=funnel.meta.nsue]{funnel}} for plotting funnel plots.

  \code{\link[=meta.nsue]{meta}} for conducting a meta-analysis.

  \code{\link[=linearHypothesis.meta.nsue]{linearHypothesis}} for changing the hypotheses.

  \code{\link{pdf}} for saving the plot in a PDF file.
}
\keyword{meta-analysis}
\keyword{plot}
