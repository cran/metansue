\name{funnel}
\alias{funnel}
\concept{plot}
\title{Funnel Plots}
\description{Draws a funnel plot.}
\usage{funnel(x, ...)}
\arguments{
  \item{x}{an object of class \code{"meta.nsue"}, or of other classes in other meta-analytic packages (e.g. \code{"meta"} or \code{"rma"}).}
  \item{...}{other arguments.}
}
\details{The function \code{funnel} is generic. See \code{\link{funnel.meta.nsue}} for plotting funnel plots for \code{"meta.nsue"} objects, or other methods from other meta-analytic packages (e.g. \code{\link[meta]{funnel.meta}} or \code{\link[metafor]{funnel.rma}}).}
\seealso{
  \code{\link{funnel.meta.nsue}} for plotting funnel plots for \code{"meta.nsue"} objects.

  \code{\link[=metabias.meta.nsue]{metabias}} for testing for funnel plot asymmetry.

  \code{\link[=forest.meta.nsue]{forest}} for plotting forest plots.

  \code{\link[=meta.nsue]{meta}} for conducting a meta-analysis.

  \code{\link{pdf}} for saving the plot in a PDF file.
}
\keyword{bias}
\keyword{meta-analysis}
\keyword{plot}
