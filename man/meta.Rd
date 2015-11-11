\name{meta}
\alias{meta}
\concept{meta-analysis}
\title{Meta-Analysis}
\description{Conduct a meta-analysis.}
\usage{meta(x, ...)}
\arguments{
  \item{x}{an object of class \code{"nsue"}, or of other classes in other meta-analytic packages.}
  \item{...}{other arguments.}
}
\details{The function \code{meta} is generic. See \code{\link{meta.nsue}} for conducting a meta-analysis for \code{"nsue"} objects.}
\seealso{
  \code{\link{meta.nsue}} for conducting a meta-analysis for \code{"nsue"} objects.

  \code{\link{smc_from_t}}, \code{\link[=smc_from_t]{smd_from_t}} and \code{\link[=smc_from_t]{z_from_r}} for creating \code{"nsue"} objects.

  \code{\link[=forest.meta.nsue]{forest}} for plotting forest plots.

  \code{\link[=funnel.meta.nsue]{funnel}} for plotting funnel plots.

  \code{\link[=metabias.meta.nsue]{metabias}} for testing for funnel plot asymmetry.

  \code{\link[=leave1out.nsue]{leave1out}} for computing leave-one-out diagnostics.

  \code{\link[=metalm.meta.nsue]{metalm}} for fitting meta-analytic models.

  \code{\link[=linearHypothesis.meta.nsue]{linearHypothesis}} for testing linear hypotheses.
}
\keyword{meta-analysis}
