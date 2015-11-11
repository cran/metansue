\name{metalm}
\alias{metalm}
\concept{regression}
\title{Fitting Meta-Analytic Linear Models}
\description{Fits meta-analytic linear models. Along with \code{\link[=linearHypothesis.meta.nsue]{linearHypothesis}}, it can be used to carry out meta-regression, meta-comparisons, analyses of variance and covariance, and etcetera.}
\usage{metalm(x, ...)}
\arguments{
  \item{x}{an object of class \code{"meta.nsue"}, or of other classes in other meta-analytic packages.}
  \item{...}{other arguments.}
}
\details{The function \code{metalm} is generic. See \code{\link{metalm.meta.nsue}} for fitting meta-analytic models for \code{"meta.nsue"} objects.}
\seealso{
  \code{\link{metalm.meta.nsue}} for fitting meta-analytic models for \code{"meta.nsue"} objects.

  \code{\link[=meta.nsue]{meta}} for conducting a meta-analysis.

  \code{\link[=linearHypothesis.meta.nsue]{linearHypothesis}} for testing linear hypotheses.

  \code{\link[=metabias.meta.nsue]{metabias}} for testing for funnel plot asymmetry.
}
\keyword{anova}
\keyword{comparison}
\keyword{meta-analysis}
\keyword{regression}
\keyword{test}
