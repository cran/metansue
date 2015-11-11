\name{linearHypothesis}
\alias{linearHypothesis}
\concept{hypothesis}
\title{Test Linear Hypothesis}
\usage{linearHypothesis(x, ...)}
\description{Function for testing a linear hypothesis. Along with \code{\link[=metalm.meta.nsue]{metalm}} (or others depending on the class, e.g. \code{"lm"}), it can be used to carry out meta-regression, meta-comparisons, analyses of variance and covariance, and etcetera.}
\arguments{
  \item{x}{an object of class \code{"meta.nsue"}, or of other classes in other packages (e.g. \code{"glm"}, \code{"lm"}, \code{"lme"}, \code{"mer"}, \code{"merMod"}, \code{"mlm"}, \code{"nlsList"}, \code{"polr"} or \code{"svyglm"}).}
  \item{...}{other arguments.}
}
\details{The function \code{linearHypothesis} is generic. See \code{\link{linearHypothesis.meta.nsue}} for testing linear hypotheses for \code{"meta.nsue"} objects, or other methods from other packages (e.g. \code{\link[car]{linearHypothesis.default}}).}
\seealso{
  \code{\link{linearHypothesis.meta.nsue}} for testing linear hypotheses for \code{"meta.nsue"} objects.

  \code{\link[=metalm.meta.nsue]{metalm}} for fitting meta-analytic models.

  \code{\link[=meta.nsue]{meta}} for conducting a meta-analysis.
}
\keyword{anova}
\keyword{hypothesis}
\keyword{regression}
\keyword{test}
