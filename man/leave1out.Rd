\name{leave1out}
\alias{leave1out}
\concept{test}
\title{Compute Leave-One-Out Diagnostics}
\description{Repeatedly fits the specified model, leaving out one study at a time.}
\usage{leave1out(x, ...)}
\arguments{
  \item{x}{an object of class \code{"nsue"}, or of other classes in other meta-analytic packages (e.g. \code{"rma"}).}
  \item{...}{other arguments.}
}
\details{The function \code{leave1out} is generic. See \code{\link{leave1out.nsue}} for computing leave-one-out diagnostics for \code{"nsue"} objects, or other methods from other meta-analytic packages (e.g. \code{\link[metafor]{leave1out.rma.uni}}).}
\seealso{
  \code{\link{leave1out.nsue}} for computing leave-one-out diagnostics for \code{"nsue"} objects.

  \code{\link{smc_from_t}}, \code{\link[=smc_from_t]{smd_from_t}}, \code{\link[=smc_from_t]{z_from_r}} and \code{\link[=smc_from_t]{r_in_smd_from_t_means_and_sds1}} for creating \code{"nsue"} objects.

  \code{\link[=meta.nsue]{meta}} for conducting a meta-analysis.
}
\keyword{jackknife}
\keyword{leave1out}
\keyword{meta-analysis}
