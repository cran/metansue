\name{leave1out}
\alias{leave1out}
\alias{leave1out.nsue}
\alias{print.leave1out.nsue}
\alias{summary.leave1out.nsue}
\concept{test}
\title{Compute Leave-One-Out Diagnostics for \dQuote{nsue} Objects}
\description{Repeatedly fits the specified model, leaving out one study at a time.}
\usage{
leave1out(x, ...)

\method{leave1out}{nsue}(x, formula = ~1, hypothesis = NULL,
    n.imp = 500, maxiter = 200, tol = 1e-06, ...)
}
\arguments{
  \item{x}{an object of class \code{"nsue"}.}
  \item{formula}{an object of class \code{"formula"}: a symbolic description of the model to be fitted.}
  \item{hypothesis}{a hypothesis, or \code{NULL} to simply test the coefficients of the model.}
  \item{n.imp}{number of imputations of NSUEs per study.}
  \item{maxiter}{maximum number of iterations in the REML estimation of \eqn{\tau^2}.}
  \item{tol}{tolerance in the REML estimation of \eqn{\tau^2}.}
  \item{...}{other arguments (currently ignored).}
}
\details{Use \code{\link{nsue}}, \code{\link[=nsue]{smc_from_t}}, \code{\link[=nsue]{smd_from_t}} or \code{\link[=nsue]{zcor_from_r}} to create objects of class \code{"nsue"}.

Models for \code{\link{meta}} and \code{leave1out} are specified symbolically. The \code{formula} is a series of terms which specify a linear predictor for \code{x}. A formula specification of the form \code{first + second} indicates a multiple regression by \code{first} and \code{second}.  A specification of the form \code{first:second} indicates the interaction of \code{first} with \code{second}. The specification \code{first*second} is the same as \code{first + second + first:second}.

Each hypothesis must be a matrix (or vector) giving linear combinations of coefficients by rows.}
\value{
 \code{leave1out} returns an object of class \code{"leave1out.nsue"}, which is a list containing a list for each iteration with the component \code{study} (the study discarded) and the component \code{meta}, which is an object of class \code{"meta.nsue"}.

  The functions \code{print} and \code{summary} may be used to print the details or a summary of the results.
}
\references{
  Radua, J., Schmidt, A., Borgwardt, S., Heinz, A., Schlagenhauf, F., McGuire, P., Fusar-Poli, P. (2015)
  Ventral striatal activation during reward processing in psychosis. A neurofunctional meta-analysis.
  \emph{JAMA Psychiatry}, \bold{72}, 1243--51, doi:10.1001/jamapsychiatry.2015.2196.

  Albajes-Eizagirre, A., Solanes, A, Radua, J. (2019)
  Meta-analysis of non-statistically significant unreported effects.
  \emph{Statistical Methods in Medical Research}, \bold{28}, 3741--54, doi:10.1177/0962280218811349.
}
\author{Joaquim Radua}
\seealso{
  \code{\link{nsue}}, \code{\link[=nsue]{smc_from_t}}, \code{\link[=nsue]{smd_from_t}} and \code{\link[=nsue]{zcor_from_r}} for creating objects of class \code{"nsue"}.

  \code{\link{meta}} for conducting a meta-analysis.
}
\examples{
t <- c(3.4, NA, NA, NA, NA, 2.8, 2.1, 3.1, 2.0, 3.4)
n <- c(40, 20, 22, 24, 18, 30, 25, 30, 16, 22)
leave1out(smc_from_t(t, n))
}
\keyword{jackknife}
\keyword{leave1out}
\keyword{meta-analysis}
