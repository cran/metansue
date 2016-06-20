\name{leave1out.nsue}
\alias{leave1out.nsue}
\alias{print.leave1out.nsue}
\alias{summary.leave1out.nsue}
\concept{test}
\title{Compute Leave-One-Out Diagnostics for \dQuote{meta.nsue} Objects}
\description{Repeatedly fits the specified model, leaving out one study at a time.}
\usage{\method{leave1out}{nsue}(x, data = data.frame(), formula = ~1, hypotheses = NULL,
    n.imp = 50, n.bins = 200, maxiter = 200, tol = 1e-06, ...)}
\arguments{
  \item{x}{an object of class \code{"nsue"}.}
  \item{data}{an optional data frame containing variables to be used as regressors in the maximum likelihood step.}
  \item{formula}{an object of class \code{"formula"}: a symbolic description of the model to be fitted.}
  \item{hypotheses}{a list of hypotheses, or \code{NULL} to simply test the coefficients of the model.}
  \item{n.imp}{number of imputations of NSUEs per study.}
  \item{n.bins}{number of bins used in the imputations.}
  \item{maxiter}{maximum number of iterations in the REML estimation of \eqn{\tau^2}.}
  \item{tol}{tolerance in the REML estimation of \eqn{\tau^2}.}
  \item{...}{other arguments (currently ignored).}
}
\details{Use \code{\link{smc_from_t}}, \code{\link[=smc_from_t]{smd_from_t}}, \code{\link[=smc_from_t]{z_from_r}} or \code{\link[=smc_from_t]{r_in_smd_from_t_means_and_sds1}} to create \code{"nsue"} objects.

Models for \code{\link{meta.nsue}}, \code{leave1out.nsue} and \code{\link{metalm.meta.nsue}} are specified symbolically. The \code{formula} is a series of terms which specify a linear predictor for \code{x}. A formula specification of the form \code{first + second} indicates a multiple regression by \code{first} and \code{second}.  A specification of the form \code{first:second} indicates the interaction of \code{first} with \code{second}. The specification \code{first*second} is the same as \code{first + second + first:second}.

Each hypothesis must be a matrix (or vector) giving linear combinations of coefficients by rows.}
\value{
 \code{leave1out.nsue} returns an object of class \code{"leave1out.nsue"}, which is a list containing a list for each iteration with the component \code{study} (the study discarded) and the component \code{meta.nsue}, which is an object of class \code{"meta.nsue"}.
 
  The functions \code{print} and \code{summary} may be used to print the details or a summary of the results.
}
\references{
  Radua, J., Schmidt, A., Borgwardt, S., Heinz, A., Schlagenhauf, F., McGuire, P., Fusar-Poli, P. (2015)
  Ventral striatal activation during reward processing in psychosis. A neurofunctional meta-analysis.
  \emph{JAMA Psychiatry}, \bold{72}, 1243--51.
}
\author{Joaquim Radua}
\seealso{
  \code{\link{smc_from_t}}, \code{\link[=smc_from_t]{smd_from_t}}, \code{\link[=smc_from_t]{z_from_r}} and \code{\link[=smc_from_t]{r_in_smd_from_t_means_and_sds1}} for creating \code{"nsue"} objects.

  \code{\link[=meta.nsue]{meta}} for conducting a meta-analysis.
}
\examples{
t <- c(3.4, NA, NA, NA, NA, 2.8, 2.1, 3.1, 2.0, 3.4)
n <- c(40, 20, 22, 24, 18, 30, 25, 30, 16, 22)
leave1out(smc_from_t(t, n))
}
\keyword{jackknife}
\keyword{leave1out}
\keyword{meta-analysis}
