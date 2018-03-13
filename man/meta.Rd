\name{meta}
\alias{meta}
\alias{meta.nsue}
\alias{coef.meta.nsue}
\alias{fitted.meta.nsue}
\alias{print.meta.nsue}
\alias{residuals.meta.nsue}
\alias{summary.meta.nsue}
\concept{meta-analysis}
\title{Meta-Analysis of Studies with Non Statistically-Significant Unreported Effects}
\description{Conduct a meta-analysis. MetaNSUE is a meta-analytic method that allows an unbiased inclusion of studies with Non Statistically-Significant Unreported Effects (NSUEs).}
\usage{
meta(x, ...)

\method{meta}{nsue}(x, formula = ~1, hypothesis = NULL,
    n.imp = 500, n.bins = 200, maxiter = 200, tol = 1e-06, ...)
}
\arguments{
  \item{x}{an object of class \code{"nsue"}.}
  \item{formula}{an object of class \code{"formula"}: a symbolic description of the model to be fitted.}
  \item{hypothesis}{a hypothesis, or NULL to test the main coefficient of the model.}
  \item{n.imp}{number of imputations of NSUEs.}
  \item{n.bins}{number of bins used in the imputations.}
  \item{maxiter}{maximum number of iterations in the REML estimation of \eqn{\tau^2}.}
  \item{tol}{tolerance in the REML estimation of \eqn{\tau^2}.}
  \item{...}{other arguments (currently ignored).}
}
\details{Use \code{\link{smc_from_t}}, \code{\link[=smc_from_t]{smd_from_t}}, \code{\link[=smc_from_t]{z_from_r}} or \code{\link[=smc_from_t]{r_in_smd_from_t_means_and_sds1}} to create \code{"nsue"} objects.

Models for \code{meta} and \code{\link{leave1out}} are specified symbolically. The \code{formula} is a series of terms which specify a linear predictor for \code{x}. A formula specification of the form \code{first + second} indicates a multiple regression by \code{first} and \code{second}.  A specification of the form \code{first:second} indicates the interaction of \code{first} with \code{second}. The specification \code{first*second} is the same as \code{first + second + first:second}.

Each hypothesis must be a matrix (or vector) giving linear combinations of coefficients by rows.}
\value{
 \code{meta} returns an object of class \code{"meta.nsue"}, which is a list containing the following components:

  \item{measure}{the effect-size measure used.}
  \item{known}{a list with the known effect sizes and their indexs.}
  \item{unknown}{a list with the imputations of NSUEs and their indexs.}
  \item{y.var}{the variances if the effect sizes.}
  \item{y2var_k1}{a constant needed to derive the variances.}
  \item{y2var_k2}{a constant needed to derive the variances.}
  \item{labels}{the labels of the studies.}
  \item{rm}{a list with the expected correlation between repeated-measures studies, a conversion matrix and the study weights.}
  \item{heterogeneity}{a list with \eqn{I^2} and Q test.}
  \item{model}{a list with the formula, matrix and coefficients of the model.}
  \item{hypothesis}{the matrixs and coefficients of the hypothesis.}

  The functions \code{print} and \code{summary} may be used to print the details or a summary of the results. The generic accessor functions \code{coefficients}, \code{fitted.values} and \code{residuals} extract various useful features of the value returned by \code{meta}.
}
\references{
  Radua, J., Schmidt, A., Borgwardt, S., Heinz, A., Schlagenhauf, F., McGuire, P., Fusar-Poli, P. (2015)
  Ventral striatal activation during reward processing in psychosis. A neurofunctional meta-analysis.
  \emph{JAMA Psychiatry}, \bold{72}, 1243--51.
}
\author{Joaquim Radua}
\seealso{
  \code{\link{smc_from_t}}, \code{\link[=smc_from_t]{smd_from_t}}, \code{\link[=smc_from_t]{z_from_r}} and \code{\link[=smc_from_t]{r_in_smd_from_t_means_and_sds1}} for creating \code{"nsue"} objects.

  \code{\link{forest}} for plotting forest plots.

  \code{\link{funnel}} for plotting funnel plots.

  \code{\link{metabias}} for testing for funnel plot asymmetry.

  \code{\link{leave1out}} for computing leave-one-out diagnostics.
}
\examples{
t <- c(3.4, NA, NA, NA, NA, 2.8, 2.1, 3.1, 2.0, 3.4)
n <- c(40, 20, 22, 24, 18, 30, 25, 30, 16, 22)
meta(smc_from_t(t, n))
}
\keyword{meta-analysis}
