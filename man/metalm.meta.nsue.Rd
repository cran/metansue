\name{metalm.meta.nsue}
\alias{metalm.meta.nsue}
\concept{regression}
\title{Fitting Meta-Analytic Linear Models for \dQuote{meta.nsue} Objects}
\description{Fits meta-analytic linear models. Along with \code{\link[=linearHypothesis.meta.nsue]{linearHypothesis}}, it can be used to carry out meta-regression, meta-comparisons, analyses of variance and covariance, and etcetera.}
\usage{\method{metalm}{meta.nsue}(x, formula, maxiter = 100, tol = 1e-06, ...)}
\arguments{
  \item{x}{an object of class \code{"meta.nsue"}.}
  \item{formula}{an object of class \code{"formula"}: a symbolic description of the model to be fitted.}
  \item{maxiter}{maximum number of iterations in the REML estimation of \eqn{\tau^2}.}
  \item{tol}{tolerance in the REML estimation of \eqn{\tau^2}.}
  \item{...}{other arguments (currently ignored).}
}
\details{Models for \code{\link{meta.nsue}}, \code{\link{leave1out.nsue}} and \code{metalm.meta.nsue} are specified symbolically. The \code{formula} is a series of terms which specify a linear predictor for \code{x}. A formula specification of the form \code{first + second} indicates a multiple regression by \code{first} and \code{second}.  A specification of the form \code{first:second} indicates the interaction of \code{first} with \code{second}. The specification \code{first*second} is the same as \code{first + second + first:second}.}
\value{\code{metalm.meta.nsue} returns an object of class \code{"meta.nsue"}.

  The functions \code{print} and \code{summary} may be used to print the details or a summary of the results. The generic accessor functions \code{coefficients}, \code{fitted.values} and \code{residuals} extract various useful features of the value returned by \code{metalm.meta.nsue}.
}
\references{
  Radua, J., Schmidt, A., Borgwardt, S., Heinz, A., Schlagenhauf, F., McGuire, P., Fusar-Poli, P. (2015)
  Ventral striatal activation during reward processing in psychosis. A neurofunctional meta-analysis.
  \emph{JAMA Psychiatry}, \bold{72}, 1243--1251.
}
\author{Joaquim Radua}
\seealso{
  \code{\link[=meta.nsue]{meta}} for conducting a meta-analysis.

  \code{\link[=linearHypothesis.meta.nsue]{linearHypothesis}} for testing linear hypotheses.

  \code{\link[=metabias.meta.nsue]{metabias}} for testing for funnel plot asymmetry.
}
\examples{
t <- c(3.4, NA, NA, NA, NA, 2.8, 2.1, 3.1, 2.0, 3.4)
n <- c(40, 20, 22, 24, 18, 30, 25, 30, 16, 22)
m <- meta(smc_from_t(t, n))
mean.age <- rnorm(10, 30, 5)
metalm(m, ~ mean.age)
}
\keyword{anova}
\keyword{comparison}
\keyword{meta-analysis}
\keyword{regression}
\keyword{test}
