\name{linearHypothesis.meta.nsue}
\alias{linearHypothesis.meta.nsue}
\alias{linearHypothesis.meta.nsue}
\concept{hypothesis}
\title{Test Linear Hypothesis for \dQuote{meta.nsue} Objects}
\description{Function for testing a linear hypothesis. Along with \code{\link[=metalm.meta.nsue]{metalm}}, it can be used to carry out meta-regression, meta-comparisons, analyses of variance and covariance, and etcetera.}
\usage{\method{linearHypothesis}{meta.nsue}(x, hypotheses, ...)}
\arguments{
  \item{x}{an object of class \code{"meta.nsue"}.}
  \item{hypotheses}{a list of hypotheses, or NULL to test the coefficients of the model.}
  \item{...}{other arguments (currently ignored).}
}
\details{Each hypothesis must be a matrix (or vector) giving linear combinations of coefficients by rows.}
\value{\code{linearHypothesis.meta.nsue} returns an object of class \code{"meta.nsue"}.

  The functions \code{print} and \code{summary} may be used to print the details or a summary of the results. The generic accessor functions \code{coefficients}, \code{fitted.values} and \code{residuals} extract various useful features of the value returned by \code{linearHypothesis.meta.nsue}.
}
\references{
  Radua, J., Schmidt, A., Borgwardt, S., Heinz, A., Schlagenhauf, F., McGuire, P., Fusar-Poli, P. (2015)
  Ventral striatal activation during reward processing in psychosis. A neurofunctional meta-analysis.
  \emph{JAMA Psychiatry}, \bold{72}, 1243--51.
}
\author{Joaquim Radua}
\seealso{
  \code{\link[=metalm.meta.nsue]{metalm}} for fitting meta-analytic models.

  \code{\link[=meta.nsue]{meta}} for conducting a meta-analysis.
}
\examples{
t <- c(3.4, NA, NA, NA, NA, 2.8, 2.1, 3.1, 2.0, 3.4)
n <- c(40, 20, 22, 24, 18, 30, 25, 30, 16, 22)
m <- meta(smc_from_t(t, n))
mean.age <- rnorm(10, 30, 5)
m <- metalm(m, ~ mean.age)
linearHypothesis(m, list(age_10years = c(0, 10)))
}
\keyword{anova}
\keyword{hypothesis}
\keyword{regression}
\keyword{test}
