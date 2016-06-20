\name{metabias.meta.nsue}
\alias{metabias.meta.nsue}
\concept{test}
\title{Test for Funnel Plot Asymmetry for \dQuote{meta.nsue} Objects}
\description{Test for funnel plot asymmetry.}
\usage{\method{metabias}{meta.nsue}(x, maxiter = 100, tol = 1e-06, ...)}
\arguments{
  \item{x}{an object of class \code{"meta.nsue"}.}
  \item{maxiter}{maximum number of iterations in the REML estimation of \eqn{\tau^2}.}
  \item{tol}{tolerance in the REML estimation of \eqn{\tau^2}.}
  \item{...}{other arguments (currently ignored).}
}
\details{The test, based on a meta-regression by the standard error, aims to detect asymmetry in the funnel plot, which may be an indication of publication bias.

However, results should be taken with caution, especially if there are too few studies (at least 10 studies were suggested by Sterne et al. (2011)), if sample sizes are too similar, if there are outliers or influential studies or subgroups in the plot, or if dichotomous outcomes have been used (some effect estimates are naturally correlated with their standard errors).

Also, publication bias cannot be excluded in case of negative results, whereas it is only one of the possible explanations in case of positive results - others being selective outcome reporting, poor methodological quality in the smaller studies, true heterogeneity (i.e. effect size truly depends on study size), and etcetera (Egger 1997).}
\value{A list with class \code{"htest"} containing the following components:
  \item{data.name}{a character string giving the name of the data.}
  \item{method}{a character string indicating that a \code{"meta.nsue"} regression test for funnel plot asymmetry was performed.}
  \item{statistic}{the value of the z-statistic.}
  \item{p.value}{the p-value for the test.}
}
\references{
  Egger, M., Smith, G.D., Schneider, M., Minder, C. (1997)
  Bias in meta-analysis detected by a simple, graphical test.
  \emph{BMJ}, \bold{315}, 629--634, doi:10.1136/bmj.315.7109.629.

  Radua, J., Schmidt, A., Borgwardt, S., Heinz, A., Schlagenhauf, F., McGuire, P., Fusar-Poli, P. (2015)
  Ventral striatal activation during reward processing in psychosis. A neurofunctional meta-analysis.
  \emph{JAMA Psychiatry}, \bold{72}, 1243--51.

  Sterne, J.A., Sutton, A.J., Ioannidis, J.P., Terrin, N., Jones, D.R., Lau, J., Carpenter, J., 
Rucker, G., Harbord, R.M., Schmid, C.H., Tetzlaff, J., Deeks, J.J., Peters, J., Macaskill, P., Schwarzer, G., Duval, S., Altman, D.G., Moher, D., Higgins, J.P. (2011)
  Recommendations for examining and interpreting funnel plot asymmetry in meta-analyses of randomised controlled trials.
  \emph{BMJ}, \bold{343}, d4002, doi:10.1136/bmj.d4002.
}
\author{Joaquim Radua}
\seealso{
  \code{\link[=funnel.meta.nsue]{funnel}} for plotting funnel plots.

  \code{\link[=meta.nsue]{meta}} for conducting a meta-analysis.

  \code{\link[=metalm.meta.nsue]{metalm}} for fitting meta-analytic models.
}
\examples{
t <- c(3.4, NA, NA, NA, NA, 2.8, 2.1, 3.1, 2.0, 3.4)
n <- c(40, 20, 22, 24, 18, 30, 25, 30, 16, 22)
m <- meta(smc_from_t(t, n))
metabias(m)
}
\keyword{bias}
\keyword{meta-analysis}
\keyword{regression}
\keyword{test}
