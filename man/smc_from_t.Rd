\name{smc_from_t}
\alias{smc_from_t}
\alias{smd_from_t}
\alias{z_from_r}
\alias{r_in_smd_from_t_means_and_sds1}
\alias{r_in_smd_from_t_means_and_sds2}
\alias{print.nsue}
\title{Calculate Effect Sizes for \dQuote{meta.nsue} Objects}
\description{These function calculate various effect sizes and output objects of class \code{"nsue"}, ready to be used by \code{\link[=meta.nsue]{meta}} and \code{\link[=leave1out.nsue]{leave1out}}.}
\usage{
smc_from_t(t, n, alpha = 0.05, labels = "study", rm.r = 0.3)
smd_from_t(t, n1, n2, alpha = 0.05, labels = "study", rm.r = 0.3)
z_from_r(r, n, alpha = 0.05, labels = "study", rm.r = 0.3)
r_in_smd_from_t_means_and_sds1(t,
    n1, mean1.pre, sd1.pre, mean1.post, sd1.post,
    n2, mean2.pre, sd2.pre, mean2.post, sd2.post,
    alpha = 0.05, labels = "study", r.range = c(0, 0.99), rm.r = 0.3)
r_in_smd_from_t_means_and_sds2(x, maxiter = 200, tol = 1e-6)
}
\arguments{
  \item{t}{a vector to specify the t-values of the studies. Use NA in studies with Non Statistically-Significant Unreported Effects (NSUEs).}
  \item{r}{a vector to specify the correlation coefficients of the studies. Use NA in studies with Non Statistically-Significant Unreported Effects (NSUEs).}
  \item{n}{a vector to specify the sample sizes of the studies.}
  \item{n1}{a vector to specify the sample sizes of the first group (e.g. patients) of  studies.}
  \item{n2}{a vector to specify the sample sizes of the second group (e.g. controls) of the studies.}
  \item{mean1.pre}{a vector to specify the means of the first group (e.g. patients) of the studies before the treatment.}
  \item{sd1.pre}{a vector to specify the standard deviations of the first group (e.g. patients) of the studies before the treatment.}
  \item{mean1.post}{a vector to specify the means of the first group (e.g. patients) of the studies after the treatment.}
  \item{sd1.post}{a vector to specify the standard deviations of the first group (e.g. patients) of the studies after the treatment.}
  \item{mean2.pre}{a vector to specify the means of the second group (e.g. patients) of the studies before the treatment.}
  \item{sd2.pre}{a vector to specify the standard deviations of the second group (e.g. patients) of the studies before the treatment.}
  \item{mean2.post}{a vector to specify the means of the second group (e.g. patients) of the studies after the treatment.}
  \item{sd2.post}{a vector to specify the standard deviations of the second group (e.g. patients) of the studies after the treatment.}
  \item{alpha}{a vector to specify the p-value thresholds used in the studies (e.g. 0.05).}
  \item{labels}{a vector to specify the names of the studies.}
  \item{r.range}{range of pre-post correlations.}
  \item{rm.r}{the expected correlation coefficient between repeated-measures.}
  \item{x}{an object of class \code{"nsue"}.}
  \item{maxiter}{maximum number of iterations in the REML estimation of \eqn{\tau^2}.}
  \item{tol}{tolerance in the REML estimation of \eqn{\tau^2}.}
}
\details{
  Use \code{smc_from_t} for calculating the standardized mean changes from the t-values of the paired Student t-tests, e.g. in repeated-measures studies analyzing the amount of change in within a group.

  Use \code{smd_from_t} for calculating the standardized mean differences from t-values of the two-sample Student t-tests, e.g. in studies comparing a quantitative (normally-distributed) variable between two groups.

  Use \code{z_from_r} for calculating the Fisher's r-to-z transformed correlations coefficients from the Pearson correlation coefficients (r), e.g. in studies examining the association between two quantitative (normally-distributed) variables.

  Use \code{r_in_smd_from_t_means_and_sds1} and \code{meta} for estimating the missing pre-post correlations in a meta-analysis of the pre-post differences, e.g. when you only have the means and standard deviations before and after a treatment. Afterwards, use \code{r_in_smd_from_t_means_and_sds2} to conduct the meta-analysis of the pre-post differences. Please see Harrison et al for details.
}
\note{
  Studies with the same label will be considered to be repeated-measures, and will thus have lower weight during the MLE step, be imputed accounting for the correlation, and combined before meta-analysis.
}
\value{
  \code{smc_from_t}, \code{smd_from_t}, \code{z_from_r} and \code{r_in_smd_from_t_means_and_sds1} return objects of class \code{"nsue"}.

  The function \code{print} may be used to print a summary of the results.

  An object of class \code{"nsue"} is a list containing several of the following components:

  \item{measure}{the effect-size measure used.}
  \item{y}{the effect-sizes.}
  \item{y_lo}{the effect-sizes corresponding to the lower statistical threshold.}
  \item{y_up}{the effect-sizes corresponding to the upper statistical threshold.}
  \item{y.var}{the variances of the effect sizes.}
  \item{y2var_k1}{a constant needed to derive the variances.}
  \item{y2var_k2}{a constant needed to derive the variances.}
  \item{labels}{the labels of the studies.}
  \item{rm}{a list with the expected correlation between repeated-measures studies.}
}
\references{
  Radua, J., Schmidt, A., Borgwardt, S., Heinz, A., Schlagenhauf, F., McGuire, P., Fusar-Poli, P. (2015)
  Ventral striatal activation during reward processing in psychosis. A neurofunctional meta-analysis.
  \emph{JAMA Psychiatry}, \bold{72}, 1243--51.
  
  Harrison, A., Fernandez de la Cruz, L., Enander, J., Radua, J., Mataix-Cols, D. (2016)
  Cognitive-behavioral therapy for body dysmorphic disorder: A systematic review and meta-analysis of randomized controlled trials.
  \emph{Clinical Psychology Review}, in Press.
}
\author{Joaquim Radua}
\seealso{
  \code{\link[=meta.nsue]{meta}} for conducting a meta-analysis.

  \code{\link[=leave1out.nsue]{leave1out}} for computing leave-one-out diagnostics.
}
\examples{
t <- c(3.4, NA, NA, NA, NA, 2.8, 2.1, 3.1, 2.0, 3.4)
n <- c(40, 20, 22, 24, 18, 30, 25, 30, 16, 22)
x <- smc_from_t(t, n)
x
meta(x)
}
\keyword{effect-size}
