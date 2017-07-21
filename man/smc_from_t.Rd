\name{smc_from_t}
\alias{smc_from_t}
\alias{smd_from_t}
\alias{z_from_r}
\alias{r_in_smd_from_t_means_and_sds1}
\alias{r_in_smd_from_t_means_and_sds2}
\alias{print.nsue}
\alias{subset.nsue}
\title{Calculate Effect Sizes for \dQuote{meta.nsue} Objects}
\description{These function calculate various effect sizes and output objects of class \code{"nsue"}, ready to be used by \code{\link{meta}} and \code{\link{leave1out}}.}
\usage{
smc_from_t(t, n, alpha = 0.05, labels = "study", rm.r = 0.3)
smd_from_t(t, n1, n2, alpha = 0.05, labels = "study", rm.r = 0.3)
z_from_r(r, n, alpha = 0.05, labels = "study", rm.r = 0.3)
r_in_smd_from_t_means_and_sds1(t,
    n1, mean1.pre, sd1.pre, mean1.post, sd1.post,
    n2, mean2.pre, sd2.pre, mean2.post, sd2.post,
    alpha = 0.05, labels = "study", r.range = c(0, 0.99), rm.r = 0.3)
r_in_smd_from_t_means_and_sds2(x, formula = ~1, hypothesis = NULL,
    maxiter = 200, tol = 1e-6)
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
  \item{formula}{an object of class \code{"formula"}: a symbolic description of the model to be fitted.}
  \item{hypothesis}{a hypothesis, or NULL to test the main coefficient of the model.}
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

  The function \code{print} may be used to print a summary of the results. The function \code{subset} returns the subset of studies that meets a condition.

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
  \code{\link{meta}} for conducting a meta-analysis.

  \code{\link{leave1out}} for computing leave-one-out diagnostics.
}
\examples{
# Standardized mean change in one sample:
t <- c(3.4, NA, NA, NA, 3.2, 2.8, 2.1, 3.1, 2.0, 3.4)
n <- c(40, 20, 22, 24, 18, 30, 25, 30, 16, 22)
smc <- smc_from_t(t, n)
m0 <- meta(smc)
smc
m0

# Standardized mean difference between two samples:
t <- c(4.8, 3.2, NA, NA, NA, 3.2, 2.0, 2.3, 2.7, 3.1)
n1 <- c(40, 20, 22, 24, 18, 30, 25, 30, 16, 22)
n2 <- c(38, 20, 22, 25, 20, 28, 25, 30, 18, 23)
smd <- smd_from_t(t, n1, n2)
m1 <- meta(smd)
smd
m1

# Pearson correlation:
r <- c(0.80, NA, NA, NA, 0.32, 0.45, 0.53, 0.67, 0.74, 0.56)
n <- c(40, 22, 13, 12, 28, 22, 27, 28, 15, 23)
r <- z_from_r(r, n)
m2 <- meta(r)
r
m2

# Pre-post standardized mean difference with missing correlations:
t <- c(3.8, 3.2, NA, NA, NA, 3.1, 2.1, 2.3, 2.7, 3.1)
n1 <- c(40, 20, 22, 24, 18, 30, 25, 30, 16, 22)
mean1.pre <- c(8.7, 19.3, 16.1, 14.2, 21.2, 10.7, 12.1, 5.9, 15.8, 13.6)
sd1.pre <- c(1.7, 3.3, 2.0, 2.0, 1.7, 2.5, 2.1, 2.1, 1.7, 2.1)
mean1.post <- c(9.7, 21.9, 17.6, 15.6, 22.9, 12.4, 13.3, 7.3, 18.0, 15.4)
sd1.post <- c(2.8, 4.2, 3.3, 3.2, 2.9, 2.6, 3.2, 3.0, 3.2, 3.4)
n2 <- c(38, 20, 22, 25, 20, 28, 25, 30, 18, 23)
mean2.pre <- c(17.3, 9.4, 12.0, 13.5, 13.9, 19.0, 9.4, 14.8, 22.8, 26.2)
sd2.pre <- c(1.8, 2.1, 1.3, 1.8, 1.7, 0.3, 2.3, 2.1, 2.5, 2.5)
mean2.post <- c(16.6, 7.4, 11.0, 12.5, 13.4, 18.8, 8.9, 14.6, 22.2, 25.2)
sd2.post <- c(2.5, 4.2, 3.6, 4.5, 3.4, 2.2, 3.7, 3.5, 3.0, 3.2)
r <- r_in_smd_from_t_means_and_sds1(t,
    n1, mean1.pre, sd1.pre, mean1.post, sd1.post,
    n2, mean2.pre, sd2.pre, mean2.post, sd2.post)
m3_cor <- meta(r)
m3 <- r_in_smd_from_t_means_and_sds2(m3_cor)
r
m3
}
\keyword{effect-size}
