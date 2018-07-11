\name{nsue}
\alias{nsue}
\alias{smc_from_t}
\alias{smd_from_t}
\alias{zcor_from_r}
\alias{print.nsue}
\alias{subset.nsue}
\title{Calculate Effect Sizes for \dQuote{meta.nsue} Objects}
\description{These functions create objects of class \code{"nsue"}, ready to be used by \code{\link{meta}} and \code{\link{leave1out}}.}
\usage{
nsue(y, y_lo = -y_up, y_up, aux, y2var, mi, backtransf = .backtransf_identity, 
     measure = "effect size", labels = "study")
smc_from_t(t, n, alpha = 0.05, labels = "study")
smd_from_t(t, n1, n2, alpha = 0.05, labels = "study")
zcor_from_r(r, n, alpha = 0.05, labels = "study")
}
\arguments{
  \item{y}{a vector to specify the effect-sizes. Use NA in studies with Non Statistically-Significant Unreported Effects (NSUEs).}
  \item{t}{a vector to specify the t-values of the studies. Use NA in studies with Non Statistically-Significant Unreported Effects (NSUEs).}
  \item{r}{a vector to specify the correlation coefficients of the studies. Use NA in studies with Non Statistically-Significant Unreported Effects (NSUEs).}
  \item{y_lo}{a vector to specify the effect-sizes corresponding to the lower statistical threshold.}
  \item{y_up}{a vector to specify the effect-sizes corresponding to the upper statistical threshold.}
  \item{aux}{a data.frame to specify information required for \code{y2var}, \code{mi} and / or \code{backtransf}.}
  \item{n}{a vector to specify the sample sizes of the studies.}
  \item{n1}{a vector to specify the sample sizes of the first group (e.g. patients) of  studies.}
  \item{n2}{a vector to specify the sample sizes of the second group (e.g. controls) of the studies.}
  \item{y2var}{a function to derive the variances of the effect sizes.}
  \item{mi}{a function to multiply impute effect sizes.}
  \item{backtransf}{a function to back-transform the effect sizes.}
  \item{measure}{a description of the effect-size measure used.}
  \item{labels}{a vector to specify the labels of the studies.}
  \item{alpha}{a vector to specify the p-value thresholds used in the studies (e.g. 0.05).}
}
\details{
  Use \code{nsue} for creating an object of class \code{"nsue"}.

  Use \code{smc_from_t} for creating an object of class \code{"nsue"} for standardized mean changes from the t-values of the paired Student t-tests, e.g. in repeated-measures studies analyzing the amount of change in within a group.

  Use \code{smd_from_t} for creating an object of class \code{"nsue"} for standardized mean differences from t-values of the two-sample Student t-tests, e.g. in studies comparing a quantitative (normally-distributed) variable between two groups.

  Use \code{zcor_from_r} for creating an object of class \code{"nsue"} for Pearson correlation coefficients (using the Fisher's transform), e.g. in studies examining the association between two quantitative (normally-distributed) variables.
}
\value{
  \code{nsue}, \code{smc_from_t}, \code{smd_from_t}, and \code{zcor_from_r} return objects of class \code{"nsue"}.

  The function \code{print} may be used to print a summary of the results. The function \code{subset} returns the subset of studies that meets a condition.

  An object of class \code{"nsue"} is a list containing the following components:

  \item{y}{the effect-sizes.}
  \item{y_lo}{the effect-sizes corresponding to the lower statistical threshold.}
  \item{y_up}{the effect-sizes corresponding to the upper statistical threshold.}
  \item{aux}{information required for \code{y2var}, \code{mi} and / or \code{backtransf}.}
  \item{y2var}{a function to derive the variances of the effect sizes.}
  \item{mi}{a function to multiply impute effect sizes.}
  \item{backtransf}{a function to back-transform the effect sizes.}
  \item{measure}{a description of the effect-size measure used.}
  \item{labels}{the labels of the studies.}

  Users can create their objects of class \code{"nsue"} for effect sizes not included in the package.
}
\references{
  Radua, J., Schmidt, A., Borgwardt, S., Heinz, A., Schlagenhauf, F., McGuire, P., Fusar-Poli, P. (2015)
  Ventral striatal activation during reward processing in psychosis. A neurofunctional meta-analysis.
  \emph{JAMA Psychiatry}, \bold{72}, 1243--51.
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
zr <- zcor_from_r(r, n)
m2 <- meta(zr)
zr
m2
}
\keyword{effect-size}
