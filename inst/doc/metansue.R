### R code from vignette source 'metansue.Rnw'

###################################################
### code chunk number 1: metansue.Rnw:152-153
###################################################
library(metansue)


###################################################
### code chunk number 2: metansue.Rnw:156-161
###################################################
t <- c(3.4, NA, NA, NA, NA, 2.8, 2.1, 3.1, 2.0, 3.4)
n <- c(40, 20, 22, 24, 18, 30, 25, 30, 16, 22)
x <- smc_from_t(t, n)
x
str(x)


###################################################
### code chunk number 3: metansue.Rnw:183-201
###################################################
t <- c(NA, -2.2, -3.0, -2.1, NA, -2.7, -2.6, -2.4, -3.2)
n1 <- c(18, 15, 22, 17, 14, 23, 29, 24, 45)
mean1.pre <- c(49.9, 49.9, 50.2, 49.8, 50.1, 50.0, 49.9, 51.0, 49.8)
sd1.pre <- c(5.1, 5.2, 6.6, 5.3, 5.5, 3.6, 3.2, 5.9, 5.9)
mean1.post <- c(65.6, 65.7, 66.0, 65.8, 66.6, 65.9, 66.0, 67.0, 65.9)
sd1.post <- c(6.1, 7.3, 5.4, 6.1, 4.6, 7.1, 4.6, 2.8, 7.1)
n2 <- c(18, 18, 23, 18, 16, 23, 27, 23, 45)
mean2.pre <- c(50.8, 50.2, 50.0, 49.7, 49.8, 50.6, 49.8, 50.9, 49.9)
sd2.pre <- c(5.2, 5.9, 4.0, 7.5, 5.7, 4.0, 5.8, 5.6, 5.4)
mean2.post <- c(69.9, 70.3, 70.4, 69.8, 70.2, 70.0, 69.7, 70.8, 69.7)
sd2.post <- c(5.2, 4.7, 4.5, 3.6, 3.7, 6.1, 4.6, 5.9, 6.2)
r <- r_in_smd_from_t_means_and_sds1(t, n1, mean1.pre, sd1.pre,
     mean1.post, sd1.post, n2, mean2.pre, sd2.pre, mean2.post, sd2.post)
r
mr <- meta(r)
mr
m <- r_in_smd_from_t_means_and_sds2(mr)
m


###################################################
### code chunk number 4: metansue.Rnw:398-404
###################################################
t <- c(3.4, NA, NA, NA, NA, 2.8, 2.1, 3.1, 2.0, 3.4)
n <- c(40, 20, 22, 24, 18, 30, 25, 30, 16, 22)
x <- smc_from_t(t, n)
m <- meta(x)
m
str(m)


###################################################
### code chunk number 5: metansue.Rnw:411-416
###################################################
t <- c(3.4, NA, NA, NA, NA, 2.8, 2.1, 3.1, 2.0, 3.4)
n <- c(40, 20, 22, 24, 18, 30, 25, 30, 16, 22)
x <- smc_from_t(t, n)
m <- meta(x)
forest(m)


###################################################
### code chunk number 6: metansue.Rnw:425-430
###################################################
t <- c(3.4, NA, NA, NA, NA, 2.8, 2.1, 3.1, 2.0, 3.4)
n <- c(40, 20, 22, 24, 18, 30, 25, 30, 16, 22)
x <- smc_from_t(t, n)
m <- meta(x)
funnel(m)


###################################################
### code chunk number 7: metansue.Rnw:437-441
###################################################
t <- c(3.4, NA, NA, NA, NA, 2.8, 2.1, 3.1, 2.0, 3.4)
n <- c(40, 20, 22, 24, 18, 30, 25, 30, 16, 22)
x <- smc_from_t(t, n)
summary(leave1out(x))


###################################################
### code chunk number 8: metansue.Rnw:448-453
###################################################
t <- c(3.4, NA, NA, NA, NA, 2.8, 2.1, 3.1, 2.0, 3.4)
n <- c(40, 20, 22, 24, 18, 30, 25, 30, 16, 22)
x <- smc_from_t(t, n)
m <- meta(x)
metabias(m)


###################################################
### code chunk number 9: metansue.Rnw:459-465
###################################################
t <- c(3.4, NA, NA, NA, NA, 2.8, 2.1, 3.1, 2.0, 3.4)
n <- c(40, 20, 22, 24, 18, 30, 25, 30, 16, 22)
x <- smc_from_t(t, n)
x
sx <- subset(x, n > 20)
sx


###################################################
### code chunk number 10: metansue.Rnw:516-524
###################################################
data <- read.table("example1.txt", header = TRUE)
y <- smd_from_t(data$t, data$n1, data$n2, labels = data$study)
m <- meta(y)
m
summary(leave1out(y))
metabias(m)
forest(m)
funnel(m)


###################################################
### code chunk number 11: metansue.Rnw:537-540
###################################################
data <- read.table("example2.txt", header = TRUE)
y <- smd_from_t(data$t, data$n1, data$n2, labels = data$study)
meta(y, ~ data$age)


###################################################
### code chunk number 12: metansue.Rnw:545-565
###################################################
data <- read.table("example1_truth.txt", header = TRUE)
y <- smd_from_t(data$t, data$n1, data$n2, labels = data$study)
meta(y)

library(metafor)
n <- nrow(data)
y <- escalc(measure = "SMD", vtype = "UB",
           m1i = data$t * sqrt(1 / data$n1 + 1 / data$n2),
           n1i = data$n1, n2i = data$n2,
           m2i = rep(0, n), sd1i = rep(1, n), sd2i = rep(1, n))
rma(y)
data <- read.table("example2_truth.txt", header = TRUE)
y <- smd_from_t(data$t, data$n1, data$n2, labels = data$study)
meta(y, ~ data$age)
n <- nrow(data)
y <- escalc(measure = "SMD", vtype = "UB",
           m1i = data$t * sqrt(1 / data$n1 + 1 / data$n2),
           n1i = data$n1, n2i = data$n2,
           m2i = rep(0, n), sd1i = rep(1, n), sd2i = rep(1, n))
rma(y, mods = data$age)


