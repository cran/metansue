.check.alpha <-
function (call, alpha, n.stud) 
{
    if (!is.numeric(alpha)) {
        .stop(call, "alpha must be a numeric vector")
    }
    if (any(is.na(alpha))) {
        .stop(call, "alpha cannot have missing values")
    }
    if (any(alpha <= 0) || any(alpha >= 1)) {
        .stop(call, "alpha cannot be <= 0 or >= 1")
    }
    if (length(alpha) == 1) {
        return(rep(alpha, n.stud))
    }
    if (length(alpha) != n.stud) {
        .stop(call, "alpha has an incorrect length")
    }
    alpha
}
.check.formula <-
function (call, formula, n.stud) 
{
    if (!inherits(formula, "formula")) {
        .stop(call, "formula must be a formula")
    }
    terms <- terms(formula)
    xnames <- attr(terms, "term.labels")
    if (length(xnames) == 0) {
        return(list(formula = "~ 1", matrix = matrix(1, n.stud), 
            labels = "(Mean)"))
    }
    formula <- paste("~", paste(xnames, collapse = " + "))
    if (!attr(terms, "intercept")) {
        warning("You have specified a regression though the origin")
        formula <- paste(formula, "- 1")
    }
    current.na.action <- options("na.action")$na.action
    options(na.action = "na.pass")
    X <- model.matrix(as.formula(formula), parent.frame(2))
    options(na.action = current.na.action)
    if (nrow(X) != n.stud) {
        .stop(call, "Independent variables of the formula have an incorrect length")
    }
    if (any(is.na(X))) {
        .stop(call, "Independent variables of the formula cannot have missing values. Impute missing values (using for example the R package 'mi'), call 'meta' for each imputation, and combine all imputations")
    }
    list(formula = formula, matrix = matrix(c(X), n.stud), labels = colnames(X))
}
.check.hypothesis <-
function (call, hypothesis, model) 
{
    n.coef <- ncol(model$matrix)
    labels <- model$labels
    if (is.null(hypothesis)) {
        if (n.coef == 1) {
            hypothesis = list(text = paste(labels[1], "=0", sep = ""), 
                matrix = matrix(1))
        }
        else {
            hypothesis = list(text = paste(labels[2], "=0", sep = ""), 
                matrix = matrix(c(0, 1, rep(0, n.coef - 2)), 
                  1))
        }
    }
    else if (is.matrix(hypothesis)) {
        if (ncol(hypothesis) != n.coef) {
            .stop(call, "Wrong number of columns in the hypothesis")
        }
        text = c()
        for (i in 1:nrow(hypothesis)) {
            text_i = ""
            for (j in 1:ncol(hypothesis)) {
                hypothesis_j = hypothesis[i, j]
                if (hypothesis_j != 0) {
                  if (hypothesis_j > 0 && nchar(text_i) > 0) {
                    text_i = paste(text_i, "+", sep = "")
                  }
                  else if (hypothesis_j < 0) {
                    text_i = paste(text_i, "-", sep = "")
                  }
                  if (abs(hypothesis_j != 1)) {
                    text_i = paste(text_i, hypothesis_j, "*", 
                      sep = "")
                  }
                  text_i = paste(text_i, labels[j], sep = "")
                }
            }
            text = c(text, paste(text_i, "=0", sep = ""))
        }
        if (nrow(hypothesis) > 1) {
            warning("All rows of the hypothesis are given the same weight in the MLE step; please adjust if required.")
        }
        hypothesis = list(text = text, matrix = hypothesis)
    }
    else if (is.numeric(hypothesis)) {
        if (length(hypothesis) != n.coef) {
            .stop(call, "Wrong vector length in the hypothesis")
        }
        text = ""
        for (j in 1:length(hypothesis)) {
            hypothesis_j = hypothesis[j]
            if (hypothesis_j != 0) {
                if (hypothesis_j > 0 && nchar(text) > 0) {
                  text = paste(text, "+", sep = "")
                }
                else if (hypothesis_j < 0) {
                  text = paste(text, "-", sep = "")
                }
                if (abs(hypothesis_j != 1)) {
                  text = paste(text, hypothesis_j, "*", sep = "")
                }
                text = paste(text, labels[j], sep = "")
            }
        }
        hypothesis = list(text = paste(text, "=0", sep = ""), 
            matrix = matrix(hypothesis, 1))
    }
    else {
        .stop(call, "Numeric vector or matrix expected in the hypothesis")
    }
    hypothesis
}
.check.labels <-
function (call, labels, n.stud) 
{
    if (!is.vector(labels) && !is.factor(labels)) {
        .stop(call, "labels must be a vector")
    }
    if (length(labels) == 1) {
        return(paste0(labels, 1:n.stud))
    }
    if (length(labels) != n.stud) {
        .stop(call, "labels has an incorrect length")
    }
    as.character(labels)
}
.check.n <-
function (call, n, min.n, n.stud) 
{
    if (!is.numeric(n)) {
        .stop(call, "n must be a numeric vector")
    }
    if (any(is.na(n))) {
        .stop(call, "n cannot have missing values")
    }
    if (any(n < min.n)) {
        .stop(call, paste("n cannot be <", min.n))
    }
    if (length(n) == 1) {
        return(rep(n, n.stud))
    }
    if (length(n) != n.stud) {
        .stop(call, "n has an incorrect length")
    }
    n
}
.d_j <-
function (x) 
{
    j <- ifelse(x < 1, -1, 1) * exp(lgamma(x/2) - 0.5 * log(x/2) - 
        lgamma((x - 1)/2))
    na.j <- which(is.na(j))
    j[na.j] <- 1 - 3/(4 * x[na.j] - 1)
    j
}
.elliptic.q <-
function (x, y, p = 0.95, col = "#cccccc", segments = 51) 
{
    center <- c(mean(x), mean(y))
    shape <- cov(cbind(x, y))
    radius <- sqrt(2 * qf(p, 2, length(x) - 1))
    angles <- (0:segments) * 2 * pi/segments
    circle <- cbind(cos(angles), sin(angles))
    choleski <- chol(shape, pivot = TRUE)
    polygon(t(center + radius * t(circle %*% choleski[, order(attr(choleski, 
        "pivot"))])), col = col, border = NA)
}
.estimate_n.mle_discard <-
function (N) 
{
    lof = 1
    while (N > switch(as.character(lof), `1` = 30, `2` = 59, 
        `3` = 77, 34 + lof * 15)) {
        lof = lof + 1
    }
    lof
}
.format.0pos <-
function (x) 
{
    formatC(x, 0, width = 2, format = "f")
}
.format.1 <-
function (x) 
{
    formatC(x, 1, width = 4, format = "f")
}
.format.1pos <-
function (x) 
{
    formatC(x, 1, width = 3, format = "f")
}
.format.2 <-
function (x) 
{
    formatC(x, 2, width = 5, format = "f")
}
.format.2pos <-
function (x) 
{
    formatC(x, 2, width = 4, format = "f")
}
.format.3 <-
function (x) 
{
    formatC(x, 3, width = 6, format = "f")
}
.format.4 <-
function (x) 
{
    formatC(x, 4, width = 7, format = "f")
}
.format.4pos <-
function (x) 
{
    formatC(x, 4, width = 6, format = "f")
}
.format.measure <-
function (measure) 
{
    switch(measure, cor = "correlation", `cor in smd` = "correlation in standardized mean difference", 
        smc = "standardized mean change", smd = "standardized mean difference")
}
.format.perc2 <-
function (x) 
{
    formatC(100 * x, 2, width = 5, format = "f")
}
.format.prob <-
function (p) 
{
    p <- formatC(p, digits = 4, format = "f")
    p[p == "0.0000"] <- "<.0001"
    p
}
.format.sign <-
function (p) 
{
    symnum(p, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", 
        "**", "*", ".", " "), na = FALSE, corr = FALSE)
}
.meta.nsue <-
function (x, model, hypothesis, n.imp, n.bins, maxiter, tol) 
{
    measure <- x$measure
    y <- x$y
    x$y = NULL
    y_lo = x$y_lo
    x$y_lo = NULL
    y_up = x$y_up
    x$y_up = NULL
    n.stud <- length(y)
    known <- which(!is.na(y) | y_up - y_lo < 1e-06)
    unknown <- setdiff(1:n.stud, known)
    y[is.na(y)] = (y_up[is.na(y)] + y_lo[is.na(y)])/2
    if (measure == "cor" || measure == "cor in smd") {
        y.var = x$y.var
        y_lo.var = y.var
        y_up.var = y.var
    }
    if (measure == "smc" || measure == "smd") {
        y2var_k1 <- x$y2var_k1
        y2var_k2 <- x$y2var_k2
        y.var = y2var_k1 + y2var_k2 * y^2
        y_lo.var = y2var_k1 + y2var_k2 * y_lo^2
        y_up.var = y2var_k1 + y2var_k2 * y_up^2
    }
    X = model$matrix
    n.coef = ncol(X)
    if (length(unknown)) {
        if (measure == "cor" || measure == "cor in smd") {
            mll_coef <- function(coef, known, known.y, known.y.var, 
                known.weights, unknown, unknown.y_lo, unknown.y_lo.se, 
                unknown.y_up, unknown.y_up.se, unknown.weights, 
                X) {
                mu <- X %*% coef
                unknown.mu <- mu[unknown]
                a = pnorm((unknown.y_up - unknown.mu)/unknown.y_up.se, 
                  log.p = TRUE)
                c = pnorm((unknown.y_lo - unknown.mu)/unknown.y_up.se, 
                  log.p = TRUE)
                sum(known.weights * (log(known.y.var) + (known.y - 
                  mu[known])^2/known.y.var))/2 - sum(unknown.weights * 
                  (a + log(-expm1(c - a))))
            }
            mll_tau2 <- function(tau2, known.err2, known.y.var, 
                known.weights, unknown.err_lo, unknown.y_lo.var, 
                unknown.err_up, unknown.y_up.var, unknown.weights) {
                if (tau2 < 0) {
                  return(Inf)
                }
                known.sigma2 <- known.y.var + tau2
                a = pnorm(unknown.err_up/sqrt(unknown.y_up.var + 
                  tau2), log.p = TRUE)
                c = pnorm(unknown.err_lo/sqrt(unknown.y_lo.var + 
                  tau2), log.p = TRUE)
                sum(known.weights * (log(known.sigma2) + known.err2/known.sigma2))/2 - 
                  sum(unknown.weights * (a + log(-expm1(c - a))))
            }
            mi1 <- function(tau2, mu, y_lo, y_up, y.var, rm.var, 
                rm.y, n.imp) {
                sigma2 <- y.var + tau2
                sigma <- sqrt(sigma2)
                if (length(rm.var)) {
                  rm.sigma <- sqrt(rm.var + tau2)
                  Sxx <- rm.sigma %*% t(rm.sigma) * (diag(1 - 
                    rm$r, length(rm.var)) + rm$r)
                  Sxy <- rm.sigma * sigma * rm$r
                  beta <- solve(Sxx) %*% Sxy
                  mus <- rm.y %*% beta
                  sigma <- sqrt(sigma2 - t(beta) %*% Sxx %*% 
                    beta)
                }
                else {
                  mus <- rep(mu, n.imp)
                }
                q <- rep(NA, n.imp)
                to_imp <- 1:n.imp
                while (length(to_imp)) {
                  q[to_imp] <- rnorm(length(to_imp), mus[to_imp], 
                    sigma)
                  to_imp <- which(q <= y_lo & q >= y_up)
                }
                q
            }
        }
        if (measure == "smc" || measure == "smd") {
            mll_coef <- function(coef, known, known.y, known.y.var, 
                known.weights, unknown, unknown.y_lo, unknown.y_lo.se, 
                unknown.y_up, unknown.y_up.se, unknown.weights, 
                X) {
                mu <- X %*% coef
                unknown.mu <- mu[unknown]
                a = pnorm((unknown.y_up - unknown.mu)/unknown.y_up.se, 
                  log.p = TRUE)
                c = pnorm((unknown.y_lo - unknown.mu)/unknown.y_lo.se, 
                  log.p = TRUE)
                sum(known.weights * (log(known.y.var) + (known.y - 
                  mu[known])^2/known.y.var))/2 - sum(unknown.weights * 
                  (a + log(-expm1(c - a))))
            }
            mll_tau2 <- function(tau2, known.err2, known.y.var, 
                known.weights, unknown.err_lo, unknown.y_lo.var, 
                unknown.err_up, unknown.y_up.var, unknown.weights) {
                if (tau2 < 0) {
                  return(Inf)
                }
                known.sigma2 <- known.y.var + tau2
                a = pnorm(unknown.err_up/sqrt(unknown.y_up.var + 
                  tau2), log.p = TRUE)
                c = pnorm(unknown.err_lo/sqrt(unknown.y_lo.var + 
                  tau2), log.p = TRUE)
                sum(known.weights * (log(known.sigma2) + known.err2/known.sigma2))/2 - 
                  sum(unknown.weights * (a + log(-expm1(c - a))))
            }
            mi2 <- function(tau2, mu, y_lo, y_up, y2var_k1, y2var_k2, 
                rm.var_k1, rm.var_k2, rm.y, n.imp, n.bins) {
                sigma2 <- y2var_k1 + y2var_k2 * mu^2 + tau2
                sigma <- sqrt(sigma2)
                if (length(rm.var_k1)) {
                  rm.sigma <- sqrt(rm.var_k1 + rm.var_k2 * mu^2 + 
                    tau2)
                  Sxx <- rm.sigma %*% t(rm.sigma) * (diag(1 - 
                    rm$r, length(rm.var_k1)) + rm$r)
                  Sxy <- rm.sigma * sigma * rm$r
                  beta <- solve(Sxx) %*% Sxy
                  mus <- rm.y %*% beta
                  sigma <- sqrt(sigma2 - t(beta) %*% Sxx %*% 
                    beta)
                }
                else {
                  mus <- rep(mu, n.imp)
                }
                width <- (y_up - y_lo)/n.bins
                y <- sample(seq(y_lo + width/2, y_up - width/2, 
                  width))
                q <- c()
                for (imp in 1:n.imp) {
                  raw_dens <- dnorm(y, mus[imp], sigma) * (y2var_k1 + 
                    y2var_k2 * y^2 + tau2)
                  pfun <- cumsum(raw_dens/sum(raw_dens))
                  p <- runif(1)
                  j <- 1
                  while (p > pfun[j]) {
                    j <- j + 1
                  }
                  q <- c(q, y[j])
                }
                q
            }
        }
    }
    rm <- x$rm
    rm.M <- t(unname(model.matrix(~0 + x$labels)))
    rm.M <- rm.M[unique(c(1:nrow(rm.M) %*% rm.M)), ]
    rm_weights = apply(apply(rm.M, 1, function(x) {
        x/(1 + (sum(x) - 1) * rm$r)
    }), 1, sum)
    rm$M <- rm.M
    rm$weights <- rm_weights
    x$rm = rm
    if (length(unknown)) {
        coef = NULL
        n.mle_discard = .estimate_n.mle_discard(n.stud)
        mle_discard = c()
        while (length(mle_discard) < n.mle_discard) {
            min_hyp = Inf
            for (i in setdiff(1:n.stud, mle_discard)) {
                sample_i = (1:n.stud)[-c(mle_discard, i)]
                known_i = sample_i[sample_i %in% known]
                unknown_i = sample_i[sample_i %in% unknown]
                if (n.coef == 1) {
                  interval = c(min(c(y[known_i], y_lo[unknown_i])), 
                    max(c(y[known_i], y_up[unknown_i])))
                  coef_i <- optimize(mll_coef, interval, known_i, 
                    y[known_i], y.var[known_i], rm_weights[known_i], 
                    unknown_i, y_lo[unknown_i], sqrt(y_lo.var[unknown_i]), 
                    y_up[unknown_i], sqrt(y_up.var[unknown_i]), 
                    rm_weights[unknown_i], X)$minimum
                }
                else {
                  initial_coef <- coef(lm.wfit(X[sample_i, ], 
                    y[sample_i], 1/y.var[sample_i]))
                  if (measure == "cor" || measure == "cor in smd") {
                    coef_i <- optim(initial_coef, mll_coef, gr = NULL, 
                      known_i, y[known_i], y.var[known_i], rm_weights[known_i], 
                      unknown_i, y_lo[unknown_i], sqrt(y_lo.var[unknown_i]), 
                      y_up[unknown_i], sqrt(y_up.var[unknown_i]), 
                      rm_weights[unknown_i], X)$par
                  }
                  if (measure == "smc" || measure == "smd") {
                    coef_i <- optim(initial_coef, mll_coef, gr = NULL, 
                      known_i, y[known_i], y.var[known_i], rm_weights[known_i], 
                      unknown_i, y_lo[unknown_i], sqrt(y_lo.var[unknown_i]), 
                      y_up[unknown_i], sqrt(y_up.var[unknown_i]), 
                      rm_weights[unknown_i], X)$par
                  }
                }
                hyp_i = sum(abs(hypothesis$matrix %*% coef_i))
                if (hyp_i < abs(min_hyp)) {
                  min_hyp = hyp_i
                  min_i = i
                  min_coef = coef_i
                }
            }
            mle_discard = c(mle_discard, min_i)
        }
        mu = X %*% min_coef
        tau2 <- optimize(mll_tau2, c(0, 999), (y[known] - mu[known])^2, 
            y.var[known], rm_weights[known], y_lo[unknown] - 
                mu[unknown], y_lo.var[unknown], y_up[unknown] - 
                mu[unknown], y_up.var[unknown], rm_weights[unknown])$minimum
        mi_y <- NULL
        for (i in unknown) {
            rm.indexs <- which(x$labels == x$labels[i] & 1:n.stud < 
                i)
            if (measure == "cor" || measure == "cor in smd") {
                rm.var <- c()
            }
            if (measure == "smc" || measure == "smd") {
                rm.var_k1 <- c()
                rm.var_k2 <- c()
            }
            rm.y = NULL
            for (j in rm.indexs) {
                is.known = !is.na(y[j])
                if (measure == "cor" || measure == "cor in smd") {
                  rm.var <- c(rm.var, y.var[j])
                }
                if (measure == "smc" || measure == "smd") {
                  rm.var_k1 <- c(rm.var_k1, y2var_k1[j])
                  rm.var_k2 <- c(rm.var_k2, y2var_k2[j])
                }
                if (is.known) {
                  rm.y <- cbind(rm.y, rep(y[j], n.imp))
                }
                else {
                  rm.y <- cbind(rm.y, mi_y[match(j, unknown), 
                    ])
                }
            }
            if (measure == "cor" || measure == "cor in smd") {
                mi_y <- rbind(mi_y, mi1(tau2, mu[i], y_lo[i], 
                  y_up[i], y.var[i], rm.var, rm.y, n.imp))
            }
            if (measure == "smc" || measure == "smd") {
                mi_y <- rbind(mi_y, mi2(tau2, mu[i], y_lo[i], 
                  y_up[i], y2var_k1[i], y2var_k2[i], rm.var_k1, 
                  rm.var_k2, rm.y, n.imp, n.bins))
            }
        }
        colnames(mi_y) <- NULL
    }
    else {
        mi_y = matrix(nrow = 0, ncol = 0)
    }
    x$known = list(i = known, y = y[known])
    x$unknown = list(i = unknown, y = mi_y)
    class(x) <- "meta.nsue"
    .meta.nsue2(x, model, hypothesis, maxiter, tol)
}
.meta.nsue2 <-
function (x, model, hypothesis, maxiter, tol) 
{
    measure = x$measure
    known = x$known$i
    unknown = x$unknown$i
    y = rep(NA, length(known) + length(unknown))
    y[known] = x$known$y
    mi_y = x$unknown$y
    if (measure == "cor" || measure == "cor in smd") {
        y.var = x$y.var
    }
    if (measure == "smc" || measure == "smd") {
        y2var_k1 <- x$y2var_k1
        y2var_k2 <- x$y2var_k2
    }
    rm = x$rm
    rm.M = rm$M
    rm_weights = rm$weights
    X = model$matrix
    n.coef = ncol(X)
    df <- nrow(rm.M) - n.coef
    rm.M2 <- t(apply(rm.M, 1, function(x) {
        x/sum(x)
    }))
    if (measure == "cor" || measure == "cor in smd") {
        y.var <- y.var/rm_weights
        ay.var <- apply(rm.M, 1, function(xx) {
            mean(y.var[which(xx == 1)])/sum(xx)
        })
    }
    if (measure == "smc" || measure == "smd") {
        y2var_k1 <- y2var_k1/rm_weights
        ay2var_k1 <- apply(rm.M, 1, function(xx) {
            mean(y2var_k1[which(xx == 1)])/sum(xx)
        })
        y2var_k2 <- y2var_k2
        ay2var_k2 <- apply(rm.M, 1, function(xx) {
            mean(y2var_k2[which(xx == 1)])
        })
    }
    aX <- rm.M2 %*% X
    mi.coef <- NULL
    mi.cov <- NULL
    mi.tau2 <- c()
    mi.qe <- c()
    mi.i2 <- c()
    mi.h2 <- c()
    for (j in 1:max(c(1, ncol(mi_y)))) {
        if (ncol(mi_y) > 0) {
            y[unknown] <- mi_y[, j]
        }
        ay <- c(rm.M2 %*% y)
        if (measure == "smc" || measure == "smd") {
            ay.var <- ay2var_k1 + ay2var_k2 * ay^2
        }
        W_fe <- diag(1/ay.var)
        P_fe <- W_fe - W_fe %*% aX %*% solve(t(aX) %*% W_fe %*% 
            aX) %*% t(aX) %*% W_fe
        tau2_j <- .tau2.reml(ay, ay.var, aX, maxiter, tol)
        W <- diag(1/(ay.var + tau2_j))
        inv_XtWX <- solve(t(aX) %*% W %*% aX)
        h2_j <- 1 + tau2_j/df * sum(diag(P_fe))
        mi.coef <- cbind(mi.coef, inv_XtWX %*% t(aX) %*% W %*% 
            ay)
        mi.cov <- cbind(mi.cov, c(inv_XtWX))
        mi.tau2 <- c(mi.tau2, tau2_j)
        mi.qe <- c(mi.qe, max(0, t(ay) %*% P_fe %*% ay))
        mi.i2 <- c(mi.i2, max(0, 1 - 1/h2_j))
        mi.h2 <- c(mi.h2, h2_j)
    }
    coef <- apply(mi.coef, 1, mean)
    cov <- .pool.cov(mi.coef, mi.cov)
    tau2 = mean(mi.tau2)
    f_df2 <- .pool.chi2(mi.qe, df)
    f <- f_df2[1]
    df2 <- f_df2[2]
    i2 = mean(mi.i2)
    h2 = mean(mi.h2)
    model$coef <- coef
    model$cov <- cov
    model$se <- sqrt(diag(cov))
    x$model = model
    x$heterogeneity <- list(tau2 = tau2, h2 = h2, i2 = i2, q = data.frame(f, 
        df1 = df, df2, p.value = 1 - pf(f, df, df2)))
    h <- hypothesis$matrix
    hcoef <- c(h %*% coef)
    hcov <- h %*% cov %*% t(h)
    hypothesis$coef <- hcoef
    if (nrow(h) == 1) {
        hse <- sqrt(hcov)
        z <- hcoef/hse
        hypothesis$se <- hse
        hypothesis$z <- z
        hypothesis$p.value <- 2 * pnorm(-abs(z))
    }
    else {
        qr <- qr(hcov)
        pivot <- qr$pivot[1:qr$rank]
        chisq <- c(t(hcoef[pivot]) %*% solve(hcov[pivot, pivot]) %*% 
            hcoef[pivot])
        df <- length(pivot)
        hypothesis$chisq <- chisq
        hypothesis$df <- df
        hypothesis$p.value <- 1 - pchisq(chisq, df)
    }
    x$hypothesis = hypothesis
    x
}
.pool.chi2 <-
function (chi2, df1) 
{
    m <- length(chi2)
    if (m == 1) {
        return(c(chi2/df1, Inf))
    }
    r <- (1 + 1/m) * var(sqrt(chi2))
    c((mean(chi2)/df1 - (m + 1)/(m - 1) * r)/(1 + r), (m - 1)/df1^(3/m) * 
        (1 + 1/r)^2)
}
.pool.cov <-
function (x, x_cov) 
{
    m <- ncol(x)
    if (m == 1) {
        return(matrix(x_cov, nrow(x)))
    }
    cov0 <- matrix(apply(x_cov, 1, mean), nrow(x))
    var0 <- diag(cov0)
    var.increase <- sqrt((var0 + (1 + 1/m) * apply(x, 1, var))/var0)
    var.increase %*% t(var.increase) * cov0
}
.pool.var <-
function (x, x_var) 
{
    if (is.vector(x)) {
        m <- length(x)
        if (m == 1) {
            return(x_var)
        }
        return(mean(x_var) + (1 + 1/length(x)) * var(x))
    }
    m <- ncol(x)
    if (m == 1) {
        return(x_var)
    }
    apply(x_var, 1, mean) + (1 + 1/ncol(x)) * apply(x, 1, var)
}
.print.heterogeneity <-
function (x) 
{
    cat("Residual heterogeneity (tau^2):", .format.4pos(x$heterogeneity$tau2), 
        "  I^2:", paste(.format.perc2(x$heterogeneity$i2), "%", 
            sep = ""), "  H^2:", .format.2pos(x$heterogeneity$h2), 
        "\n")
    cat("F-statistic (heterogeneity):", .format.2pos(x$heterogeneity$q$f), 
        "on", x$heterogeneity$q$df1, "and", .format.1pos(x$heterogeneity$q$df2), 
        "df  Pr(>F):", .format.prob(x$heterogeneity$q$p.value), 
        "\n")
}
.print.hypothesis <-
function (x) 
{
    coef <- x$hypothesis$coef
    nrow = length(coef)
    p.value <- x$hypothesis$p.value
    prob <- .format.prob(p.value)
    sign <- .format.sign(p.value)
    if (nrow == 1) {
        measure <- x$measure
        se <- x$hypothesis$se
        cat("One-row hypothesis:\n")
        ci.low <- coef + qnorm(0.025) * se
        ci.up <- coef + qnorm(0.975) * se
        if (measure == "cor" || measure == "cor in smd") {
            hypothesis <- cbind(cbind(.format.3(tanh(coef)), 
                .format.4(x$hypothesis$z), prob, .format.3(tanh(cbind(ci.low, 
                  ci.up)))), sign)
            colnames(hypothesis) <- c("Corr", "z value", "Pr(>|z|)", 
                "CI(low)", "CI(up)", "")
        }
        if (measure == "smc" || measure == "smd") {
            hypothesis <- cbind(cbind(.format.4(coef), .format.4(x$hypothesis$z), 
                prob, .format.4(cbind(ci.low, ci.up))), sign)
            colnames(hypothesis) <- c("Estimate", "z value", 
                "Pr(>|z|)", "CI(low)", "CI(up)", "")
        }
    }
    else {
        cat("Multi-row hypothesis:\n")
        hypothesis <- cbind(.format.4(coef), c(.format.2pos(x$hypothesis$chisq), 
            rep("", nrow - 1)), c(.format.0pos(x$hypothesis$df), 
            rep("", nrow - 1)), c(prob, rep("", nrow - 1)), c(sign, 
            rep("", nrow - 1)))
        colnames(hypothesis) <- c("Estimate", "chisq", "df", 
            "Pr(>chisq)", "")
    }
    rownames(hypothesis) <- x$hypothesis$text
    print(hypothesis, quote = FALSE, right = TRUE, print.gap = 2)
}
.print.model <-
function (x) 
{
    cat("Model:\n")
    table <- cbind(.format.4(x$model$coef), .format.4pos(x$model$se))
    colnames(table) <- c("Estimate", "Std. Error")
    rownames(table) <- x$model$labels
    print(table, quote = FALSE, right = TRUE, print.gap = 2)
}
.print.sign <-
function () 
{
    cat("---\n")
    cat("Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
}
.r_in_smd_from_sds <-
function (var.diff, df1, var1.sum, sd1.prod, df2, var2.sum, sd2.prod, 
    r.min, r.max) 
{
    r <- (df1 * (var1.sum - var.diff) + df2 * (var2.sum - var.diff))/(2 * 
        (df1 * sd1.prod + df2 * sd2.prod))
    r[r < r.min] = r.min
    r[r > r.max] = r.max
    r
}
.r_in_smd_from_t_means_and_sds2 <-
function (x, model, hypothesis, maxiter, tol) 
{
    x$measure <- "smd"
    smd <- x$smd[x$known$i, ]
    x$known$y <- smd$j * smd$diff/sqrt((smd$df1 * (smd$var1.sum - 
        2 * smd$sd1.prod * tanh(x$known$y)) + smd$df2 * (smd$var2.sum - 
        2 * smd$sd2.prod * tanh(x$known$y)))/(smd$df1 + smd$df2))
    smd <- x$smd[x$unknown$i, ]
    if (length(x$unknown$y)) {
        x$unknown$y <- smd$j * smd$diff/sqrt((smd$df1 * (smd$var1.sum - 
            2 * smd$sd1.prod * tanh(x$unknown$y)) + smd$df2 * 
            (smd$var2.sum - 2 * smd$sd2.prod * tanh(x$unknown$y)))/(smd$df1 + 
            smd$df2))
    }
    x$y2var_k1 = x$smd$y2var_k1
    x$y2var_k2 = x$smd$y2var_k2
    x$y.var <- NULL
    x$smd <- NULL
    .meta.nsue2(x, model, hypothesis, maxiter, tol)
}
.signif.up <-
function (x, digits = 6) 
{
    power <- 10^(round(digits) - 1 - floor(log10(abs(x))))
    ceiling(x * power)/power
}
.stop <-
function (call, message) 
{
    cat("\n")
    print(call)
    stop(paste0(message, "\n "), call. = FALSE)
}
.tau2.reml <-
function (y, y.var, X, maxiter, tol) 
{
    tau2 = var(y) - mean(y.var)
    if (tau2 < 0) {
        tau2 = 0
    }
    for (iter in 1:maxiter) {
        old_tau2 <- tau2
        W <- diag(1/(y.var + tau2))
        P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% 
            W
        tau2 <- max(0, tau2 + solve(sum(diag(P %*% P))) %*% (t(y) %*% 
            P %*% P %*% y - sum(diag((P)))))
        if (abs(tau2 - old_tau2) < tol) {
            break
        }
    }
    tau2
}
.warning <-
function (message) 
{
    warning(message, call. = FALSE)
}
coef.meta.nsue <-
function (object, ...) 
{
    table <- cbind(object$model$coef, object$model$se)
    colnames(table) <- c("Estimate", "Std. Error")
    table
}
fitted.meta.nsue <-
function (object, ...) 
{
    object$model$matrix %*% object$model$coef
}
forest <-
function (x, ...) 
UseMethod("forest")
forest.meta.nsue <-
function (x, ...) 
{
    if (!inherits(x, "meta.nsue")) {
        .stop(match.call(), "The argument must be a 'meta.nsue' object")
    }
    measure = x$measure
    known = x$known$i
    unknown = x$unknown$i
    unknown.n.stud <- length(unknown)
    n.stud <- length(known) + unknown.n.stud
    if (nrow(x$rm$M) < n.stud) {
        .warning("This plot shows repeated measures as separate studies")
    }
    if (length(x$hypothesis$coef) > 1) {
        .warning("This plot only shows the first row of the hypothesis")
    }
    labels <- c(x$labels[unknown], x$labels[known], x$hypothesis$text[1])
    pos.y <- c(n.stud + 2 - c(unknown, known), 0)
    if (unknown.n.stud) {
        y <- apply(x$unknown$y, 1, mean)
        y.low <- apply(x$unknown$y, 1, function(x) {
            quantile(x, 0.025)
        })
        y.upp <- apply(x$unknown$y, 1, function(x) {
            quantile(x, 0.975)
        })
        if (measure == "cor" || measure == "cor in smd") {
            y.se <- sqrt(x$y.var[unknown])
        }
        if (measure == "smc" || measure == "smd") {
            y.se <- sqrt(.pool.var(x$unknown$y, x$y2var_k1[unknown] + 
                x$y2var_k2[unknown] * x$unknown$y^2))
        }
    }
    else {
        y <- y.low <- y.upp <- y.se <- c()
    }
    y <- c(y, x$known$y, x$hypothesis$coef[1])
    if (measure == "cor" || measure == "cor in smd") {
        y.se <- c(y.se, sqrt(x$y.var[known]), x$hypothesis$se[1])
    }
    if (measure == "smc" || measure == "smd") {
        y.se <- c(y.se, sqrt(x$y2var_k1[known] + x$y2var_k2[known] * 
            x$known$y^2), x$hypothesis$se[1])
    }
    ci.low <- y + qnorm(0.025) * y.se
    ci.upp <- y + qnorm(0.975) * y.se
    if (measure == "cor" || measure == "cor in smd") {
        y <- tanh(y)
        if (unknown.n.stud) {
            y.low <- tanh(y.low)
            y.upp <- tanh(y.upp)
        }
        ci.low <- tanh(ci.low)
        ci.upp <- tanh(ci.upp)
    }
    lwd <- 1/y.se
    lwd <- sqrt(9 + 216 * (lwd - min(lwd))/(max(lwd) - min(lwd)))
    ci.text <- paste0(.format.2(y), " [ ", .format.2(ci.low), 
        ", ", .format.2(ci.upp), " ] ", .format.sign(2 * pnorm(-abs(y/y.se))))
    plot.new()
    xlim <- c(-2.5 - max(strwidth(labels, units = "inches")), 
        max(strwidth(ci.text, units = "inches")) + 2.5)
    ylim <- c(-2, n.stud + 1)
    plot.window(xlim = xlim, ylim = ylim)
    xthr <- .signif.up(max(abs(c(quantile(ci.low, 0.1), quantile(ci.upp, 
        0.9)))), 1)
    lines(rep(0, 2), c(n.stud + 1.5, -1), col = "#bbbbbb", lty = 1)
    lines(c(-2, 2), rep(-1, 2), col = "#bbbbbb", lty = 1)
    for (pos.x in -2:2) {
        lines(rep(pos.x, 2), c(-1, -1.3), col = "#bbbbbb", lty = 1)
        text(pos.x, -1.5, .format.2(pos.x/2 * xthr), pos = 1, 
            col = "#bbbbbb")
    }
    for (i in 1:(n.stud + 1)) {
        pos.y_i <- pos.y[i]
        y_i <- y[i]
        ci.low_i <- ci.low[i]
        ci.upp_i <- ci.upp[i]
        if (i > unknown.n.stud) {
            col <- "#000000"
        }
        else {
            y.low_i <- y.low[i]
            y.upp_i <- y.upp[i]
            if (y.upp_i > -xthr && y.low_i < xthr) {
                lines(c(max(y.low_i/xthr * 2, -2), min(y.upp_i/xthr * 
                  2, 2)), rep(pos.y_i, 2), lwd = lwd[i], col = "#dddddd")
            }
            col <- "#a7a7a7"
        }
        if (y_i > -xthr && y_i < xthr) {
            lines(rep(y_i/xthr * 2, 2), rep(pos.y_i, 2), lwd = lwd[i], 
                col = col)
        }
        if (ci.upp_i > -xthr && ci.low_i < xthr) {
            lines(c(max(ci.low_i/xthr * 2, -2), min(ci.upp_i/xthr * 
                2, 2)), rep(pos.y_i, 2), lend = 2, col = col)
            if (ci.low_i > -xthr) {
                lines(rep(ci.low_i/xthr * 2, 2), pos.y_i + c(-0.15, 
                  0.15), lend = 2, col = col)
            }
            if (ci.upp_i < xthr) {
                lines(rep(ci.upp_i/xthr * 2, 2), pos.y_i + c(-0.15, 
                  0.15), lend = 2, col = col)
            }
        }
        text(-2.1, pos.y_i, labels[i], pos = 2, col = col)
        text(2.1, pos.y_i, ci.text[i], pos = 4, col = col)
    }
    width <- round(diff(xlim))
    height <- round(diff(ylim)/3)
    cat("\n")
    cat("Use pdf(filename, width, height) before calling forest to save it.\n")
    cat("The optimal width and height of this plot is ~", width, 
        " x ~", height, " inches.\n", sep = "")
    cat("\n")
    invisible(list(optimal.width = width, optimal.height = height))
}
funnel <-
function (x, ...) 
UseMethod("funnel")
funnel.meta.nsue <-
function (x, ...) 
{
    if (!inherits(x, "meta.nsue")) {
        .stop(match.call(), "The argument must be a 'meta.nsue' object")
    }
    measure = x$measure
    known <- x$known$i
    known.n.stud <- length(known)
    unknown <- x$unknown$i
    unknown.n.stud <- length(unknown)
    if (nrow(x$rm$M) < known.n.stud + unknown.n.stud) {
        .warning("This analysis does not take repeated measures into account")
    }
    fitted <- fitted(x)
    known.res <- x$known$y - fitted[known]
    if (measure == "cor" || measure == "cor in smd") {
        known.se <- sqrt(x$y.var[known])
    }
    if (measure == "smc" || measure == "smd") {
        known.se <- sqrt(x$y2var_k1[known] + x$y2var_k2[known] * 
            x$known$y^2)
    }
    if (unknown.n.stud) {
        unknown.fitted <- fitted[unknown]
        unknown.res <- apply(x$unknown$y, 1, mean) - unknown.fitted
        if (measure == "cor" || measure == "cor in smd") {
            unknown.se <- sqrt(x$y.var[unknown])
        }
        if (measure == "smc" || measure == "smd") {
            unknown.se <- sqrt(apply(x$y2var_k1[unknown] + x$y2var_k2[unknown] * 
                x$unknown$y^2, 1, mean))
        }
        max.se <- .signif.up(max(c(known.se, unknown.se)), 1)
    }
    else {
        max.se <- .signif.up(max(known.se), 1)
    }
    ci <- qnorm(0.975) * max.se
    plot(NA, NA, type = "n", xlim = 1.3 * c(-ci, ci), ylim = c(max.se, 
        0), lty = 2, frame.plot = FALSE, xlab = "Residual effect size", 
        ylab = "Standard error")
    ci.x <- c(-ci, 0, ci)
    ci.y <- c(max.se, 0, max.se)
    polygon(c(ci.x, rep(1.3 * ci, 2), rep(-1.3 * ci, 2)), c(ci.y, 
        max.se, 0, 0, max.se), col = "#fcfcfc", border = "#dddddd")
    if (unknown.n.stud) {
        for (i in 1:unknown.n.stud) {
            if (measure == "cor" || measure == "cor in smd") {
                lines(c(max(quantile(x$unknown$y[i, ] - unknown.fitted[i], 
                  0.025), -1.3 * ci), min(quantile(x$unknown$y[i, 
                  ] - unknown.fitted[i], 0.975), 1.3 * ci)), 
                  rep(sqrt(x$y.var[unknown][i]), 2), lwd = 7, 
                  col = "#dddddd")
            }
            if (measure == "smc" || measure == "smd") {
                .elliptic.q(x$unknown$y[i, ] - unknown.fitted[i], 
                  sqrt(x$y2var_k1[unknown][i] + x$y2var_k2[unknown][i] * 
                    x$unknown$y[i, ]^2), col = "#dddddd")
            }
        }
    }
    lines(ci.x, ci.y, lty = 2)
    lines(c(0, 0), c(max.se, 0), lty = 2)
    if (unknown.n.stud) {
        for (i in 1:unknown.n.stud) {
            lines(rep(unknown.res[i], 2), rep(unknown.se[i], 
                2), lwd = 7, col = "#a7a7a7")
        }
    }
    for (i in 1:known.n.stud) {
        lines(rep(known.res[i], 2), rep(known.se[i], 2), lwd = 7, 
            col = "#000000")
    }
    cat("\n")
    cat("Use pdf(filename) before calling funnel to save it.\n")
    cat("\n")
}
leave1out <-
function (x, ...) 
UseMethod("leave1out")
leave1out.nsue <-
function (x, formula = ~1, hypothesis = NULL, n.imp = 500, n.bins = 200, 
    maxiter = 200, tol = 1e-06, ...) 
{
    call <- match.call()
    y <- x$y
    measure <- x$measure
    n.stud <- length(y)
    model <- .check.formula(call, formula, n.stud)
    hypothesis <- .check.hypothesis(call, hypothesis, model)
    if (n.imp < 2) {
        .stop(call, "The number of imputations must be at least 2")
    }
    if (length(unique(x$labels)) < n.stud) {
        .warning("This analysis understand repeated measures as separate studies")
    }
    nsue_i <- x
    model_i <- model
    obj <- list()
    for (i in 1:n.stud) {
        nsue_i$y <- x$y[-i]
        nsue_i$y_lo <- x$y_lo[-i]
        nsue_i$y_up <- x$y_up[-i]
        if (measure == "cor" || measure == "cor in smd") {
            nsue_i$y.var <- x$y.var[-i]
        }
        if (measure == "smc" || measure == "smd") {
            nsue_i$y2var_k1 <- x$y2var_k1[-i]
            nsue_i$y2var_k2 <- x$y2var_k2[-i]
        }
        if (measure == "cor in smd") {
            nsue_i$smd = x$smd[-i, ]
        }
        nsue_i$labels <- x$labels[-i]
        class(nsue_i) <- "nsue"
        model_i$matrix <- as.matrix(model$matrix[-i, ])
        obj[[i]] <- list(study = x$labels[i], meta.nsue = .meta.nsue(nsue_i, 
            model_i, hypothesis, n.imp, n.bins, maxiter, tol))
    }
    class(obj) <- "leave1out.nsue"
    obj
}
meta <-
function (x, ...) 
UseMethod("meta")
meta.nsue <-
function (x, formula = ~1, hypothesis = NULL, n.imp = 500, n.bins = 200, 
    maxiter = 200, tol = 1e-06, ...) 
{
    call <- match.call()
    if (!inherits(x, "nsue")) {
        .stop(call, "Use an smc_from_t, smd_from_t or r_from_z call as the first (nsue) argument.")
    }
    n.stud <- length(x$y)
    model <- .check.formula(call, formula, n.stud)
    hypothesis <- .check.hypothesis(call, hypothesis, model)
    if (n.imp < 2) {
        .stop(call, "The number of imputations must be at least 2")
    }
    .meta.nsue(x, model, hypothesis, n.imp, n.bins, maxiter, 
        tol)
}
metabias <-
function (x, ...) 
UseMethod("metabias")
metabias.meta.nsue <-
function (x, maxiter = 100, tol = 1e-06, ...) 
{
    if (!inherits(x, "meta.nsue")) {
        .stop(match.call(), "The argument must be a 'meta.nsue' object")
    }
    if (nrow(x$rm$M) < length(x$known$i) + length(x$unknown$i)) {
        .warning("This analysis does not take repeated measures into account")
    }
    measure <- x$measure
    known <- x$known$i
    unknown <- x$unknown$i
    mi_y = x$unknown$y
    X <- x$model$matrix
    n.coef_j <- ncol(x$model$matrix) + 1
    y <- c()
    y[known] <- x$known$y
    se_coef <- c()
    se_coef_var <- c()
    for (j in 1:max(c(1, ncol(mi_y)))) {
        if (ncol(mi_y) > 0) {
            y[unknown] <- mi_y[, j]
        }
        if (measure == "cor" || measure == "cor in smd") {
            y.var <- x$y.var
        }
        if (measure == "smc" || measure == "smd") {
            y.var <- x$y2var_k1 + x$y2var_k2 * y^2
        }
        X_j <- cbind(X, sqrt(y.var))
        W <- diag(1/(y.var + .tau2.reml(y, y.var, X_j, maxiter, 
            tol)))
        inv_XtWX <- solve(t(X_j) %*% W %*% X_j)
        se_coef <- c(se_coef, (inv_XtWX %*% t(X_j) %*% W %*% 
            y)[n.coef_j])
        se_coef_var <- c(se_coef_var, diag(inv_XtWX)[n.coef_j])
    }
    z <- mean(se_coef)/sqrt(.pool.var(se_coef, se_coef_var))
    names(z) <- "z"
    p <- 2 * pnorm(-abs(z))
    x <- list(method = "'meta.nsue' regression test for funnel plot asymmetry", 
        data.name = as.character(match.call()[2]), statistic = z, 
        p.value = p)
    class(x) <- "htest"
    x
}
plot.meta.nsue <-
function (x, ...) 
{
    if (!inherits(x, "meta.nsue")) {
        .stop(match.call(), "The argument must be a 'meta.nsue' object")
    }
    forest.meta.nsue(x)
}
print.leave1out.nsue <-
function (x, ...) 
{
    if (!inherits(x, "leave1out.nsue")) {
        .stop(match.call(), "The argument must be a 'leave1out.nsue' object")
    }
    cat("\n")
    cat("Meta-analysis description:\n")
    cat("- Measure:", .format.measure(x[[1]]$meta.nsue$measure), 
        "\n")
    cat("- Model: measure", x[[1]]$meta.nsue$model$formula, "\n")
    cat("- Hypothesis: ", paste(x$hypothesis$text, collapse = " & "), 
        "\n")
    cat("\n")
    for (i in 1:length(x)) {
        cat("\n")
        cat("Discarded study:", x[[i]]$study, "\n")
        cat("\n")
        .print.heterogeneity(x[[i]]$meta.nsue)
        cat("\n")
        .print.model(x[[i]]$meta.nsue)
        cat("\n")
        .print.hypothesis(x[[i]]$meta.nsue)
        cat("\n")
    }
    .print.sign()
    cat("\n")
}
print.meta.nsue <-
function (x, ...) 
{
    if (!inherits(x, "meta.nsue")) {
        .stop(match.call(), "The argument must be a 'meta.nsue' object")
    }
    cat("\n")
    cat("Meta-analysis description:\n")
    cat("- Measure:", .format.measure(x$measure), "\n")
    cat("- Known measures:", length(x$known$i), "\n")
    if (length(x$known$i) == 0) {
        .warning("No known measures!")
    }
    cat("- Non-statistically-significant unknown measures:", 
        length(x$unknown$i), "\n")
    rm.n.stud = nrow(x$rm$M)
    if (rm.n.stud < length(x$known$i) + length(x$unknown$i)) {
        cat("- Measures after combination of repeated-measures:", 
            rm.n.stud, "\n")
    }
    cat("- Imputations:", ncol(x$unknown$y), "\n")
    cat("- Model: measure", x$model$formula, "\n")
    cat("- Hypothesis: ", paste(x$hypothesis$text, collapse = " & "), 
        "\n")
    cat("\n")
    .print.heterogeneity(x)
    cat("\n")
    .print.model(x)
    cat("\n")
    .print.hypothesis(x)
    cat("\n")
    .print.sign()
    cat("\n")
}
print.nsue <-
function (x, ...) 
{
    cat("\n")
    cat("'nsue' object description:\n")
    cat("- Measure:", .format.measure(x$measure), "\n")
    known.n.stud = sum(!is.na(x$y))
    unknown.n.stud = sum(is.na(x$y))
    cat("- Known measures:", known.n.stud, "\n")
    if (known.n.stud == 0) {
        .warning("No known measures!")
    }
    cat("- Non-statistically-significant unknown measures:", 
        sum(is.na(x$y)), "\n")
    cat("\n")
}
r_in_smd_from_t_means_and_sds1 <-
function (t, n1, mean1.pre, sd1.pre, mean1.post, sd1.post, n2, 
    mean2.pre, sd2.pre, mean2.post, sd2.post, alpha = 0.05, labels = "study", 
    r.range = c(0, 0.99), rm.r = 0.3) 
{
    call <- match.call()
    if (missing(t) || missing(n1) || missing(n2)) {
        .stop(call, "You must specify t, n1, mean1.pre, sd1.pre, mean1.post, sd1.post, n2, mean2.pre, sd2.pre, mean2.post and sd2.post")
    }
    if (!is.numeric(t)) {
        .stop(call, "t is not a numeric vector")
    }
    n.stud <- length(t)
    if (!n.stud) {
        .stop(call, "No studies to meta-analyze")
    }
    n1 <- .check.n(call, n1, 2, n.stud)
    if (!is.numeric(mean1.pre)) {
        .stop(call, "mean1.pre is not a numeric vector")
    }
    if (!is.numeric(sd1.pre) || any(sd1.pre < 0)) {
        .stop(call, "sd1.pre is not a positive numeric vector")
    }
    if (!is.numeric(mean1.post)) {
        .stop(call, "mean1.post is not a numeric vector")
    }
    if (!is.numeric(sd1.post) || any(sd1.post < 0)) {
        .stop(call, "sd1.post is not a positive numeric vector")
    }
    n2 <- .check.n(call, n2, 2, n.stud)
    if (!is.numeric(mean2.pre)) {
        .stop(call, "mean2.pre is not a numeric vector")
    }
    if (!is.numeric(sd2.pre) || any(sd2.pre < 0)) {
        .stop(call, "sd2.pre is not a positive numeric vector")
    }
    if (!is.numeric(mean2.post)) {
        .stop(call, "mean2.post is not a numeric vector")
    }
    if (!is.numeric(sd2.post) || any(sd2.post < 0)) {
        .stop(call, "sd2.post is not a positive numeric vector")
    }
    alpha <- .check.alpha(call, alpha, n.stud)
    labels <- .check.labels(call, labels, n.stud)
    if (!is.numeric(r.range) || r.range[1] > r.range[2] || r.range[1] < 
        -1 || r.range[2] > 1) {
        .stop(call, "Incorrect r.range")
    }
    if (!is.numeric(rm.r) || rm.r < -1 || rm.r > 1) {
        .stop(call, "Incorrect rm.r")
    }
    for (i in 1:n.stud) {
        if ((is.na(t[i]) && is.na(alpha[i])) || is.na(n1[i]) || 
            is.na(n2[i])) {
            .stop("Not enough information in study", labels[i])
        }
    }
    n <- n1 + n2
    df <- n - 2
    inv_n1_n2 <- 1/n1 + 1/n2
    j <- .d_j(df)
    df1 <- n1 - 1
    df2 <- n2 - 1
    diff <- mean1.post - mean1.pre - mean2.post + mean2.pre
    k_t2var <- diff^2/inv_n1_n2
    var1.sum <- sd1.pre^2 + sd1.post^2
    var2.sum <- sd2.pre^2 + sd2.post^2
    sd1.prod <- sd1.pre * sd1.post
    sd2.prod <- sd2.pre * sd2.post
    r.min <- r.range[1]
    r.max <- r.range[2]
    obj <- list(measure = "cor in smd", y = atanh(.r_in_smd_from_sds(k_t2var/t^2, 
        df1, var1.sum, sd1.prod, df2, var2.sum, sd2.prod, r.min, 
        r.max)), y_lo = atanh(rep(r.min, n.stud)), y_up = atanh(rep(r.max, 
        n.stud)), y.var = (n2/n)^2/(n1 - 3) + (n1/n)^2/(n2 - 
        3), smd = data.frame(diff, df1, var1.sum, sd1.prod, df2, 
        var2.sum, sd2.prod, j = j, y2var_k1 = inv_n1_n2, y2var_k2 = 1 - 
            (df - 2)/(df * j^2)), labels = labels, rm = list(r = rm.r))
    class(obj) <- "nsue"
    obj
}
r_in_smd_from_t_means_and_sds2 <-
function (x, formula = ~1, hypothesis = NULL, maxiter = 200, 
    tol = 1e-06) 
{
    call <- match.call()
    is_meta = inherits(x, "meta.nsue")
    is_leave1out = inherits(x, "leave1out.nsue")
    if (!is_meta && !is_leave1out) {
        .stop(call, "The argument must be a 'meta.nsue' or 'leave1out.nsue' object")
    }
    if (is_leave1out) {
        n.stud = length(x[[1]]$meta.nsue$known$i) + length(x[[1]]$meta.nsue$unknown$i)
    }
    else {
        n.stud = length(x$known$i) + length(x$unknown$i)
    }
    model <- .check.formula(call, formula, n.stud)
    hypothesis <- .check.hypothesis(call, hypothesis, model)
    if (is_leave1out) {
        for (i in 1:length(x)) {
            x[[i]]$meta.nsue = .r_in_smd_from_t_means_and_sds2(x[[i]]$meta.nsue, 
                model, hypothesis, maxiter, tol)
        }
        return(x)
    }
    .r_in_smd_from_t_means_and_sds2(x, model, hypothesis, maxiter, 
        tol)
}
residuals.meta.nsue <-
function (object, ...) 
{
    fitted <- fitted(object)
    known.i <- object$known$i
    residuals <- c()
    residuals[known.i] <- object$known$y - fitted[known.i]
    unknown.i <- object$unknown$i
    if (length(object$unknown$i) > 0) {
        residuals[unknown.i] <- apply(object$unknown$y, 1, mean) - 
            fitted[unknown.i]
    }
    residuals
}
smc_from_t <-
function (t, n, alpha = 0.05, labels = "study", rm.r = 0.3) 
{
    call <- match.call()
    if (missing(t) || missing(n)) {
        .stop(call, "You must specify t and n")
    }
    if (!is.numeric(t)) {
        .stop(call, "t is not a numeric vector")
    }
    n.stud <- length(t)
    if (!n.stud) {
        .stop(call, "No studies to meta-analyze")
    }
    n <- .check.n(call, n, 3, n.stud)
    alpha <- .check.alpha(call, alpha, n.stud)
    labels <- .check.labels(call, labels, n.stud)
    if (!is.numeric(rm.r) || rm.r < -1 || rm.r > 1) {
        .stop(call, "Incorrect rm.r")
    }
    for (i in 1:n.stud) {
        if ((is.na(t[i]) && is.na(alpha[i])) || is.na(n[i])) {
            stop("Not enough information in study", labels[i])
        }
    }
    df <- n - 1
    inv_n <- 1/n
    j <- .d_j(df)
    k_t2d <- j * sqrt(inv_n)
    y_up <- k_t2d * qt(1 - alpha/2, df)
    obj <- list(measure = "smc", y = k_t2d * t, y_lo = -y_up, 
        y_up = y_up, y2var_k1 = inv_n, y2var_k2 = 1 - (df - 2)/(df * 
            j^2), labels = labels, rm = list(r = rm.r))
    class(obj) <- "nsue"
    obj
}
smd_from_t <-
function (t, n1, n2, alpha = 0.05, labels = "study", rm.r = 0.3) 
{
    call <- match.call()
    if (missing(t) || missing(n1) || missing(n2)) {
        .stop(call, "You must specify t, n1 and n2")
    }
    if (!is.numeric(t)) {
        .stop(call, "t is not a numeric vector")
    }
    n.stud <- length(t)
    if (!n.stud) {
        .stop(call, "No studies to meta-analyze")
    }
    n1 <- .check.n(call, n1, 2, n.stud)
    n2 <- .check.n(call, n2, 2, n.stud)
    alpha <- .check.alpha(call, alpha, n.stud)
    labels <- .check.labels(call, labels, n.stud)
    if (!is.numeric(rm.r) || rm.r < -1 || rm.r > 1) {
        .stop(call, "Incorrect rm.r")
    }
    for (i in 1:n.stud) {
        if ((is.na(t[i]) && is.na(alpha[i])) || is.na(n1[i]) || 
            is.na(n2[i])) {
            stop("Not enough information in study", labels[i])
        }
    }
    n <- n1 + n2
    df <- n - 2
    inv_n1_n2 <- 1/n1 + 1/n2
    j <- .d_j(df)
    k_t2d <- j * sqrt(inv_n1_n2)
    y_up = k_t2d * qt(1 - alpha/2, df)
    obj <- list(measure = "smd", y = k_t2d * t, y_lo = -y_up, 
        y_up = y_up, y2var_k1 = inv_n1_n2, y2var_k2 = 1 - (df - 
            2)/(df * j^2), labels = labels, rm = list(r = rm.r))
    class(obj) <- "nsue"
    obj
}
subset.nsue <-
function (x, subset, ...) 
{
    call <- match.call()
    if (!inherits(x, "nsue")) {
        .stop(call, "The argument must be a 'nsue' object")
    }
    if (!is.logical(subset)) {
        .stop(call, "subset must be logical")
    }
    if (length(subset) != length(x$y)) {
        .stop(call, "wrong length")
    }
    measure <- x$measure
    selected = which(subset)
    x$y <- x$y[selected]
    x$y_lo <- x$y_lo[selected]
    x$y_up <- x$y_up[selected]
    if (measure == "cor" || measure == "cor in smd") {
        x$y.var <- x$y.var[selected]
    }
    if (measure == "smc" || measure == "smd") {
        x$y2var_k1 <- x$y2var_k1[selected]
        x$y2var_k2 <- x$y2var_k2[selected]
    }
    x$labels <- x$labels[selected]
    if (measure == "cor in smd") {
        x$smd <- x$smd[selected, ]
    }
    x
}
summary.leave1out.nsue <-
function (object, ...) 
{
    if (!inherits(object, "leave1out.nsue")) {
        .stop(match.call(), "The argument must be a 'leave1out.nsue' object")
    }
    cat("\n")
    cat("Meta-analysis model:", object[[1]]$meta.nsue$measure, 
        object[[1]]$meta.nsue$model$formula, "\n")
    cat("\n")
    for (i in 1:length(object)) {
        cat("Discarded study:", object[[i]]$study, "\n")
        .print.hypothesis(object[[i]]$meta.nsue)
        cat("\n")
    }
    .print.sign()
    cat("\n")
    invisible(object)
}
summary.meta.nsue <-
function (object, ...) 
{
    if (!inherits(object, "meta.nsue")) {
        .stop(match.call(), "The argument must be a 'meta.nsue' object")
    }
    cat("\n")
    cat("Meta-analysis model:", object$measure, object$model$formula, 
        "\n")
    cat("\n")
    .print.hypothesis(object)
    cat("\n")
    .print.sign()
    cat("\n")
    invisible(object)
}
z_from_r <-
function (r, n, alpha = 0.05, labels = "study", rm.r = 0.3) 
{
    call <- match.call()
    if (missing(r) || missing(n)) {
        .stop(call, "You must specify r and n")
    }
    if (!is.numeric(r)) {
        .stop(call, "r is not a numeric vector")
    }
    if (any(r < -1, na.rm = TRUE) || any(r > 1, na.rm = TRUE)) {
        .stop(call, "r cannot be <= -1 or >= 1")
    }
    n.stud <- length(r)
    if (!n.stud) {
        .stop(call, "No studies to meta-analyze")
    }
    n <- .check.n(call, n, 4, n.stud)
    alpha <- .check.alpha(call, alpha, n.stud)
    labels <- .check.labels(call, labels, n.stud)
    if (!is.numeric(rm.r) || rm.r < -1 || rm.r > 1) {
        .stop(call, "Incorrect rm.r")
    }
    for (i in 1:n.stud) {
        if ((is.na(r[i]) && is.na(alpha[i])) || is.na(n[i])) {
            stop("Not enough information in study", labels[i])
        }
    }
    y_up = atanh((1 + (n - 2)/qt(alpha/2, n - 2)^2)^-0.5)
    obj <- list(measure = "cor", y = atanh(r), y_lo = -y_up, 
        y_up = y_up, y.var = 1/(n - 3), labels = labels, rm = list(r = rm.r))
    class(obj) <- "nsue"
    obj
}
