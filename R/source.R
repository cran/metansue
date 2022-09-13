.backtransf_identity <-
function (aux, y) 
{
    y
}
.backtransf_tanh <-
function (aux, y) 
{
    tanh(y)
}
.check.alpha <-
function (call, alpha, n.stud) 
{
    if (!is.numeric(alpha)) {
        .stop(call, "alpha must be a numeric vector")
    }
    if (any(is.na(alpha))) {
        .stop(call, "alpha must have no missing values")
    }
    if (any(alpha <= 0) || any(alpha >= 1)) {
        .stop(call, "alpha values must be between 0 and 1")
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
        .stop(call, "Independent variables of the formula must have no missing values. Impute missing values (using for example the R package 'mi'), call 'meta' for each imputation, and combine all imputations")
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
.check.n <-
function (call, n, min.n, n.stud) 
{
    if (!is.numeric(n) || any(n < min.n) || any(is.na(n))) {
        .stop(call, "n must be a numeric vector wit hvalues > min.n and no missing values")
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
    symnum(p, cutpoints = c(0, 0.001, 0.01, 0.050000000000000003, 
        0.10000000000000001, 1), symbols = c("***", "**", "*", 
        ".", " "), na = FALSE, corr = FALSE)
}
.meta.nsue <-
function (x, model, hypothesis, n.imp, maxiter, tol) 
{
    y <- x$y
    y_lo = x$y_lo
    y_up = x$y_up
    x$y = NULL
    x$y_lo = NULL
    x$y_up = NULL
    n.stud <- length(y)
    known <- which(!is.na(y) | y_up - y_lo < 9.9999999999999995e-07)
    unknown <- setdiff(1:n.stud, known)
    y[is.na(y)] = (y_up[is.na(y)] + y_lo[is.na(y)])/2
    aux <- x$aux
    y2var = x$y2var
    y.var = y2var(aux, y)
    y_lo.var = y2var(aux, y_lo)
    y_up.var = y2var(aux, y_up)
    X = model$matrix
    n.coef = ncol(X)
    df <- nrow(X) - n.coef
    if (length(unknown)) {
        if (n.stud == 1) {
            .warning("Only one unknown effect and no known effects!")
            mu = 0
        }
        else {
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
                    if (interval[1] == interval[2]) {
                      coef_i = interval[1]
                    }
                    else {
                      coef_i <- optimize(.mll_coef, interval, 
                        known_i, y[known_i], y.var[known_i], 
                        unknown_i, y_lo[unknown_i], sqrt(y_lo.var[unknown_i]), 
                        y_up[unknown_i], sqrt(y_up.var[unknown_i]), 
                        X)$minimum
                    }
                  }
                  else {
                    initial_coef <- coef(lm.wfit(X[sample_i, 
                      ], y[sample_i], 1/y.var[sample_i]))
                    coef_i <- optim(initial_coef, .mll_coef, 
                      gr = NULL, known_i, y[known_i], y.var[known_i], 
                      unknown_i, y_lo[unknown_i], sqrt(y_lo.var[unknown_i]), 
                      y_up[unknown_i], sqrt(y_up.var[unknown_i]), 
                      X)$par
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
        }
        tau2 <- optimize(.mll_tau2, c(0, 999), (y[known] - mu[known])^2, 
            y.var[known], y_lo[unknown] - mu[unknown], y_lo.var[unknown], 
            y_up[unknown] - mu[unknown], y_up.var[unknown])$minimum
        mi.y <- NULL
        for (i in unknown) {
            mi.y <- rbind(mi.y, x$mi(aux, mu[i], tau2, y_lo[i], 
                y_up[i], y2var, i, n.imp))
        }
        colnames(mi.y) <- NULL
    }
    else {
        mi.y = matrix(nrow = 0, ncol = 0)
    }
    x$known = list(i = known, y = y[known])
    x$unknown = list(i = unknown, y = mi.y)
    hyp <- hypothesis$matrix
    mi.coef <- NULL
    mi.cov <- NULL
    mi.tau <- c()
    mi.i <- c()
    mi.qe <- c()
    for (j in 1:max(c(1, ncol(mi.y)))) {
        if (ncol(mi.y) > 0) {
            y[unknown] <- mi.y[, j]
            y.var[unknown] <- y2var(aux, mi.y[, j], unknown)
        }
        if (n.stud == 1) {
            W_fe_j <- matrix(1/y.var)
        }
        else {
            W_fe_j <- diag(1/y.var)
        }
        P_fe_j <- W_fe_j - W_fe_j %*% X %*% solve(t(X) %*% W_fe_j %*% 
            X) %*% t(X) %*% W_fe_j
        if (n.stud == 1) {
            tau2_j <- 0
            W_j <- as.matrix(1/y.var)
        }
        else {
            tau2_j <- .tau2.reml(y, y.var, X, maxiter, tol)
            W_j <- diag(1/(y.var + tau2_j))
        }
        cov_j <- solve(t(X) %*% W_j %*% X)
        coef_j <- cov_j %*% t(X) %*% W_j %*% y
        if (n.stud == 1) {
            h2_j <- 1
        }
        else {
            h2_j <- 1 + tau2_j/df * sum(diag(P_fe_j))
        }
        i2_j <- 1 - 1/h2_j
        qe_j <- max(0, t(y) %*% P_fe_j %*% y)
        mi.coef <- cbind(mi.coef, coef_j)
        mi.cov <- cbind(mi.cov, c(cov_j))
        mi.tau <- c(mi.tau, sqrt(tau2_j))
        mi.i <- c(mi.i, sqrt(i2_j))
        mi.qe <- c(mi.qe, qe_j)
    }
    coef <- .pool(mi.coef)
    cov <- .pool.cov(mi.coef, mi.cov)
    tau2 <- .pool(mi.tau)^2
    i2 <- .pool(mi.i)^2
    if (length(mi.qe) == 1) {
        qe <- mi.qe
    }
    else {
        if (n.stud == 1) {
            qe <- 0
        }
        else {
            f_df2 <- .pool.chi2(mi.qe, df)
            qe <- qchisq(pf(f_df2[1], df, f_df2[2], log.p = TRUE), 
                df, log.p = TRUE)
        }
    }
    hcoef <- c(hyp %*% coef)
    hcov <- hyp %*% cov %*% t(hyp)
    if (nrow(hyp) == 1) {
        hse <- sqrt(hcov)
        hz <- hcoef/hse
    }
    else {
        hqr <- qr(hcov)
        hpivot <- hqr$pivot[1:hqr$rank]
        hchisq <- c(t(hcoef[hpivot]) %*% solve(hcov[hpivot, hpivot]) %*% 
            hcoef[hpivot])
        hdf <- length(hpivot)
    }
    model$coef <- coef
    model$cov <- cov
    model$se <- sqrt(diag(cov))
    x$model = model
    x$heterogeneity <- list(tau2 = tau2, h2 = 1/(1 - i2), i2 = i2, 
        qe = qe, df = df, p.value = 1 - pchisq(qe, df))
    hypothesis$coef <- hcoef
    if (nrow(hyp) == 1) {
        hypothesis$se <- hse
        hypothesis$z <- hz
        hypothesis$p.value <- 2 * pnorm(-abs(hz))
    }
    else {
        hypothesis$chisq <- hchisq
        hypothesis$df <- hdf
        hypothesis$p.value <- 1 - pchisq(hchisq, df)
    }
    x$hypothesis = hypothesis
    class(x) <- "meta.nsue"
    x
}
.mi_g <-
function (aux, mu, tau2, y_lo, y_up, y2var, selected, n.imp) 
{
    n.bins <- 500
    mus <- rep(mu, n.imp)
    sigma <- sqrt(y2var(aux, mu, selected) + tau2)
    width <- (y_up - y_lo)/n.bins
    y <- sample(seq(y_lo + width/2, y_up - width/2, width))
    q <- c()
    for (imp in 1:n.imp) {
        raw_dens <- dnorm(y, mus[imp], sigma) * (y2var(aux, y, 
            selected) + tau2)
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
.mi_z <-
function (aux, mu, tau2, y_lo, y_up, y2var, selected, n.imp) 
{
    mus <- rep(mu, n.imp)
    sigma <- sqrt(y2var(aux, mu, selected) + tau2)
    q <- rep(NA, n.imp)
    to_imp <- 1:n.imp
    while (length(to_imp)) {
        q[to_imp] <- rnorm(length(to_imp), mus[to_imp], sigma)
        to_imp <- which(q <= y_lo & q >= y_up)
    }
    q
}
.mll_coef <-
function (coef, known, known.y, known.y.var, unknown, unknown.y_lo, 
    unknown.y_lo.se, unknown.y_up, unknown.y_up.se, X) 
{
    mu <- X %*% coef
    unknown.mu <- mu[unknown]
    a <- pnorm((unknown.y_up - unknown.mu)/unknown.y_up.se, log.p = TRUE)
    c <- pnorm((unknown.y_lo - unknown.mu)/unknown.y_up.se, log.p = TRUE)
    sum(log(known.y.var) + (known.y - mu[known])^2/known.y.var)/2 - 
        sum(a + log(-expm1(c - a)))
}
.mll_tau2 <-
function (tau2, known.err2, known.y.var, unknown.err_lo, unknown.y_lo.var, 
    unknown.err_up, unknown.y_up.var) 
{
    if (tau2 < 0) {
        return(Inf)
    }
    known.sigma2 <- known.y.var + tau2
    a <- pnorm(unknown.err_up/sqrt(unknown.y_up.var + tau2), 
        log.p = TRUE)
    c <- pnorm(unknown.err_lo/sqrt(unknown.y_lo.var + tau2), 
        log.p = TRUE)
    sum(log(known.sigma2) + known.err2/known.sigma2)/2 - sum(a + 
        log(-expm1(c - a)))
}
.pool <-
function (x) 
{
    if (is.vector(x)) {
        return(mean(x))
    }
    apply(x, 1, mean)
}
.pool.chi2 <-
function (chi2, df1) 
{
    m <- length(chi2)
    r <- (1 + 1/m) * var(sqrt(chi2))
    c((mean(chi2)/df1 - (m + 1)/(m - 1) * r)/(1 + r), (m - 1)/df1^(3/m) * 
        (1 + 1/r)^2)
}
.pool.cov <-
function (x, x_cov) 
{
    n.coef <- nrow(x)
    n.imp <- ncol(x)
    if (n.imp == 1) {
        return(matrix(x_cov, n.coef))
    }
    matrix(apply(x_cov, 1, mean), n.coef) + (1 + 1/n.imp) * cov(t(x))
}
.pool.var <-
function (x, x_var) 
{
    if (is.vector(x)) {
        n.imp <- length(x)
        if (n.imp == 1) {
            return(x_var)
        }
        return(mean(x_var) + (1 + 1/n.imp) * var(x))
    }
    n.imp <- ncol(x)
    if (n.imp == 1) {
        return(x_var)
    }
    apply(x_var, 1, mean) + (1 + 1/n.imp) * apply(x, 1, var)
}
.print.heterogeneity <-
function (x) 
{
    cat("Residual heterogeneity:", " tau^2:", .format.4pos(x$heterogeneity$tau2), 
        "  I^2:", paste(.format.perc2(x$heterogeneity$i2), "%", 
            sep = ""), "  H^2:", .format.2pos(x$heterogeneity$h2), 
        "\n")
    cat("Q-statistic:", .format.2pos(x$heterogeneity$qe), "on", 
        x$heterogeneity$df, "df  Pr(>Q):", .format.prob(x$heterogeneity$p.value), 
        "\n")
    cat("Note: we strongly suggest focusing more on I^2 than on Pr(>Q)\n")
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
        se <- c(x$hypothesis$se)
        cat("One-row hypothesis:\n")
        hypothesis <- cbind(.format.4(x$backtransf(x$aux, coef)), 
            .format.4(x$hypothesis$z), prob, .format.4(x$backtransf(x$aux, 
                c(coef) + cbind(qnorm(0.025000000000000001), 
                  qnorm(0.97499999999999998)) * se)), sign)
        colnames(hypothesis) <- c("Estimate", "z value", "Pr(>|z|)", 
            "CI(low)", "CI(up)", "")
    }
    else {
        cat("Multi-row hypothesis:\n")
        hypothesis <- cbind(.format.4(x$backtransf(x$aux, coef)), 
            c(.format.2pos(x$hypothesis$chisq), rep("", nrow - 
                1)), c(.format.0pos(x$hypothesis$df), rep("", 
                nrow - 1)), c(prob, rep("", nrow - 1)), c(sign, 
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
.y2var_smc <-
function (aux, y, selected = 1:length(y)) 
{
    1/aux$n[selected] + (1 - (aux$df[selected] - 2)/(aux$df[selected] * 
        aux$j[selected]^2)) * y^2
}
.y2var_smd <-
function (aux, y, selected = 1:length(y)) 
{
    1/aux$n1[selected] + 1/aux$n2[selected] + (1 - (aux$df[selected] - 
        2)/(aux$df[selected] * aux$j[selected]^2)) * y^2
}
.y2var_zr <-
function (aux, y, selected = 1:length(y)) 
{
    if (is.vector(y)) {
        return(1/(aux$n[selected] - 3))
    }
    n.imp <- ncol(y)
    matrix(rep(1/(aux$n[selected] - 3), n.imp), ncol = n.imp)
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
function (x, width, ...) 
UseMethod("forest")
forest.meta.nsue <-
function (x, width = NULL, ...) 
{
    if (!inherits(x, "meta.nsue")) {
        .stop(match.call(), "The argument must be a 'meta.nsue' object")
    }
    known = x$known$i
    unknown = x$unknown$i
    unknown.n.stud <- length(unknown)
    n.stud <- length(known) + unknown.n.stud
    if (length(x$hypothesis$coef) > 1) {
        .warning("This plot only shows the first row of the hypothesis")
    }
    labels <- c(x$labels[unknown], x$labels[known], x$hypothesis$text[1])
    pos.y <- c(n.stud + 2 - c(unknown, known), 0)
    if (unknown.n.stud) {
        y.low <- x$backtransf(x$aux, apply(x$unknown$y, 1, function(x) {
            quantile(x, 0.025000000000000001)
        }))
        y.upp <- x$backtransf(x$aux, apply(x$unknown$y, 1, function(x) {
            quantile(x, 0.97499999999999998)
        }))
        y_untransf <- c(.pool(x$unknown$y), x$known$y, x$hypothesis$coef[1])
        y.se <- c(sqrt(.pool.var(x$unknown$y, x$y2var(x$aux, 
            x$unknown$y, unknown))), sqrt(x$y2var(x$aux, x$known$y, 
            known)), x$hypothesis$se[1])
    }
    else {
        y_untransf <- c(x$known$y, x$hypothesis$coef[1])
        y.se <- c(sqrt(x$y2var(x$aux, x$known$y, known)), x$hypothesis$se[1])
    }
    y <- x$backtransf(x$aux, y_untransf)
    ci.low <- x$backtransf(x$aux, y_untransf + qnorm(0.025000000000000001) * 
        y.se)
    ci.upp <- x$backtransf(x$aux, y_untransf + qnorm(0.97499999999999998) * 
        y.se)
    lwd <- 1/y.se
    lwd <- sqrt(9 + 216 * (lwd - min(lwd))/(max(lwd) - min(lwd)))
    ci.text <- paste0(.format.2(y), " [ ", .format.2(ci.low), 
        ", ", .format.2(ci.upp), " ] ", .format.sign(2 * pnorm(-abs(y_untransf/y.se))))
    if (is.null(width)) {
        width <- .signif.up(max(abs(c(quantile(ci.low, 0.10000000000000001), 
            quantile(ci.upp, 0.90000000000000002)))), 1)
    }
    plot.new()
    xlim <- c(-2.5 - max(strwidth(labels, units = "inches")), 
        max(strwidth(ci.text, units = "inches")) + 2.5)
    ylim <- c(-2, n.stud + 1)
    plot.window(xlim = xlim, ylim = ylim)
    lines(rep(0, 2), c(n.stud + 1.5, -1), col = "#bbbbbb", lty = 1)
    lines(c(-2, 2), rep(-1, 2), col = "#bbbbbb", lty = 1)
    for (pos.x in -2:2) {
        lines(rep(pos.x, 2), c(-1, -1.3), col = "#bbbbbb", lty = 1)
        text(pos.x, -1.5, .format.2(pos.x/2 * width), pos = 1, 
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
            if (y.upp_i > -width && y.low_i < width) {
                lines(c(max(y.low_i/width * 2, -2), min(y.upp_i/width * 
                  2, 2)), rep(pos.y_i, 2), lwd = lwd[i], col = "#dddddd")
            }
            col <- "#a7a7a7"
        }
        if (y_i > -width && y_i < width) {
            lines(rep(y_i/width * 2, 2), rep(pos.y_i, 2), lwd = lwd[i], 
                col = col)
        }
        if (ci.upp_i > -width && ci.low_i < width) {
            lines(c(max(ci.low_i/width * 2, -2), min(ci.upp_i/width * 
                2, 2)), rep(pos.y_i, 2), lend = 2, col = col)
            if (ci.low_i > -width) {
                lines(rep(ci.low_i/width * 2, 2), pos.y_i + c(-0.14999999999999999, 
                  0.14999999999999999), lend = 2, col = col)
            }
            if (ci.upp_i < width) {
                lines(rep(ci.upp_i/width * 2, 2), pos.y_i + c(-0.14999999999999999, 
                  0.14999999999999999), lend = 2, col = col)
            }
        }
        text(-2.1000000000000001, pos.y_i, labels[i], pos = 2, 
            col = col)
        text(2.1000000000000001, pos.y_i, ci.text[i], pos = 4, 
            col = col)
    }
    width <- round(diff(xlim))
    height <- round((diff(ylim) + 5)/2.5)
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
    known <- x$known$i
    known.n.stud <- length(known)
    unknown <- x$unknown$i
    unknown.n.stud <- length(unknown)
    fitted <- fitted(x)
    known.res <- x$known$y - fitted[known]
    known.se <- sqrt(x$y2var(x$aux, x$known$y, known))
    if (unknown.n.stud) {
        n.imp <- ncol(x$unknown$y)
        unknown.res <- x$unknown$y - fitted[unknown]
        unknown.se <- sqrt(x$y2var(x$aux, x$unknown$y, unknown))
        max.se <- .signif.up(max(c(known.se, unknown.se)), 1)
    }
    else {
        max.se <- .signif.up(max(known.se), 1)
    }
    ci <- qnorm(0.97499999999999998) * max.se
    plot(NA, NA, type = "n", xlim = 1.3 * c(-ci, ci), ylim = c(max.se, 
        0), lty = 2, frame.plot = FALSE, xlab = "Residual effect size", 
        ylab = "Standard error")
    ci.x <- c(-ci, 0, ci)
    ci.y <- c(max.se, 0, max.se)
    polygon(c(ci.x, rep(1.3 * ci, 2), rep(-1.3 * ci, 2)), c(ci.y, 
        max.se, 0, 0, max.se), col = "#fcfcfc", border = "#dddddd")
    lines(ci.x, ci.y, lty = 2)
    lines(c(0, 0), c(max.se, 0), lty = 2)
    n.imp <- if (unknown.n.stud) {
        for (i in 1:unknown.n.stud) {
            for (j in round(quantile(order(unknown.res[i, ]), 
                seq(0.01, 0.98999999999999999, 0.01)))) {
                lines(rep(unknown.res[i, j], 2), rep(unknown.se[i, 
                  j], 2), lwd = 21, col = rgb(0, 0, 0, 0.0030000000000000001))
                lines(rep(unknown.res[i, j], 2), rep(unknown.se[i, 
                  j], 2), lwd = 14, col = rgb(0, 0, 0, 0.0050000000000000001))
                lines(rep(unknown.res[i, j], 2), rep(unknown.se[i, 
                  j], 2), lwd = 7, col = rgb(0, 0, 0, 0.01))
            }
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
function (x, formula = ~1, hypothesis = NULL, n.imp = 500, maxiter = 200, 
    tol = 9.9999999999999995e-07, ...) 
{
    call <- match.call()
    y <- x$y
    n.stud <- length(y)
    model <- .check.formula(call, formula, n.stud)
    hypothesis <- .check.hypothesis(call, hypothesis, model)
    if (n.imp < 2) {
        .stop(call, "The number of imputations must be at least 2")
    }
    nsue_i <- x
    model_i <- model
    obj <- list()
    for (i in 1:n.stud) {
        nsue_i$y <- x$y[-i]
        nsue_i$y_lo <- x$y_lo[-i]
        nsue_i$y_up <- x$y_up[-i]
        nsue_i$aux <- x$aux[-i, ]
        nsue_i$labels <- x$labels[-i]
        class(nsue_i) <- "nsue"
        model_i$matrix <- as.matrix(model$matrix[-i, ])
        obj[[i]] <- list(study = x$labels[i], meta = .meta.nsue(nsue_i, 
            model_i, hypothesis, n.imp, maxiter, tol))
    }
    class(obj) <- "leave1out.nsue"
    obj
}
meta <-
function (x, ...) 
UseMethod("meta")
meta.nsue <-
function (x, formula = ~1, hypothesis = NULL, n.imp = 500, maxiter = 200, 
    tol = 9.9999999999999995e-07, ...) 
{
    call <- match.call()
    if (!inherits(x, "nsue")) {
        .stop(call, "Use an nsue, smc_from_t, smd_from_t or r_from_z call as the first (nsue) argument.")
    }
    n.stud <- length(x$y)
    model <- .check.formula(call, formula, n.stud)
    hypothesis <- .check.hypothesis(call, hypothesis, model)
    if (n.imp < 2) {
        .stop(call, "The number of imputations must be at least 2")
    }
    .meta.nsue(x, model, hypothesis, n.imp, maxiter, tol)
}
metabias <-
function (x, ...) 
UseMethod("metabias")
metabias.meta.nsue <-
function (x, maxiter = 100, tol = 9.9999999999999995e-07, ...) 
{
    if (!inherits(x, "meta.nsue")) {
        .stop(match.call(), "The argument must be a 'meta.nsue' object")
    }
    known <- x$known$i
    unknown <- x$unknown$i
    n.stud <- length(known) + length(unknown)
    y <- rep(NA, n.stud)
    y[known] <- x$known$y
    aux <- x$aux
    y2var = x$y2var
    y.var <- y2var(aux, y)
    X <- x$model$matrix
    n.coef <- ncol(x$model$matrix) + 1
    mi.y <- x$unknown$y
    mi.hcoef <- c()
    mi.hvar <- c()
    for (j in 1:max(c(1, ncol(mi.y)))) {
        if (ncol(mi.y) > 0) {
            y[unknown] <- mi.y[, j]
            y.var[unknown] <- y2var(aux, mi.y[, j], unknown)
        }
        X_j <- cbind(X, sqrt(y.var))
        tau2_j <- .tau2.reml(y, y.var, X_j, maxiter, tol)
        W_j <- diag(1/(y.var + tau2_j))
        cov_j <- solve(t(X_j) %*% W_j %*% X_j)
        coef_j <- cov_j %*% t(X_j) %*% W_j %*% y
        mi.hcoef <- c(mi.hcoef, coef_j[n.coef])
        mi.hvar <- c(mi.hvar, cov_j[n.coef, n.coef])
    }
    hcoef <- .pool(mi.hcoef)
    hvar <- .pool.var(mi.hcoef, mi.hvar)
    hz <- hcoef/sqrt(hvar)
    names(hz) <- "z"
    x <- list(method = "'meta.nsue' regression test for funnel plot asymmetry", 
        data.name = as.character(match.call()[2]), statistic = hz, 
        p.value = 2 * pnorm(-abs(hz)))
    class(x) <- "htest"
    x
}
nsue <-
function (y, y_lo = -y_up, y_up, aux, y2var, mi, backtransf = .backtransf_identity, 
    measure = "effect size", labels = "study") 
{
    if (!is.numeric(y)) {
        .stop(call, "y must be a numeric vector")
    }
    n.stud <- length(y)
    if (!is.numeric(y_lo)) {
        .stop(call, "y_lo must be a numeric vector")
    }
    if (length(y_lo) != n.stud) {
        .stop(call, "y_lo has an incorrect length")
    }
    if (!is.numeric(y_up)) {
        .stop(call, "y_up must be a numeric vector")
    }
    if (length(y_up) != n.stud) {
        .stop(call, "y_up has an incorrect length")
    }
    if (!is.data.frame(aux)) {
        .stop(call, "aux must be a data.frame")
    }
    if (nrow(aux) != n.stud) {
        .stop(call, "aux has an incorrect number of rows")
    }
    if (!is.function(y2var)) {
        .stop(call, "y2var must be a function")
    }
    if (!is.function(mi)) {
        .stop(call, "mi must be a function")
    }
    if (!is.function(backtransf)) {
        .stop(call, "backtransf must be a function")
    }
    if (!is.numeric(y)) {
        .stop(call, "y must be a numeric vector")
    }
    if (!is.vector(measure) && !is.factor(measure) && length(measure) != 
        1) {
        .stop(call, "measure must be a vector of length 1")
    }
    if (!is.vector(labels) && !is.factor(labels)) {
        .stop(call, "labels must be a vector")
    }
    if (length(labels) == 1) {
        labels = paste0(labels, 1:n.stud)
    }
    else if (length(labels) != length(unique(labels))) {
        .stop(call, "labels must be unique")
    }
    else if (length(labels) != n.stud) {
        .stop(call, "labels has an incorrect length")
    }
    obj <- list(y = y, y_lo = y_lo, y_up = y_up, aux = aux, y2var = y2var, 
        mi = mi, backtransf = backtransf, measure = as.character(measure), 
        labels = as.character(labels))
    class(obj) <- "nsue"
    obj
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
    cat("- Measure:", x[[1]]$meta$measure, "\n")
    cat("- Model: measure", x[[1]]$meta$model$formula, "\n")
    cat("- Hypothesis: ", paste(x$hypothesis$text, collapse = " & "), 
        "\n")
    cat("\n")
    for (i in 1:length(x)) {
        cat("\n")
        cat("Discarded study:", x[[i]]$study, "\n")
        cat("\n")
        .print.heterogeneity(x[[i]]$meta)
        cat("\n")
        .print.model(x[[i]]$meta)
        cat("\n")
        .print.hypothesis(x[[i]]$meta)
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
    cat("- Measure:", x$measure, "\n")
    cat("- Known effects:", length(x$known$i), "\n")
    if (length(x$known$i) == 0) {
        .warning("No known effects!")
    }
    cat("- Non-statistically significant unknown effects:", length(x$unknown$i), 
        "\n")
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
    cat("- Measure:", x$measure, "\n")
    known.n.stud = sum(!is.na(x$y))
    unknown.n.stud = sum(is.na(x$y))
    cat("- Known effects:", known.n.stud, "\n")
    if (known.n.stud == 0) {
        .warning("No known effects!")
    }
    cat("- Non-statistically significant unknown effects:", sum(is.na(x$y)), 
        "\n")
    cat("\n")
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
        residuals[unknown.i] <- .pool(object$unknown$y) - fitted[unknown.i]
    }
    residuals
}
smc_from_t <-
function (t, n, alpha = 0.050000000000000003, labels = "study") 
{
    call <- match.call()
    if (missing(t) || missing(n)) {
        .stop(call, "You must specify t and n")
    }
    if (!is.numeric(t)) {
        .stop(call, "t must be a numeric vector")
    }
    n.stud <- length(t)
    if (!n.stud) {
        .stop(call, "No studies to meta-analyze")
    }
    n <- .check.n(call, n, 3, n.stud)
    alpha <- .check.alpha(call, alpha, n.stud)
    for (i in 1:n.stud) {
        if ((is.na(t[i]) && is.na(alpha[i])) || is.na(n[i])) {
            stop("Not enough information in study", labels[i])
        }
    }
    df <- n - 1
    j <- .d_j(df)
    nsue(y = j * sqrt(1/n) * t, y_up = j * sqrt(1/n) * qt(1 - 
        alpha/2, df), aux = data.frame(n, df, j), y2var = .y2var_smc, 
        mi = .mi_g, measure = "Standardized mean change (Hedges corrected)", 
        labels = labels)
}
smd_from_t <-
function (t, n1, n2, alpha = 0.050000000000000003, labels = "study") 
{
    call <- match.call()
    if (missing(t) || missing(n1) || missing(n2)) {
        .stop(call, "You must specify t, n1 and n2")
    }
    if (!is.numeric(t)) {
        .stop(call, "t must be a numeric vector")
    }
    n.stud <- length(t)
    if (!n.stud) {
        .stop(call, "No studies to meta-analyze")
    }
    n1 <- .check.n(call, n1, 2, n.stud)
    n2 <- .check.n(call, n2, 2, n.stud)
    alpha <- .check.alpha(call, alpha, n.stud)
    for (i in 1:n.stud) {
        if ((is.na(t[i]) && is.na(alpha[i])) || is.na(n1[i]) || 
            is.na(n2[i])) {
            stop("Not enough information in study", labels[i])
        }
    }
    df <- n1 + n2 - 2
    j <- .d_j(df)
    nsue(y = j * sqrt(1/n1 + 1/n2) * t, y_up = j * sqrt(1/n1 + 
        1/n2) * qt(1 - alpha/2, df), aux = data.frame(n1, n2, 
        df, j), y2var = .y2var_smd, mi = .mi_g, measure = "Standardized mean difference (Hedges corrected)", 
        labels = labels)
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
    selected <- which(subset)
    x$y <- x$y[selected]
    x$y_lo <- x$y_lo[selected]
    x$y_up <- x$y_up[selected]
    x$aux <- x$aux[selected, ]
    x$labels <- x$labels[selected]
    x
}
summary.leave1out.nsue <-
function (object, ...) 
{
    if (!inherits(object, "leave1out.nsue")) {
        .stop(match.call(), "The argument must be a 'leave1out.nsue' object")
    }
    cat("\n")
    cat("Meta-analysis model:", object[[1]]$meta$measure, object[[1]]$meta$model$formula, 
        "\n")
    cat("\n")
    for (i in 1:length(object)) {
        cat("Discarded study:", object[[i]]$study, "\n")
        .print.hypothesis(object[[i]]$meta)
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
zcor_from_r <-
function (r, n, alpha = 0.050000000000000003, labels = "study") 
{
    call <- match.call()
    if (missing(r) || missing(n)) {
        .stop(call, "You must specify r and n")
    }
    if (!is.numeric(r) || any(r < -1, na.rm = TRUE) || any(r > 
        1, na.rm = TRUE)) {
        .stop(call, "r must be a numeric vector with values between -1 and 1")
    }
    n.stud <- length(r)
    if (!n.stud) {
        .stop(call, "No studies to meta-analyze")
    }
    n <- .check.n(call, n, 4, n.stud)
    alpha <- .check.alpha(call, alpha, n.stud)
    for (i in 1:n.stud) {
        if ((is.na(r[i]) && is.na(alpha[i])) || is.na(n[i])) {
            stop("Not enough information in study", labels[i])
        }
    }
    nsue(y = atanh(r), y_up = atanh(1/sqrt(1 + (n - 2)/qt(alpha/2, 
        n - 2)^2)), y2var = .y2var_zr, aux = data.frame(n), mi = .mi_z, 
        backtransf = .backtransf_tanh, measure = "Pearson correlation coefficient (using Fisher's transform)", 
        labels = labels)
}
