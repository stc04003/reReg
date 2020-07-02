#' General estimating equation when the scale-change model is used to for regression process
#'
#' Returns:
#' alpha: regression coefficietion with length 2 * p; p is the number of covaraites
#' conv: convergence code
#' log.muZ: log of mu_Z
#' zi: Z estimates for id i
#' @noRd
reSC <- function(DF, eqType, solver, a0, b0, wgt = NULL) {
    df0 <- subset(DF, event == 0)
    df1 <- subset(DF, event == 1)
    rownames(df0) <- rownames(df1) <- NULL
    m <- aggregate(event ~ id, data = DF, sum)[,2]
    xi <- as.matrix(df1[,-c(1:6)])
    yi <- rep(df0$time2, m)
    ti <- df1$time2
    if (is.null(wgt)) {
        Wi <- rep(1, length(m))
        wi <- rep(Wi, m)
    } else {
        if (length(wgt) != length(m)) stop("Weight length mismatch")
        Wi <- wgt
        wi <- rep(wgt, m)
    }
    if (eqType == "logrank") U1 <- function(a) as.numeric(relog(a, xi, ti, yi, wi))
    if (eqType == "gehan") U1 <- function(a) as.numeric(regehan(a, xi, ti, yi, wi))
    fit.a <- eqSolve(a0, U1, solver)
    ahat <- fit.a$par
    texa <- log(ti) + xi %*% ahat
    yexa <- log(yi) + xi %*% ahat
    T0 <- sort(unique(c(texa, yexa)))
    rate <- c(reRate(texa, yexa, wi, T0))
    yexa2 <- c(log(df0$time2) + as.matrix(df0[,-c(1:6)]) %*% ahat)
    Lam0 <- exp(-rate)
    Lam <- Lam0[pmax(1, findInterval(yexa2, T0))]
    R <- m / Lam
    R <- ifelse(R > 1e5, (m + .01) / (Lam + .01), R)
    Xi <- as.matrix(cbind(1, df0[,-c(1:6)]))
    U2 <- function(b) as.numeric(re2(b, R, Xi, Wi))
    fit.b <- eqSolve(c(0, b0), U2, solver)
    list(alpha = c(fit.a$par, fit.b$par[-1] + fit.a$par),
         aconv = c(fit.a$convergence, fit.b$convergence),
         log.muZ = fit.b$par[1],
         zi = R / exp(fit.b$par[-1]),
         Lam0 = approxfun(T0, Lam0, yleft = min(Lam0), yright = max(Lam0)))
}

reAR <- function(DF, eqType, solver, a0, b0) {
    df0 <- subset(DF, event == 0)
    df1 <- subset(DF, event == 1)
    rownames(df0) <- rownames(df1) <- NULL
    m <- aggregate(event ~ id, data = DF, sum)[,2]
    xi <- as.matrix(df1[,-c(1:6)])
    yi <- rep(df0$time2, m)
    ti <- df1$time2    
    if (is.null(wgt)) {
        Wi <- rep(1, length(m))
        wi <- rep(Wi, m)
    } else {
        if (length(wgt) != nrow(Xi)) stop("Weight length mismatch")
        Wi <- wgt
        wi <- rep(wgt, m)
    }
    if (eqType == "logrank") U1 <- function(a) as.numeric(relog(a, xi, ti, yi, wi))
    if (eqType == "gehan") U1 <- function(a) as.numeric(regehan(a, xi, ti, yi, wi))
    fit.a <- eqSolve(a0, U1, solver)
    ahat <- fit.a$par
    texa <- log(ti) + xi %*% ahat
    yexa <- log(yi) + xi %*% ahat
    T0 <- sort(unique(c(texa, yexa)))
    rate <- c(reRate(texa, yexa, wi, T0))
    yexa2 <- c(log(df0$time2) + as.matrix(df0[,-c(1:6)]) %*% ahat)
    Lam0 <- exp(-rate)
    Lam <- Lam0[pmax(1, findInterval(yexa2, T0))]
    R <- m / Lam
    R <- ifelse(R > 1e5, (m + .01) / (Lam + .01), R)
    zi <- R / exp(fit.a$par[-1])    
    list(alpha = fit.a$par,
         aconv = fit.a$convergence,
         log.muZ = log(mean(zi)), zi = zi,
         Lam0 = approxfun(T0, Lam0, yleft = min(Lam0), yright = max(Lam0)))
}

reCox <- function(DF, eqType, solver, a0, b0, wgt = NULL) {
    df0 <- subset(DF, event == 0)
    df1 <- subset(DF, event == 1)
    rownames(df0) <- rownames(df1) <- NULL
    m <- aggregate(event ~ id, data = DF, sum)[,2]
    yi <- rep(df0$time2, m)
    ti <- df1$time2
    if (is.null(wgt)) {
        Wi <- rep(1, length(m))
        wi <- rep(Wi, m)
    } else {
        if (length(wgt) != nrow(Xi)) stop("Weight length mismatch")
        Wi <- wgt
        wi <- rep(wgt, m)
    }
    T0 <- sort(unique(c(ti, yi)))
    rate <- c(reRate(ti, yi, wi, T0))
    yi2 <- as.numeric(df0$time2)
    Lam0 <- exp(-rate)
    Lam <- Lam0[pmax(1, findInterval(yi2, T0))]
    R <- m / Lam
    R <- ifelse(R > 1e5, (m + .01) / (Lam + .01), R)
    Xi <- as.matrix(cbind(1, df0[,-c(1:6)]))
    U1 <- function(b) as.numeric(re2(b, R, Xi, Wi))
    fit.a <- eqSolve(c(0, a0), U1, solver)
    list(alpha = fit.a$par[-1],
         aconv = fit.a$convergence,
         log.muZ = fit.a$par[1],
         zi = R / exp(fit.a$par[-1]),
         Lam0 = approxfun(T0, Lam0, yleft = min(Lam0), yright = max(Lam0)))
}

reAM <- function(DF, eqType, solver, a0, b0, wgt = NULL) {
    df0 <- subset(DF, event == 0)
    df1 <- subset(DF, event == 1)
    rownames(df0) <- rownames(df1) <- NULL
    m <- aggregate(event ~ id, data = DF, sum)[,2]
    xi <- as.matrix(df0[,-c(1:6)])
    yi <- df0$time2
    ti <- df1$time2
    if (is.null(wgt)) {
        Wi <- rep(1, length(m))
        wi <- rep(Wi, m)
    } else {
        if (length(wgt) != nrow(Xi)) stop("Weight length mismatch")
        Wi <- wgt
        wi <- rep(wgt, m)
    }
    U1 <- function(a) as.numeric(am1(a, ti, yi, wi, xi, m))
    fit.a <- eqSolve(a0, U1, solver)
    ahat <- fit.a$par
    texa <- log(ti) + as.matrix(df1[,-c(1:6)]) %*% ahat
    yexa <- rep(log(yi) + xi %*% ahat, m)
    T0 <- sort(unique(c(texa, yexa)))
    rate <- c(reRate(texa, yexa, rep(wi, m), T0))
    Lam0 <- exp(-rate)
    Lam <- Lam0[pmax(1, findInterval(yi, T0))]
    R <- m / Lam
    R <- ifelse(R > 1e5, (m + .01) / (Lam + .01), R)
    list(alpha = fit.a$par,
         aconv = fit.a$convergence,
         log.muZ = log(mean(R)), zi = R,
         Lam0 = approxfun(T0, Lam0, yleft = min(Lam0), yright = max(Lam0)))
}

#' General estimating equation when the scale-change model is used for terminal event
#'
#' Returns:
#' beta: regression coefficietion with length p or 2 * p depending on model;
#'       p is the number of covaraites
#' conv: convergence code
#' log.muZ: log of mu_Z
#' zi: Z estimates for id i
#' @noRd

temSC <- function(DF, eqType, solver, a0, b0, zi, wgt = NULL) {
    df0 <- subset(dat, event == 0)
    rownames(df0) <- NULL
    xi <- as.matrix(df0[,-c(1:6)])    
    di <- df0$status
    yi <- df0$Time
    p <- ncol(xi)
    if (is.null(wgt)) {
        wi <- rep(1, nrow(Xi))
    } else {
        if (length(wgt) != nrow(xi)) stop("Weight length mismatch")
        wi <- wgt
    }
    if (eqType == "logrank")
        U1 <- function(x) as.numeric(temScLog(x[1:p], x[1:p + p], xi, yi, zi, di, wi))
    if (eqType == "gehan") 
        U1 <- function(x) as.numeric(temScGehan(x[1:p], x[1:p + p], xi, yi, zi, di, wi))
    fit.a <- eqSolve(c(a0, b0), U1, solver)
    rate <- c(temHaz(fit.a$par[1:p + p], fit.a$par[1:p], xi, yi, zi, di, wi, sort(unique(yi))))
    Lam0 <- exp(-rate)    
    list(beta = fit.a$par,
         bconv = fit.a$convergence,
         Haz0 = approxfun(T0, Lam0, yleft = min(Lam0), yright = max(Lam0)))
}

temAM <- function(DF, eqType, solver, a0, b0, zi, wgt = NULL) {
    df0 <- subset(dat, event == 0)
    rownames(df0) <- NULL
    xi <- as.matrix(df0[,-c(1:6)])    
    di <- df0$status
    yi <- df0$Time
    p <- ncol(xi)
    if (is.null(wgt)) {
        wi <- rep(1, nrow(Xi))
    } else {
        if (length(wgt) != nrow(xi)) stop("Weight length mismatch")
        wi <- wgt
    }
    if (eqType == "logrank")
        U1 <- function(x) as.numeric(temScLog(x, x, xi, yi, zi, di, wi))
    if (eqType == "gehan") 
        U1 <- function(x) as.numeric(temScGehan(x, x, xi, yi, zi, di, wi))
    fit.a <- eqSolve(b0, U1, solver)
    rate <- c(temHaz(fit.a$par, fit.a$par, xi, yi, zi, di, wi, sort(unique(yi))))
    Lam0 <- exp(-rate)    
    list(beta = fit.a$par,
         bconv = fit.a$convergence,
         Haz0 = approxfun(T0, Lam0, yleft = min(Lam0), yright = max(Lam0)))
}

temCox <- function(DF, eqType, solver, a0, b0, zi, wgt = NULL) {
    df0 <- subset(dat, event == 0)
    rownames(df0) <- NULL
    xi <- as.matrix(df0[,-c(1:6)])    
    di <- df0$status
    yi <- df0$Time
    p <- ncol(xi)
    if (is.null(wgt)) {
        wi <- rep(1, nrow(Xi))
    } else {
        if (length(wgt) != nrow(xi)) stop("Weight length mismatch")
        wi <- wgt
    }
    if (eqType == "logrank")
        U1 <- function(x) as.numeric(temScLog(rep(0, p), x, xi, yi, zi, di, wi))
    if (eqType == "gehan") 
        U1 <- function(x) as.numeric(temScGehan(rep(0, p), x, xi, yi, zi, di, wi))
    fit.a <- eqSolve(b0, U1, solver)
    rate <- c(temHaz(rep(0, p), fit.a$par, xi, yi, zi, di, wi, sort(unique(yi))))
    Lam0 <- exp(-rate)    
    list(beta = fit.a$par,
         bconv = fit.a$convergence,
         Haz0 = approxfun(T0, Lam0, yleft = min(Lam0), yright = max(Lam0)))    
}

temAR <- function(DF, eqType, solver, a0, b0, zi, wgt = NULL) {
    df0 <- subset(dat, event == 0)
    rownames(df0) <- NULL
    xi <- as.matrix(df0[,-c(1:6)])    
    di <- df0$status
    yi <- df0$Time
    p <- ncol(xi)
    if (is.null(wgt)) {
        wi <- rep(1, nrow(Xi))
    } else {
        if (length(wgt) != nrow(xi)) stop("Weight length mismatch")
        wi <- wgt
    }
    if (eqType == "logrank")
        U1 <- function(x) as.numeric(temScLog(x, rep(0, p), xi, yi, zi, di, wi))
    if (eqType == "gehan") 
        U1 <- function(x) as.numeric(temScGehan(x, rep(0, p), xi, yi, zi, di, wi))
    fit.a <- eqSolve(b0, U1, solver)
    rate <- c(temHaz(fit.a$par, rep(0, p), xi, yi, zi, di, wi, sort(unique(yi))))
    Lam0 <- exp(-rate)    
    list(beta = fit.a$par,
         bconv = fit.a$convergence,
         Haz0 = approxfun(T0, Lam0, yleft = min(Lam0), yright = max(Lam0)))
}
