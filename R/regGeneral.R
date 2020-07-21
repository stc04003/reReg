#' General estimating equation when the scale-change model is used to for regression process
#'
#' Returns:
#' alpha: regression coefficietion with length 2 * p; p is the number of covaraites
#' conv: convergence code
#' log.muZ: log of mu_Z
#' zi: Z estimates for id i
#'
#' @importFrom utils tail
#' @noRd
reSC <- function(DF, eqType, solver, a0, wgt = NULL) {
    df0 <- DF[DF$event == 0,]
    df1 <- DF[DF$event == 1,]
    rownames(df0) <- rownames(df1) <- NULL
    m <- aggregate(event ~ id, data = DF, sum)[,2]
    xi <- as.matrix(df1[,-c(1:6)])
    p <- ncol(xi)
    yi <- df0$time2
    yii <- rep(yi, m)
    ti <- df1$time2
    if (is.null(eqType)) {
        ## Used for variance estimation; wgt assume to be a n by p matrix
        texa <- log(ti) + xi %*% a0[1:p]
        yexa <- log(yii) + xi %*% a0[1:p]
        yexa2 <- c(log(yi) + as.matrix(df0[,-c(1:6)]) %*% a0[1:p])
        rate <- apply(wgt, 2, function(e) reRate(texa, yexa, rep(e, m), yexa2))
        rate <- apply(rate, 1, quantile, c(.025, .975))
        Lam <- exp(-rate)
        ind <- !duplicated(yexa2)
        Lam.lower <- approxfun(exp(yexa2)[ind], Lam[2, ind],
                               yleft = min(Lam[2, ]), yright = max(Lam[2, ]))
        Lam.upper <- approxfun(exp(yexa2)[ind], Lam[1, ind],
                               yleft = min(Lam[1, ]), yright = max(Lam[1, ]))
        return(list(Lam0.lower = Lam.lower, Lam0.upper = Lam.upper))
    }
    if (is.null(wgt)) {
        Wi <- rep(1, length(m))
        wi <- rep(Wi, m)
    } else {
        if (length(wgt) != length(m)) stop("Weight length mismatch")
        Wi <- wgt
        wi <- rep(wgt, m)
    }
    if (eqType == "logrank") U1 <- function(a) as.numeric(reLog(a, xi, ti, yii, wi))
    if (eqType == "gehan") U1 <- function(a) as.numeric(reGehan(a, xi, ti, yii, wi))
    Xi <- as.matrix(cbind(1, df0[,-c(1:6)]))
    U2 <- function(b) as.numeric(re2(b, R, Xi, Wi))
    if (is.null(solver)) {
        texa <- log(ti) + xi %*% a0[1:p]
        yexa <- log(yii) + xi %*% a0[1:p]
        yexa2 <- c(log(yi) + as.matrix(df0[,-c(1:6)]) %*% a0[1:p])
        rate <- c(reRate(texa, yexa, wi, yexa2))
        Lam <- exp(-rate)
        R <- m / Lam
        R <- ifelse(R > 1e5, (m + .01) / (Lam + .01), R)
        zi <- R / exp(Xi[,-1, drop = FALSE] %*% tail(a0, p))
        return(list(value = c(U1(a0[1:p]), U2(a0[-(1:p)])), zi = zi))
    } else {
        fit.a <- eqSolve(a0[1:p], U1, solver)
        ahat <- fit.a$par
        texa <- log(ti) + xi %*% ahat
        yexa <- log(yii) + xi %*% ahat
        yexa2 <- c(log(yi) + as.matrix(df0[,-c(1:6)]) %*% ahat)
        rate <- c(reRate(texa, yexa, wi, yexa2))
        Lam <- exp(-rate)
        R <- m / Lam
        R <- ifelse(R > 1e5, (m + .01) / (Lam + .01), R)
        fit.b <- eqSolve(a0[-(1:p)], U2, solver)
        ind <- !duplicated(yexa2)
        return(list(alpha = c(fit.a$par, fit.b$par[-1] + fit.a$par),
                    aconv = c(fit.a$convergence, fit.b$convergence),
                    log.muZ = fit.b$par[1],
                    zi = R / exp(Xi[,-1, drop = FALSE] %*% fit.b$par[-1]),
                    Lam0 = approxfun(exp(yexa2)[ind], Lam[ind],
                                     yleft = min(Lam), yright = max(Lam))))
    }
}

reAR <- function(DF, eqType, solver, a0, wgt = NULL) {
    df0 <- DF[DF$event == 0,]
    df1 <- DF[DF$event == 1,]
    rownames(df0) <- rownames(df1) <- NULL
    m <- aggregate(event ~ id, data = DF, sum)[,2]
    xi <- as.matrix(df1[,-c(1:6)])
    yi <- df0$time2
    yii <- rep(yi, m)
    ti <- df1$time2
    if (is.null(eqType)) {
        texa <- log(ti) + xi %*% a0
        yexa <- log(yii) + xi %*% a0
        yexa2 <- c(log(yi) + as.matrix(df0[,-c(1:6)]) %*% a0)
        rate <- apply(wgt, 2, function(e) reRate(texa, yexa, rep(e, m), yexa2))
        rate <- apply(rate, 1, quantile, c(.025, .975))
        Lam <- exp(-rate)
        ind <- !duplicated(yexa2)
        Lam.lower <- approxfun(exp(yexa2)[ind], Lam[2, ind],
                               yleft = min(Lam[2,]), yright = max(Lam[2,]))
        Lam.upper <- approxfun(exp(yexa2)[ind], Lam[1, ind],
                               yleft = min(Lam[1,]), yright = max(Lam[1,]))
        return(list(Lam0.lower = Lam.lower, Lam0.upper = Lam.upper))
    }
    if (is.null(wgt)) {
        Wi <- rep(1, length(m))
        wi <- rep(Wi, m)
    } else {
        if (length(wgt) != length(m)) stop("Weight length mismatch")
        Wi <- wgt
        wi <- rep(wgt, m)
    }
    if (eqType == "logrank") U1 <- function(a) as.numeric(reLog(a, xi, ti, yii, wi))
    if (eqType == "gehan") U1 <- function(a) as.numeric(reGehan(a, xi, ti, yii, wi))
    if (is.null(solver)) {
        texa <- log(ti) + xi %*% a0
        yexa <- log(yii) + xi %*% a0
        yexa2 <- c(log(yi) + as.matrix(df0[,-c(1:6)]) %*% a0)
        rate <- c(reRate(texa, yexa, wi, yexa2))
        Lam <- exp(-rate)
        R <- m / Lam
        R <- ifelse(R > 1e5, (m + .01) / (Lam + .01), R)
        zi <- Wi * R * exp(as.matrix(df0[,-c(1:6)]) %*% a0)
        return(list(value = U1(a0), zi = zi))
    } else {
        fit.a <- eqSolve(a0, U1, solver)
        ahat <- fit.a$par
        texa <- log(ti) + xi %*% ahat
        yexa <- log(yii) + xi %*% ahat
        yexa2 <- c(log(yi) + as.matrix(df0[,-c(1:6)]) %*% ahat)
        rate <- c(reRate(texa, yexa, wi, yexa2))
        Lam <- exp(-rate)
        R <- m / Lam
        R <- ifelse(R > 1e5, (m + .01) / (Lam + .01), R)
        zi <- R * exp(as.matrix(df0[,-c(1:6)]) %*% fit.a$par)
        return(list(alpha = fit.a$par,
                    aconv = fit.a$convergence,
                    log.muZ = log(mean(zi)), zi = zi,
                    Lam0 = approxfun(exp(yexa2)[!duplicated(yexa2)], Lam[!duplicated(yexa2)],
                                     yleft = min(Lam), yright = max(Lam))))
    }
}

reCox <- function(DF, eqType, solver, a0, wgt = NULL) {
    df0 <- DF[DF$event == 0,]
    df1 <- DF[DF$event == 1,]
    rownames(df0) <- rownames(df1) <- NULL
    m <- aggregate(event ~ id, data = DF, sum)[,2]
    ## yi <- rep(df0$time2, m)
    yi <- df0$time2
    ti <- df1$time2
    if (is.null(eqType)) {
        rate <- apply(wgt, 2, function(e) reRate(ti, rep(yi, m), rep(e, m), yi))
        rate <- apply(rate, 1, quantile, c(.025, .975))
        Lam <- exp(-rate)
        ind <- !duplicated(yi)
        Lam.lower <- approxfun(yi[ind], Lam[2, ind],
                               yleft = min(Lam[2,]), yright = max(Lam[2,]))
        Lam.upper <- approxfun(yi[ind], Lam[1, ind],
                               yleft = min(Lam[1,]), yright = max(Lam[1,]))
        return(list(Lam0.lower = Lam.lower, Lam0.upper = Lam.upper))
    }
    if (is.null(wgt)) {
        Wi <- rep(1, length(m))
        wi <- rep(Wi, m)
    } else {
        if (length(wgt) != length(m)) stop("Weight length mismatch")
        Wi <- wgt
        wi <- rep(wgt, m)
    }
    ## T0 <- sort(unique(c(ti, yi)))
    ## rate <- c(reRate(ti, yi, wi, T0))
    ## yi2 <- as.numeric(df0$time2)
    ## Lam0 <- exp(-rate)
    ## Lam <- Lam0[pmax(1, findInterval(yi2, T0))]
    rate <- c(reRate(ti, rep(yi, m), wi, yi))
    Lam <- exp(-rate)
    R <- m / Lam
    R <- ifelse(R > 1e5, (m + .01) / (Lam + .01), R)
    Xi <- as.matrix(cbind(1, df0[,-c(1:6)]))
    U1 <- function(b) as.numeric(re2(b, R, Xi, Wi))
    if (is.null(solver)) { 
        return(list(value = U1(a0),
                    ## zi = R / exp(Xi[,-1, drop = FALSE] %*% a0[-1])))
                    zi = Wi * R / exp(Xi[,-1, drop = FALSE] %*% a0[-1])))
    } else {
        fit.a <- eqSolve(a0, U1, solver)
        return(list(alpha = fit.a$par[-1],
                    aconv = fit.a$convergence,
                    log.muZ = fit.a$par[1],
                    zi = R / exp(Xi[,-1, drop = FALSE] %*% fit.a$par[-1]),
                    Lam0 = approxfun(yi[!duplicated(yi)], Lam[!duplicated(yi)],
                                     yleft = min(Lam), yright = max(Lam))))
    }
}

reAM <- function(DF, eqType, solver, a0, wgt = NULL) {
    df0 <- DF[DF$event == 0,]
    df1 <- DF[DF$event == 1,]
    rownames(df0) <- rownames(df1) <- NULL
    m <- aggregate(event ~ id, data = DF, sum)[,2]
    xi <- as.matrix(df0[,-c(1:6)])
    yi <- df0$time2
    ti <- df1$time2
    if (is.null(eqType)) {
        texa <- log(ti) + as.matrix(df1[,-c(1:6)]) %*% a0
        yexa <- log(yi) + xi %*% a0
        rate <- apply(wgt, 2, function(e) reRate(texa, rep(yexa, m), rep(e, m), yexa))
        rate <- apply(rate, 1, quantile, c(.025, .975))
        Lam <- exp(-rate)
        ind <- !duplicated(yexa)
        Lam.lower <- approxfun(exp(yexa)[ind], Lam[2, ind],
                               yleft = min(Lam[2,]), yright = max(Lam[2,]))
        Lam.upper <- approxfun(exp(yexa)[ind], Lam[1, ind],
                               yleft = min(Lam[1,]), yright = max(Lam[1,]))
        return(list(Lam0.lower = Lam.lower, Lam0.upper = Lam.upper))
    }
    if (is.null(wgt)) {
        Wi <- rep(1, length(m))
    } else {
        if (length(wgt) != length(m)) stop("Weight length mismatch")
        Wi <- wgt
    }
    U1 <- function(a) as.numeric(am1(a, ti, yi, Wi, xi, m))
    if (is.null(solver)) {
        texa <- log(ti) + as.matrix(df1[,-c(1:6)]) %*% a0
        yexa <- log(yi) + xi %*% a0
        rate <- c(reRate(texa, rep(yexa, m), rep(Wi, m), yexa))
        Lam <- exp(-rate)
        R <- m / Lam
        R <- ifelse(R > 1e5, (m + .01) / (Lam + .01), R)
        return(list(value = U1(a0), zi = R))
    } else {
        fit.a <- eqSolve(a0, U1, solver)
        ahat <- fit.a$par
        texa <- log(ti) + as.matrix(df1[,-c(1:6)]) %*% ahat
        yexa <- log(yi) + xi %*% ahat
        rate <- c(reRate(texa, rep(yexa, m), rep(Wi, m), yexa))
        Lam <- exp(-rate)
        R <- m / Lam
        R <- ifelse(R > 1e5, (m + .01) / (Lam + .01), R)
        return(list(alpha = fit.a$par,
                    aconv = fit.a$convergence,
                    log.muZ = log(mean(R)), zi = R,
                    Lam0 = approxfun(exp(yexa)[!duplicated(yexa)], Lam[!duplicated(yexa)],
                                     yleft = min(Lam), yright = max(Lam))))
    }
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

temSC <- function(DF, eqType, solver, b0, zi, wgt = NULL) {
    df0 <- DF[DF$event == 0,]
    rownames(df0) <- NULL
    xi <- as.matrix(df0[,-c(1:6)])    
    di <- df0$terminal
    yi <- df0$time2
    p <- ncol(xi)
    if (is.null(eqType)) {
        yi <- log(yi) + xi %*% b0[1:p + p]
        yi2 <- sort(unique(yi))
        Haz <- apply(wgt, 2, function(e)
            temHaz(b0[1:p + p], b0[1:p], xi, yi, zi / mean(zi), di, e, yi2))
        Haz <- apply(Haz, 1, quantile, c(.025, .975))
        ind <- !duplicated(exp(yi2))
        Haz.lower <- approxfun(exp(yi2)[ind], Haz[1,ind], yleft = min(Haz[1,]), yright = max(Haz[1,]))
        Haz.upper <- approxfun(exp(yi2)[ind], Haz[2,ind], yleft = min(Haz[2,]), yright = max(Haz[2,]))
        return(list(Haz0.lower = Haz.lower, Haz0.upper = Haz.upper))
    }
    if (is.null(wgt)) {
        wi <- rep(1, nrow(xi))
    } else {
        if (length(wgt) != nrow(xi)) stop("Weight length mismatch")
        wi <- wgt
    }
    if (eqType == "logrank")
        U1 <- function(x) as.numeric(temScLog(x[1:p], x[1:p + p], xi, yi, zi, di, wi))
    if (eqType == "gehan") 
        U1 <- function(x) as.numeric(temScGehan(x[1:p], x[1:p + p], xi, yi, zi, di, wi))
    if (is.null(solver)) return(U1(b0))
    else {
        fit.a <- eqSolve(b0, U1, solver)
        yi <- log(yi) + xi %*% fit.a$par[1:p + p]
        yi2 <- sort(unique(yi))
        ind <- !duplicated(exp(yi2))
        Haz <- c(temHaz(fit.a$par[1:p + p], fit.a$par[1:p], xi, yi, zi / mean(zi), di, wi, yi2))
        return(list(beta = fit.a$par,
                    bconv = fit.a$convergence,
                    Haz0 = approxfun(exp(yi2)[ind], Haz[ind], yleft = min(Haz), yright = max(Haz))))
    }
}

temAM <- function(DF, eqType, solver, b0, zi, wgt = NULL) {
    df0 <- DF[DF$event == 0,]
    rownames(df0) <- NULL
    xi <- as.matrix(df0[,-c(1:6)])    
    di <- df0$terminal
    yi <- df0$time2
    p <- ncol(xi)
    if (is.null(eqType)) {
        yi <- log(yi) + xi %*% b0
        yi2 <- sort(unique(yi))
        Haz <- apply(wgt, 2, function(e) temHaz(b0, b0, xi, yi, zi / mean(zi), di, e, yi2))
        Haz <- apply(Haz, 1, quantile, c(.025, .975))
        Haz.lower <- approxfun(exp(yi2), Haz[1,], yleft = min(Haz[1,]), yright = max(Haz[1,]))
        Haz.upper <- approxfun(exp(yi2), Haz[2,], yleft = min(Haz[2,]), yright = max(Haz[2,]))
        return(list(Haz0.lower = Haz.lower, Haz0.upper = Haz.upper))
    }
    if (is.null(wgt)) {
        wi <- rep(1, nrow(xi))
    } else {
        if (length(wgt) != nrow(xi)) stop("Weight length mismatch")
        wi <- wgt
    }
    if (eqType == "logrank")
        U1 <- function(x) as.numeric(temLog(x, x, xi, yi, zi, di, wi))
    if (eqType == "gehan") 
        U1 <- function(x) as.numeric(temGehan(x, x, xi, yi, zi, di, wi))
    if (is.null(solver)) return(U1(b0))
    else {
        fit.a <- eqSolve(b0, U1, solver)
        yi <- log(yi) + xi %*% fit.a$par
        yi2 <- sort(unique(yi))
        Haz <- c(temHaz(fit.a$par, fit.a$par, xi, yi, zi / mean(zi), di, wi, yi2))
        ind <- !duplicated(exp(yi2))
        return(list(beta = fit.a$par,
                    bconv = fit.a$convergence,
                    Haz0 = approxfun(exp(yi2)[ind], Haz[ind], yleft = min(Haz), yright = max(Haz))))
    }
}

temCox <- function(DF, eqType, solver, b0, zi, wgt = NULL) {
    df0 <- DF[DF$event == 0,]
    rownames(df0) <- NULL
    xi <- as.matrix(df0[,-c(1:6)])    
    di <- df0$terminal
    yi <- df0$time2
    p <- ncol(xi)
    if (is.null(eqType)) {
        yi2 <- sort(unique(yi))
        Haz <- apply(wgt, 2, function(e) temHaz(rep(0, p), b0, xi, yi, zi / mean(zi), di, e, yi2))
        Haz <- apply(Haz, 1, quantile, c(.025, .975))
        Haz.lower <- approxfun(yi2, Haz[1,], yleft = min(Haz[1,]), yright = max(Haz[1,]))
        Haz.upper <- approxfun(yi2, Haz[2,], yleft = min(Haz[2,]), yright = max(Haz[2,]))
        return(list(Haz0.lower = Haz.lower, Haz0.upper = Haz.upper))
    }
    if (is.null(wgt)) {
        wi <- rep(1, nrow(xi))
    } else {
        if (length(wgt) != nrow(xi)) stop("Weight length mismatch")
        wi <- wgt
    }
    if (eqType == "logrank")
        U1 <- function(x) as.numeric(temLog(rep(0, p), x, xi, yi, zi, di, wi))
    if (eqType == "gehan") 
        U1 <- function(x) as.numeric(temGehan(rep(0, p), x, xi, yi, zi, di, wi))
    if (is.null(solver)) return(U1(b0))
    else {
        fit.a <- eqSolve(b0, U1, solver)
        yi2 <- sort(unique(yi))
        Haz <- c(temHaz(rep(0, p), fit.a$par, xi, yi, zi / mean(zi), di, wi, yi2))
        return(list(beta = fit.a$par,
                    bconv = fit.a$convergence,
                    Haz0 = approxfun(sort(unique(yi2)), Haz, yleft = min(Haz), yright = max(Haz))))
    }
}

temAR <- function(DF, eqType, solver, b0, zi, wgt = NULL) {
    df0 <- DF[DF$event == 0,]
    rownames(df0) <- NULL
    xi <- as.matrix(df0[,-c(1:6)])    
    di <- df0$terminal
    yi <- df0$time2
    p <- ncol(xi)
    if (is.null(eqType)) {
        yi <- log(yi) + xi %*% b0
        yi2 <- sort(unique(yi))
        Haz <- apply(wgt, 2, function(e) temHaz(b0, rep(0, p), xi, yi, zi / mean(zi), di, e, yi2))
        Haz <- apply(Haz, 1, quantile, c(.025, .975))
        Haz.lower <- approxfun(exp(yi2), Haz[1,], yleft = min(Haz[1,]), yright = max(Haz[1,]))
        Haz.upper <- approxfun(exp(yi2), Haz[2,], yleft = min(Haz[2,]), yright = max(Haz[2,]))
        return(list(Haz0.lower = Haz.lower, Haz0.upper = Haz.upper))
    }
    if (is.null(wgt)) {
        wi <- rep(1, nrow(xi))
    } else {
        if (length(wgt) != nrow(xi)) stop("Weight length mismatch")
        wi <- wgt
    }
    if (eqType == "logrank")
        U1 <- function(x) as.numeric(temLog(x, rep(0, p), xi, yi, zi, di, wi))
    if (eqType == "gehan") 
        U1 <- function(x) as.numeric(temGehan(x, rep(0, p), xi, yi, zi, di, wi))
    if (is.null(solver)) return(U1(b0))
    else {
        fit.a <- eqSolve(b0, U1, solver)
        yi <- log(yi) + xi %*% fit.a$par
        yi2 <- sort(unique(yi))
        Haz <- c(temHaz(fit.a$par, rep(0, p), xi, yi, zi / mean(zi), di, wi, yi2))
        ind <- !duplicated(exp(yi2))
        list(beta = fit.a$par,
             bconv = fit.a$convergence,
             Haz0 = approxfun(exp(yi2)[ind], Haz[ind], yleft = min(Haz), yright = max(Haz)))
    }
}
