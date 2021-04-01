#' General estimating equation when the scale-change model is used to for regression process
#'
#' Returns:
#' alpha: regression coefficietion with length 2 * p; p is the number of covaraites
#' conv: convergence code
#' log.muZ: log of mu_Z
#' zi: Z estimates for id i
#'
#' @param DF is the data.frame from reReg
#' @param eqType is either logrank or gehan; if eqType = NULL, then do non-parametric
#' @param solver is the solver name; if solver = NULL, evaluate the estimating equations for
#' sandwich estimator
#' @param par1 is \alpha from Xu et al. (2019)
#' @param par2 is \theta from Xu et al. (2019)
#' @importFrom utils tail
#' @noRd
reSC <- function(DF, eqType, solver, par1, par2, Lam0 = NULL, w1 = NULL, w2 = NULL) {
    df0 <- DF[DF$event == 0,]
    df1 <- DF[DF$event > 0,]
    rownames(df0) <- rownames(df1) <- NULL
    m <- aggregate(event > 0 ~ id, data = DF, sum)[,2]
    xi <- as.matrix(df1[,-c(1:6)])
    p <- ncol(xi)
    yi <- df0$time2
    yii <- rep(yi, m)
    ti <- df1$time2
    if (is.null(eqType)) {
        ## Used for variance estimation; wgt assume to be a n by p matrix
        texa <- log(ti) + xi %*% par1
        yexa <- log(yii) + xi %*% par1
        yexa2 <- c(log(yi) + as.matrix(df0[,-c(1:6)]) %*% par1)
        rate <- apply(w1, 2, function(e) reRate(texa, yexa, rep(e, m), yexa2))
        rate <- apply(rate, 1, quantile, c(.025, .975))
        Lam <- exp(-rate)
        ind <- !duplicated(yexa2)
        Lam.lower <- approxfun(exp(yexa2)[ind], Lam[2, ind],
                               yleft = min(Lam[2, ]), yright = max(Lam[2, ]))
        Lam.upper <- approxfun(exp(yexa2)[ind], Lam[1, ind],
                               yleft = min(Lam[1, ]), yright = max(Lam[1, ]))
        return(list(Lam0.lower = Lam.lower, Lam0.upper = Lam.upper))
    }
    if (is.null(w1)) w1 <- rep(1, length(m))
    if (is.null(w2)) w2 <- rep(1, length(m))
    if (eqType == "logrank") U1 <- function(a) as.numeric(reLog(a, xi, ti, yii, rep(w1, m)))
    if (eqType == "gehan") U1 <- function(a) as.numeric(reGehan(a, xi, ti, yii, rep(w1, m)))
    Xi <- as.matrix(cbind(1, df0[,-c(1:6)]))
    U2 <- function(b) as.numeric(re2(b, R, Xi, w1))
    if (is.null(solver)) {
        texa <- log(ti) + xi %*% par1
        yexa <- log(yii) + xi %*% par1
        yexa2 <- c(log(yi) + as.matrix(df0[,-c(1:6)]) %*% par1)
        if (is.null(Lam0)) {
            rate <- c(reRate(texa, yexa, rep(w1, m), yexa2))
            Lam <- exp(-rate)
        } else {
            Lam <- Lam0(exp(yexa2))
        }
        R <- w2 * (m + 0.01) / (Lam + 0.01)
        zi <- R / exp(Xi[,-1, drop = FALSE] %*% par2[-1])
        return(list(value = c(U1(par1), re2(par2, R, Xi, rep(1, length(m)))), zi = zi))
    } else {
        fit.a <- eqSolve(par1, U1, solver)
        ahat <- fit.a$par
        texa <- log(ti) + xi %*% ahat
        yexa <- log(yii) + xi %*% ahat
        yexa2 <- c(log(yi) + as.matrix(df0[,-c(1:6)]) %*% ahat)
        rate <- c(reRate(texa, yexa, rep(w1, m), yexa2))
        Lam <- exp(-rate)
        R <- (m + 0.01) / (Lam + 0.01)
        ## R <- m / Lam
        ## R <- ifelse(R > 1e5, (m + .01) / (Lam + .01), R)
        fit.b <- eqSolve(par2, U2, solver)
        ind <- !duplicated(yexa2)
        return(list(par1 = fit.a$par,
                    par2 = fit.b$par,
                    par1.conv = fit.a$convergence,
                    par2.conv = fit.b$convergence,
                    ## alpha = c(fit.a$par, fit.b$par[-1] + fit.a$par),
                    ## aconv = c(fit.a$convergence, fit.b$convergence),
                    log.muZ = fit.b$par[1],
                    zi = R / exp(Xi[,-1, drop = FALSE] %*% fit.b$par[-1]),
                    Lam0 = function(x)
                        approx(x = yexa2[ind], y = Lam[ind], xout = log(x),
                               yleft = min(Lam), yright = max(Lam))$y))
                    ## Lam0 = approxfun(exp(yexa2)[ind], Lam[ind],
                    ##                  yleft = min(Lam), yright = max(Lam))))
    }
}

reAR <- function(DF, eqType, solver, par1, Lam0 = NULL, w1 = NULL) {
    df0 <- DF[DF$event == 0,]
    df1 <- DF[DF$event > 0,]
    rownames(df0) <- rownames(df1) <- NULL
    m <- aggregate(event > 0 ~ id, data = DF, sum)[,2]
    xi <- as.matrix(df1[,-c(1:6)])
    yi <- df0$time2
    yii <- rep(yi, m)
    ti <- df1$time2
    if (is.null(eqType)) {
        texa <- log(ti) + xi %*% par1
        yexa <- log(yii) + xi %*% par1
        yexa2 <- c(log(yi) + as.matrix(df0[,-c(1:6)]) %*% par1)
        rate <- apply(w1, 2, function(e) reRate(texa, yexa, rep(e, m), yexa2))
        rate <- apply(rate, 1, quantile, c(.025, .975))
        Lam <- exp(-rate)
        ind <- !duplicated(yexa2)
        Lam.lower <- approxfun(exp(yexa2)[ind], Lam[2, ind],
                               yleft = min(Lam[2,]), yright = max(Lam[2,]))
        Lam.upper <- approxfun(exp(yexa2)[ind], Lam[1, ind],
                               yleft = min(Lam[1,]), yright = max(Lam[1,]))
        return(list(Lam0.lower = Lam.lower, Lam0.upper = Lam.upper))
    }
    if (is.null(w1)) w1 <- rep(1, length(m))
    if (eqType == "logrank") U1 <- function(a) as.numeric(reLog(a, xi, ti, yii, rep(w1, m)))
    if (eqType == "gehan") U1 <- function(a) as.numeric(reGehan(a, xi, ti, yii, rep(w1, m)))
    if (is.null(solver)) {
        texa <- log(ti) + xi %*% par1
        yexa2 <- c(log(yi) + as.matrix(df0[,-c(1:6)]) %*% par1)
        if (is.null(Lam0)) {
            yexa <- log(yii) + xi %*% par1
            rate <- c(reRate(texa, yexa, rep(w1, m), yexa2))
            Lam <- exp(-rate)
        } else {
            Lam <- Lam0(exp(yexa2))
        }
        R <- (m + 0.01) / (Lam + 0.01)
        zi <- R * exp(as.matrix(df0[,-c(1:6)]) %*% par1)
        return(list(value = U1(par1) / length(m), zi = zi))
    } else {
        fit.a <- eqSolve(par1, U1, solver)
        ahat <- fit.a$par
        texa <- log(ti) + xi %*% ahat
        yexa <- log(yii) + xi %*% ahat
        yexa2 <- c(log(yi) + as.matrix(df0[,-c(1:6)]) %*% ahat)
        rate <- c(reRate(texa, yexa, rep(w1, m), yexa2))
        Lam <- exp(-rate)
        R <- (m + 0.01) / (Lam + 0.01)
        ## R <- m / Lam
        ## R <- ifelse(R > 1e5, (m + .01) / (Lam + .01), R)
        zi <- R * exp(as.matrix(df0[,-c(1:6)]) %*% fit.a$par)
        return(list(par1 = fit.a$par,
                    par1.conv = fit.a$convergence,
                    log.muZ = log(mean(zi)), zi = zi,
                    Lam0 = function(x)
                        approx(x = yexa2[!duplicated(yexa2)], y = Lam[!duplicated(yexa2)],
                               xout = log(x),
                               yleft = min(Lam), yright = max(Lam))$y))
    ## approxfun(exp(yexa2)[!duplicated(yexa2)], Lam[!duplicated(yexa2)],
    ##                                  yleft = min(Lam), yright = max(Lam))))
    }
}

#' @param par1 is \gamma in Huang et al (2004)
#' @param eqType is either logrank or gehan; if Null this only returns baseline estimate 
#' @param solver specifies the solver; if NULL evaluate the estimating equation instead of solving it
#' @param Lam0 is the estiamted cumulative baseline rate function; if NULL calculate here
#' 
#' @noRd
reCox <- function(DF, eqType, solver, par1, Lam0 = NULL, w1 = NULL) {
    df0 <- DF[DF$event == 0,]
    df1 <- DF[DF$event > 0,]
    rownames(df0) <- rownames(df1) <- NULL
    m <- aggregate(event > 0 ~ id, data = DF, sum)[,2]
    ## yi <- rep(df0$time2, m)
    yi <- df0$time2
    ti <- df1$time2
    t0 <- sort(unique(c(ti, yi)))
    if (is.null(eqType)) {
        rate <- apply(w1, 2, function(e) reRate(ti, rep(yi, m), rep(e, m), t0))
        rate <- apply(rate, 1, quantile, c(.025, .975))
        Lam <- exp(-rate)
        Lam.lower <- approxfun(t0, Lam[2,], yleft = min(Lam[2,]), yright = max(Lam[2,]))
        Lam.upper <- approxfun(t0, Lam[1,], yleft = min(Lam[1,]), yright = max(Lam[1,]))
        return(list(Lam0.lower = Lam.lower, Lam0.upper = Lam.upper))
    }
    if (is.null(w1)) w1 <- rep(1, length(m))
    if (is.null(Lam0)) {
        rate <- c(reRate(ti, rep(yi, m), rep(w1, m), t0))
        Lam0 <- exp(-rate)
        Lam <- Lam0[findInterval(yi, t0)]
    } else {
        Lam <- Lam0(yi)
    }
    ## R <- (m + 0.01) / (Lam + 0.01)
    R <- m / Lam
    Xi <- as.matrix(cbind(1, df0[,-c(1:6)]))
    U1 <- function(b) as.numeric(re2(b, R, Xi, w1))
    if (is.null(solver)) { 
        return(list(value = U1(par1),
                    zi = R / exp(Xi[,-1, drop = FALSE] %*% par1[-1])))
    } else {
        fit.a <- eqSolve(par1, U1, solver)
        return(list(par1 = fit.a$par, ## alpha = fit.a$par[-1],
                    par1.conv = fit.a$convergence,
                    log.muZ = fit.a$par[1],
                    zi = R / exp(Xi[,-1, drop = FALSE] %*% fit.a$par[-1]),
                    Lam0 = approxfun(t0, Lam0, yleft = min(Lam0), yright = max(Lam0))))
    }
}

reAM <- function(DF, eqType, solver, par1, Lam0 = NULL, w1 = NULL) {
    df0 <- DF[DF$event == 0,]
    df1 <- DF[DF$event > 0,]
    rownames(df0) <- rownames(df1) <- NULL
    m <- aggregate(event > 0 ~ id, data = DF, sum)[,2]
    xi <- as.matrix(df0[,-c(1:6)])
    yi <- df0$time2
    ti <- df1$time2
    if (is.null(eqType)) {
        texa <- log(ti) + as.matrix(df1[,-c(1:6)]) %*% par1
        yexa <- log(yi) + xi %*% par1
        rate <- apply(w1, 2, function(e) reRate(texa, rep(yexa, m), rep(e, m), yexa))
        rate <- apply(rate, 1, quantile, c(.025, .975))
        Lam <- exp(-rate)
        ind <- !duplicated(yexa)
        Lam.lower <- approxfun(exp(yexa)[ind], Lam[2, ind],
                               yleft = min(Lam[2,]), yright = max(Lam[2,]))
        Lam.upper <- approxfun(exp(yexa)[ind], Lam[1, ind],
                               yleft = min(Lam[1,]), yright = max(Lam[1,]))
        return(list(Lam0.lower = Lam.lower, Lam0.upper = Lam.upper))
    }
    if (is.null(w1)) w1 <- rep(1, length(m))
    U1 <- function(a) as.numeric(am1(a, ti, yi, w1, xi, m))
    if (is.null(solver)) {
        texa <- log(ti) + as.matrix(df1[,-c(1:6)]) %*% par1
        yexa <- log(yi) + xi %*% par1
        if (is.null(Lam0)) {
            rate <- c(reRate(texa, rep(yexa, m), rep(w1, m), yexa))
            Lam <- exp(-rate)
        } else {
            Lam <- Lam0(exp(yexa))
        }
        R <- (m + 0.01) / (Lam + 0.01)
        return(list(value = as.numeric(U1(par1)), zi = R))
    } else {
        fit.a <- eqSolve(par1, U1, solver)
        ahat <- fit.a$par
        texa <- log(ti) + as.matrix(df1[,-c(1:6)]) %*% ahat
        yexa <- log(yi) + xi %*% ahat
        rate <- c(reRate(texa, rep(yexa, m), rep(w1, m), yexa))
        Lam <- exp(-rate)
        R <- (m + 0.01) / (Lam + 0.01)
        return(list(par1 = fit.a$par,
                    par1.conv = fit.a$convergence,
                    log.muZ = log(mean(R)), zi = R,
                    Lam0 = function(x)
                        approx(x = yexa[!duplicated(yexa)], y = Lam[!duplicated(yexa)],
                               xout = log(x), yleft = min(Lam), yright = max(Lam))$y))
        ## Lam0 = approxfun(exp(yexa)[!duplicated(yexa)], Lam[!duplicated(yexa)],
        ##                  yleft = min(Lam), yright = max(Lam))))
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

#' @param par3 is \eta
#' @param par4 is \theta
temSC <- function(DF, eqType, solver, par3, par4, zi, wgt = NULL) {
    df0 <- DF[DF$event == 0,]
    rownames(df0) <- NULL
    xi <- as.matrix(df0[,-c(1:6)])    
    di <- df0$terminal
    yi <- df0$time2
    p <- ncol(xi)
    if (is.null(eqType)) {
        yi <- log(yi) + xi %*% par4
        yi2 <- sort(unique(yi))
        Haz <- apply(wgt, 2, function(e)
            temHaz(par3, par4, xi, yi, zi / mean(zi), di, e, yi2))
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
    if (is.null(solver)) return(U1(c(par3, par4)))
    else {
        fit.a <- eqSolve(c(par3, par4), U1, solver)
        yi <- log(yi) + xi %*% fit.a$par[1:p + p]
        yi2 <- sort(unique(yi))
        ind <- !duplicated(exp(yi2))
        Haz <- c(temHaz(fit.a$par[1:p + p], fit.a$par[1:p], xi, yi, zi / mean(zi), di, wi, yi2))
        return(list(par3 = fit.a$par[1:p],
                    par4 = fit.a$par[1:p + p],
                    par3.conv = fit.a$convergence,
                    par4.conv = fit.a$convergence,
                    Haz0 = approxfun(exp(yi2)[ind], Haz[ind], yleft = min(Haz), yright = max(Haz))))
    }
}

temAM <- function(DF, eqType, solver, par3, zi, wgt = NULL) {
    df0 <- DF[DF$event == 0,]
    rownames(df0) <- NULL
    xi <- as.matrix(df0[,-c(1:6)])    
    di <- df0$terminal
    yi <- df0$time2
    p <- ncol(xi)
    if (is.null(eqType)) {
        yi <- log(yi) + xi %*% par3
        yi2 <- sort(unique(yi))
        Haz <- apply(wgt, 2, function(e) temHaz(par3, par3, xi, yi, zi / mean(zi), di, e, yi2))
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
    if (is.null(solver)) return(U1(par3))
    else {
        fit.a <- eqSolve(par3, U1, solver)
        yi <- log(yi) + xi %*% fit.a$par
        yi2 <- sort(unique(yi))
        Haz <- c(temHaz(fit.a$par, fit.a$par, xi, yi, zi / mean(zi), di, wi, yi2))
        ind <- !duplicated(exp(yi2))
        return(list(par3 = fit.a$par,
                    par3.conv = fit.a$convergence,
                    Haz0 = approxfun(exp(yi2)[ind], Haz[ind], yleft = min(Haz), yright = max(Haz))))
    }
}

temCox <- function(DF, eqType, solver, par3, zi, wgt = NULL) {
    df0 <- DF[DF$event == 0,]
    rownames(df0) <- NULL
    xi <- as.matrix(df0[,-c(1:6)])    
    di <- df0$terminal
    yi <- df0$time2
    p <- ncol(xi)
    ## zi <- zi / mean(zi)
    if (is.null(eqType)) {
        yi2 <- sort(unique(yi))
        Haz <- apply(wgt, 2, function(e) temHaz(rep(0, p), par3, xi, yi, zi / mean(zi), di, e, yi2))
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
    if (is.null(solver)) return(U1(par3))
    else {
        fit.a <- eqSolve(par3, U1, solver)
        yi2 <- sort(unique(yi))
        Haz <- c(temHaz(rep(0, p), fit.a$par, xi, yi, zi / mean(zi), di, wi, yi2))
        return(list(par3 = fit.a$par,
                    par3.conv = fit.a$convergence,
                    Haz0 = approxfun(sort(unique(yi2)), Haz, yleft = min(Haz), yright = max(Haz))))
    }
}

temAR <- function(DF, eqType, solver, par3, zi, wgt = NULL) {
    df0 <- DF[DF$event == 0,]
    rownames(df0) <- NULL
    xi <- as.matrix(df0[,-c(1:6)])    
    di <- df0$terminal
    yi <- df0$time2
    p <- ncol(xi)
    if (is.null(eqType)) {
        yi <- log(yi) + xi %*% par3
        yi2 <- sort(unique(yi))
        Haz <- apply(wgt, 2, function(e) temHaz(par3, rep(0, p), xi, yi, zi / mean(zi), di, e, yi2))
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
    if (is.null(solver)) return(U1(par3))
    else {
        fit.a <- eqSolve(par3, U1, solver)
        yi <- log(yi) + xi %*% fit.a$par
        yi2 <- sort(unique(yi))
        Haz <- c(temHaz(fit.a$par, rep(0, p), xi, yi, zi / mean(zi), di, wi, yi2))
        ind <- !duplicated(exp(yi2))
        list(par3 = fit.a$par,
             par3.conv = fit.a$convergence,
             Haz0 = approxfun(exp(yi2)[ind], Haz[ind], yleft = min(Haz), yright = max(Haz)))
    }
}
