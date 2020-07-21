## globalVariables("DF") ## global variables for reReg

##############################################################################
## Functions for different methods
## stdErr is estimated with resampling if method = sc or am.xc,
##        bootstrap otherwise
##############################################################################

regFit.am.GL <- function(DF, engine, stdErr) {
    DF0 <- DF[DF$event == 0,]
    p <- ncol(DF0) - 6
    alpha <- beta <- gamma <- rep(0, p)
    Y <- log(DF0$time2)
    X <- as.matrix(DF0[,-(1:6)])
    status <- DF0$terminal
    n <- nrow(DF0)
    log.est <- function(b) {
        .C("log_ns_est", as.double(b), as.double(Y), as.double(X), as.double(status),
           as.integer(rep(1, n)), as.integer(n), as.integer(p), as.integer(n),
           as.double(rep(1, n)), as.double(rep(1, n)), 
           result = double(p), PACKAGE = "reReg")$result
    }
    fit.b <- eqSolve(engine@b0, log.est, engine@solver)
    bhat <- fit.b$par
    m <- aggregate(event ~ id, data = DF, sum)[,2]
    index <- c(1, cumsum(m)[-n] + 1)
    ghoshU2 <- function(a) {
        d <- max(X %*% (a - bhat))
        tij <- log(DF$time2) - as.matrix(DF[,-(1:6)]) %*% a
        tij <- tij[DF$event == 1]
        ## yi <- Y - X %*% a - d ## Correct version
        yi <- Y - X %*% bhat - d ## Paper version
        if (sum(tij < rep(yi, m)) == 0) return(1e5)
        else
            .C("glU2", as.integer(n), as.integer(p), as.integer(index - 1), as.integer(m),
               as.double(yi), as.double(tij), as.double(X), as.double(rep(1, n)), result = double(p),
               PACKAGE = "reReg")$result
    }
    fit.a <- eqSolve(engine@a0, ghoshU2, engine@solver)
    fit.b$par <- -fit.b$par
    fit.a$par <- -fit.a$par
    out <- list(alpha = fit.a$par, aconv = fit.a$convergence,
                beta = fit.b$par, bconv = fit.b$convergence, muZ = NA)
    out$recType <- engine@recType
    out$temType <- engine@temType
    return(out)
}

regFit.am.GL.resampling <- function(DF, engine, stdErr) {
    res <- regFit(DF, engine, NULL)
    DF0 <- DF[DF$event == 0,]
    p <- ncol(DF0) - 6
    alpha <- beta <- gamma <- rep(0, p)
    Y <- log(DF0$time2)
    X <- as.matrix(DF0[,-(1:6)])
    status <- DF0$terminal
    n <- nrow(DF0)
    m <- aggregate(event ~ id, data = DF, sum)[,2]
    index <- c(1, cumsum(m)[-n] + 1)
    B <- stdErr@B
    E <- matrix(rexp(n * B), nrow = n)
    Z <- matrix(rnorm(p * B), nrow = p)
    Sn <- function(a, b, w, r = "both") {
        s2 <- .C("log_ns_est", as.double(b), as.double(Y), as.double(X), as.double(status),
                 as.integer(rep(1, n)), as.integer(n), as.integer(p), as.integer(n),
                 as.double(w), as.double(rep(1, n)), 
                 result = double(p), PACKAGE = "reReg")$result
        d <- max(X %*% (a - b))
        tij <- log(DF$time2) - as.matrix(DF[,-(1:6)]) %*% a
        tij <- tij[DF$event == 1]
        ## yi <- Y - X %*% a - d ## Correct version
        yi <- Y - X %*% b - d ## Paper version
        if (sum(tij < rep(yi, m)) == 0) s1 <- NA
        else s1 <- .C("glU2", as.integer(n), as.integer(p), as.integer(index - 1), as.integer(m),
                      as.double(yi), as.double(tij), as.double(X), as.double(w), result = double(p),
                      PACKAGE = "reReg")$result
        if (r == "s1") return(s1)
        if (r == "s2") return(s2)
        return(c(s1, s2))
    }
    V <- var(t(apply(E, 2, function(x) Sn(-res$alpha, -res$beta, x))))
    V1 <- V[1:p, 1:p]
    V2 <- V[1:p + p, 1:p + p]
    lmfit1 <- t(apply(Z, 2, function(x) Sn(-res$alpha + x / sqrt(n), -res$beta, rep(1, n), "s1")))
    lmfit2 <- t(apply(Z, 2, function(x) Sn(-res$alpha, -res$beta + x / sqrt(n), rep(1, n), "s2")))
    if (p == 1) {
        J1 <- coef(lm(sqrt(n) * c(lmfit1) ~ c(Z)))[-1]
        J2 <- coef(lm(sqrt(n) * c(lmfit2) ~ c(Z)))[-1]
    } else {        
        J1 <- coef(lm(sqrt(n) * lmfit1 ~ t(Z)))[-1,]
        J2 <- coef(lm(sqrt(n) * lmfit2 ~ t(Z)))[-1,]
    }
    if (qr(J1)$rank == p) aVar <- solve(J1) %*% V1 %*% t(solve(J1))
    else aVar <- ginv(J1) %*% V1 %*% t(ginv(J1))
    if (qr(J2)$rank == p) bVar <- solve(J2) %*% V2 %*% t(solve(J2))
    else bVar <- ginv(J2) %*% V2 %*% t(ginv(J2))    
    aSE <- sqrt(diag(aVar))
    bSE <- sqrt(diag(bVar))
    out <- c(res, list(alphaSE = aSE, betaSE = bSE, alphaVar = aVar, betaVar = bVar))
    out$recType <- engine@recType
    out$temType <- engine@temType
    return(out)
}

#' @importFrom survival cluster
#' @importFrom stats vcov
regFit.cox.LWYY <- function(DF, engine, stdErr) {
    id <- DF$id
    event <- DF$event
    X <- as.matrix(DF[,-c(1:6)])
    p <- ncol(X)
    T <- DF$time2
    T0 <- unlist(lapply(split(T, id), function(x) c(0, x[-length(x)])))
    fit.coxph <- coxph(Surv(T0, T, event) ~ X + cluster(id))
    out <- list(alpha = coef(fit.coxph), alphaSE = sqrt(diag(vcov(fit.coxph))))
    out$recType <- engine@recType
    out$temType <- engine@temType
    return(out)
}

#' This is also the ARF in Luo
#' @importFrom survival basehaz
#' @noRd
regFit.cox.GL <- function(DF, engine, stdErr) {
    id <- DF$id
    event <- DF$event
    X <- as.matrix(DF[,-c(1:6)])
    p <- ncol(X)
    T <- DF$time2
    mt <- aggregate(event ~ id, data = DF, sum)$event
    Y <- rep(DF$time2[event == 0], mt + 1)
    X0 <- X[event == 0,,drop = FALSE]
    fit.coxph <- coxph(Surv(T[event == 0], DF$terminal[event == 0]) ~ X0)
    cumHaz <- basehaz(fit.coxph)
    ## cumHaz$hazard <- cumHaz$hazard / max(cumHaz$hazard)
    wgt <- sapply(exp(X0 %*% coef(fit.coxph)), function(x)
        approxfun(cumHaz$time, exp(-cumHaz$hazard * x), yleft = 1, yright = min(exp(-cumHaz$hazard * x)),
                  method = "constant")(T))
    wgt <- 1 / wgt ## ifelse(wgt == 0, 1 / sort(c(wgt))[2], 1 / wgt)
    wgt <- ifelse(wgt > 1e5, 1e5, wgt)
    out <- dfsane(par = engine@a0[1:ncol(X)], fn = coxGLeq, wgt = wgt, 
                  X = as.matrix(X[!event, ]),
                  Y = Y[!event], T = ifelse(T == Y, 1e5, T), cl = mt + 1,
                  alertConvergence = FALSE, quiet = TRUE, control = list(trace = FALSE))
    out <- list(alpha = out$par,
                beta = coef(fit.coxph),
                betaSE = sqrt(diag(vcov(fit.coxph))))
    out$recType <- engine@recType
    out$temType <- engine@temType
    return(out)}

## #' @importFrom rlang is_empty
regFit.general <- function(DF, engine, stdErr) {
    if (is.na(match(engine@solver, c("dfsane", "BBsolve", "optim", "BBoptim")))) {
        print("Warning: Unidentified solver; BB::dfsane is used.")
        engine@solver <- "dfsane"
    }
    out <- s1(engine@recType, DF, engine@eqType, engine@solver, engine@a0)
    if (engine@temType != ".")
        out <- c(out, s2(engine@temType, DF, engine@eqType, engine@solver, engine@b0, out$zi))
    out$recType <- engine@recType
    out$temType <- engine@temType
    return(out)
}

s1 <- function(type, DF, eqType, solver, a0, wgt = NULL) {
    if (type == "sc") return(reSC(DF, eqType, solver, a0, wgt))
    if (type == "cox") return(reCox(DF, eqType, solver, a0, wgt))
    if (type == "am") return(reAM(DF, eqType, solver, a0, wgt))
    if (type == "ar") return(reAR(DF, eqType, solver, a0, wgt))
}
s2 <- function(type, DF, eqType, solver, b0, zi, wgt = NULL) {
    if (type == "sc") return(temSC(DF, eqType, solver, b0, zi, wgt))
    if (type == "cox") return(temCox(DF, eqType, solver, b0, zi, wgt))
    if (type == "am") return(temAM(DF, eqType, solver, b0, zi, wgt))
    if (type == "ar") return(temAR(DF, eqType, solver, b0, zi, wgt))
}

regFit.general.resampling <- function(DF, engine, stdErr) {
    if (is.na(match(engine@solver, c("dfsane", "BBsolve", "optim", "BBoptim")))) {
        print("Warning: Unidentified solver; BB::dfsane is used.")
        engine@solver <- "dfsane"
    }
    res <- regFit(DF, engine, NULL)
    n <- length(unique(DF$id))
    B <- stdErr@B
    E1 <- matrix(rexp(n * B), n)
    E2 <- matrix(rexp(n * B), n)
    p <- ncol(DF) - 6
    a0 <- res$alpha
    b0 <- res$beta
    if (engine@recType == "sc")
        a0 <- c(res$alpha[1:p], res$log.muZ, res$alpha[1:p + p] - res$alpha[1:p])
    if (engine@recType == "cox")
        a0 <- c(res$log.muZ, res$alpha)    
    tmpV <- sapply(1:B, function(ee) {
        tmp <- s1(engine@recType, DF, engine@eqType, NULL, a0, E1[,ee])
        c(tmp$value, s2(engine@temType, DF, engine@eqType, NULL, b0, E2[,ee] * tmp$zi, E2[,ee]))
    })
    V <- var(t(tmpV))
    Z <- matrix(rnorm(ncol(V) * B), B)
    na <- length(a0)
    nb <- length(b0)
    L <- apply(Z, 1, function(zz) {
        tmp <- s1(engine@recType, DF, engine@eqType, NULL, a0 + zz[1:na] / sqrt(n))
        c(tmp$value, s2(engine@temType, DF, engine@eqType, NULL, b0 + tail(zz, nb) / sqrt(n), tmp$zi))
    })
    L <- t(L)
    J <- solve(t(Z) %*% Z) %*% t(Z) %*% (sqrt(n) * L)
    aVar <- varMat <- solve(J[1:na, 1:na]) %*% V[1:na, 1:na] %*% t(solve(J[1:na, 1:na]))
    if (engine@recType == "cox") aVar <- varMat[-1, -1]
    if (engine@recType == "sc") {
        aVar <- varMat[-(p + 1), -(p + 1)]
    }
    res <- c(res, list(alphaSE = sqrt(diag(as.matrix(aVar))), alphaVar = aVar, varMat = varMat))
    if (nb > 0) {
        ind2 <- tail(1:nrow(J), nb)
        J2 <- solve(t(Z[,ind2]) %*% Z[,ind2]) %*% t(Z[,ind2]) %*% (sqrt(n) * L[,ind2])
        bVar <- solve(J2) %*% V[ind2, ind2] %*% t(solve(J2))
        res$betaSE <- sqrt(diag(as.matrix(bVar)))
        res$betaVar <- bVar
    } 
    return(res)
}

##############################################################################
# Variance estimation 
##############################################################################
regFit.Engine.Bootstrap <- function(DF, engine, stdErr) {
    res <- regFit(DF, engine, NULL)
    ## engine@a0 <- res$alpha
    ## engine@b0 <- res$beta
    id <- DF$id
    event <- DF$event
    status <- DF$terminal
    X <- as.matrix(DF[,-c(1:6)])    
    n <- length(unique(id))
    p <- ncol(X)
    T <- DF$time2
    mt <- aggregate(event ~ id, data = DF, sum)$event
    clsz <- mt + 1
    Y <- rep(DF$time2[event == 0], mt + 1)
    cluster <- unlist(sapply(mt + 1, function(x) 1:x))
    B <- stdErr@B
    uID <- unique(DF$id) # unique(DF$ID)
    if (stdErr@parallel) {
        cl <- makeCluster(stdErr@parCl)
        clusterExport(cl = cl,
                      varlist = c("DF", "engine"),
                      envir = environment())
        out <- parSapply(cl, 1:B, function(x) {
            sampled.id <- sample(unique(id), n, TRUE)
            ind <- unlist(sapply(sampled.id, function(x) which(id == x)))
            DF2 <- DF[ind,]
            DF2$id <- rep(1:n, clsz[sampled.id])
            tmp <- regFit(DF2, engine, NULL)
            return(c(tmp$alpha, tmp$beta))
        })
        stopCluster(cl)
        betaMatrix <- t(out)
        convergence <- apply(betaMatrix, 1, function(x)
            1 * (x %*% x > 1e3 * c(res$alpha, res$beta) %*% c(res$alpha, res$beta)))
    } else {
        betaMatrix <- matrix(0, B, p * 2)
        convergence <- rep(0, B)
            for (i in 1:B) {
            sampled.id <- sample(unique(id), n, TRUE)
            ind <- unlist(sapply(sampled.id, function(x) which(id == x)))
            DF2 <- DF[ind,]
            DF2$id <- rep(1:n, clsz[sampled.id])
            tmp <- regFit(DF2, engine, NULL)
            betaMatrix[i,] <- c(tmp$alpha, tmp$beta)
            convergence[i] <- 1 * (betaMatrix[i,] %*% betaMatrix[i,] >
                                   1e3 * c(res$alpha, res$beta) %*% c(res$alpha, res$beta))
        }
    }
    converged <- which(convergence == 0)
    if (sum(convergence != 0) > 0) {
        print("Warning: Some bootstrap samples failed to converge")
        tmp <- apply(betaMatrix, 1, function(x) x %*% x)
        converged <- (1:B)[- which(tmp %in% boxplot(tmp, plot = FALSE)$out)]        
    }
    if (all(convergence != 0) || sum(convergence == 0) == 1) {
        print("Warning: some bootstrap samples failed to converge")
        converged <- 1:B
    }
    betaVar <- var(betaMatrix[converged, ], na.rm = TRUE)
    betaSE <- sqrt(diag(as.matrix(betaVar)))
    c(res, list(alphaSE = betaSE[1:p], betaSE = betaSE[1:p + p],
                alphaVar = as.matrix(betaVar[1:p, 1:p]),
                betaVar = as.matrix(betaVar[1:p + p, 1:p + p]),
                SEmat = betaMatrix, B = length(converged)))
}

##############################################################################
# Nonparametric (~1)
##############################################################################

## ~1
npFit <- function(DF, B = 0) {
    df0 <- DF[DF$event == 0,]
    df1 <- DF[DF$event == 1,]
    rownames(df0) <- rownames(df1) <- NULL
    m <- aggregate(event ~ id, data = DF, sum)[, 2]
    yi <- df0$time2
    ti <- df1$time2
    zi <- wi <- rep(1, length(ti))
    di <- df0$terminal
    xi <- as.matrix(df0[,-c(1:6)])
    p <- ncol(xi)
    yi2 <- sort(unique(yi))
    if (B > 0) {
        n <- length(unique(DF$id))
        E1 <- matrix(rexp(n * B), n)
        E2 <- matrix(rexp(n * B), n)
        rate <- apply(E1, 2, function(e) reRate(ti, rep(yi, m), rep(e, m), yi))
        rate <- apply(rate, 1, quantile, c(.025, .975))
        Lam <- exp(-rate)
        Haz <- apply(E2, 2, function(e) temHaz(rep(0, p), rep(0, p), xi, yi, zi, di, e, yi2))
        Haz <- apply(Haz, 1, quantile, c(.025, .975))
        return(list(Lam0.lower = approxfun(yi[!duplicated(yi)], Lam[2, !duplicated(yi)],
                                           yleft = min(Lam[2,]), yright = max(Lam[2,])),
                    Lam0.upper = approxfun(yi[!duplicated(yi)], Lam[1, !duplicated(yi)],
                                           yleft = min(Lam[1,]), yright = max(Lam[1,])),
                    Haz0.lower = approxfun(yi2, Haz[1,],
                                           yleft = min(Haz[1,]), yright = max(Haz[1,])),
                    Haz0.upper = approxfun(yi2, Haz[2,],
                                           yleft = min(Haz[2,]), yright = max(Haz[2,]))))
    } else {    
        rate <- c(reRate(ti, rep(yi, m), wi, yi))
        Lam <- exp(-rate)
        Lam0 <- approxfun(yi[!duplicated(yi)], Lam[!duplicated(yi)],
                          yleft = min(Lam), yright = max(Lam))
        Haz <- c(temHaz(rep(0, p), rep(0, p), xi, yi, zi, di, wi, yi2))
        Haz0 <- approxfun(yi2, Haz, yleft = min(Haz), yright = max(Haz))
        return(list(Lam0 = Lam0, Haz0 = Haz0))
    }    
}

npFitSE <- function(DF, recType, temType, a0, b0, zi, B) {
    n <- length(unique(DF$id))
    E1 <- matrix(rexp(n * B), n)
    E2 <- matrix(rexp(n * B), n)
    c(s1(recType, DF, NULL, NULL, a0, E1),
      s2(temType, DF, NULL, NULL, b0, zi, E2))
}

##############################################################################
# Class Definition
##############################################################################

setClass("Engine",
         representation(tol = "numeric", a0 = "numeric", b0 = "numeric",
                        baseSE = "logical", 
                        solver = "character", eqType = "character", 
                        recType = "character", temType = "character"),
         prototype(eqType = "logrank", tol = 1e-7, a0 = 0, b0 = 0, baseSE = FALSE, solver = "dfsane"),
         contains = "VIRTUAL")
setClass("general", contains = "Engine")
setClass("cox.LWYY", contains = "Engine")
setClass("cox.HW", contains = "Engine")
setClass("am.XCHWY", contains = "Engine")
setClass("am.GL", contains = "Engine")
setClass("sc.XCYH", representation(muZ = "numeric"),
         prototype(muZ = 0), contains = "Engine")
setClass("cox.GL",
         representation(wgt = "matrix"), prototype(wgt = matrix(0)), contains = "Engine")

setClass("stdErr",
         representation(B = "numeric", parallel = "logical", parCl = "numeric"),
         prototype(B = 100, parallel = FALSE, parCl = parallel::detectCores() / 2),
         contains = "VIRTUAL")

setClass("bootstrap", contains = "stdErr")
setClass("resampling", contains = "stdErr")


##############################################################################
# Method Dispatch
##############################################################################
setGeneric("regFit", function(DF, engine, stdErr) {standardGeneric("regFit")})

setMethod("regFit", signature(engine = "general", stdErr = "NULL"), regFit.general)
setMethod("regFit", signature(engine = "general", stdErr = "resampling"), regFit.general.resampling)
setMethod("regFit", signature(engine = "cox.LWYY", stdErr = "NULL"), regFit.cox.LWYY)
setMethod("regFit", signature(engine = "cox.LWYY", stdErr = "bootstrap"), regFit.cox.LWYY)
setMethod("regFit", signature(engine = "cox.LWYY", stdErr = "resampling"), regFit.cox.LWYY)
setMethod("regFit", signature(engine = "cox.GL", stdErr = "NULL"), regFit.cox.GL)
setMethod("regFit", signature(engine = "cox.GL", stdErr = "resampling"), regFit.cox.GL)
setMethod("regFit", signature(engine = "am.GL", stdErr = "NULL"), regFit.am.GL)
setMethod("regFit", signature(engine = "Engine", stdErr = "bootstrap"),
          regFit.Engine.Bootstrap)
setMethod("regFit", signature(engine = "am.GL", stdErr = "resampling"),
          regFit.am.GL.resampling)


#' Fits Semiparametric Regression Models for Recurrent Event Data
#'
#' Fits a general (joint) semiparametric regression model for the recurrent event data,
#' where the rate function of the underlying recurrent event process and
#' the hazard function of the terminal event can be specified as a Cox-type model,
#' an accelerated mean model, an accelerated rate model, or a generalized scale-change model.
#' See details for model specifications.
#'
#'
#' Suppose the recurrent event process and the failure events are
#' observed in the time interval \eqn{t\in[0,\tau]},
#' for some constant \eqn{\tau}.
#' We formulate the recurrent event rate function, \eqn{\lambda(t)},
#' and the terminal event hazard function, \eqn{h(t)}, 
#' in the form of
#' \deqn{\lambda(t) = Z \lambda_0(te^{X^\top\alpha}) e^{X^\top\beta}, h(t) = Z h_0(te^{X^\top\eta})e^{X^\top\theta},}
#' where \eqn{\lambda_0(t)} is the baseline rate function,
#' \eqn{h_0(t)} is the baseline hazard function,
#' \eqn{X} is a \eqn{n} by \eqn{p} covariate matrix and \eqn{\alpha},
#' \eqn{Z} is an unobserved shared frailty variable, and
#' \eqn{(\alpha, \eta)} and \eqn{(\beta, \theta)} correspond to the shape and size parameters of the
#' rate function and the hazard function, respectively.
#' The model includes several popular semiparametric models as special cases,
#' which can be specified via the \code{method} argument with the rate function
#' and hazard function separated by "\code{|}".
#' For examples,
#' Wang, Qin and Chiang (2001) (\eqn{\alpha = \eta = \theta = 0})
#' can be called with \code{method = "cox|."};
#' Huang and Wang (2004) (\eqn{\alpha = \eta = 0})
#' can be called with \code{method = "cox|cox"};
#' Xu et al. (2017) (\eqn{\alpha = \beta} and \eqn{\eta = \theta})
#' can be called with \code{method = "am|am"};
#' Xu et al. (2019) (\eqn{\eta = \theta = 0}) can be called with \code{method = "sc|."}.
#' Users can mix the models depending on the application. For example,
#' \code{method = "cox|ar"} postulate a Cox proportional model for the
#' recurrent event rate function and an accelerated rate model for
#' the terminal event hazard function (\eqn{\alpha = \theta = 0}).
#' If only one method is specified without an "\code{|}",
#' it is used for both the rate function and the hazard function.
#' For example, specifying \code{method = "cox"} is equivalent to \code{method = "cox|cox"}.
#' Some methods that assumes \code{Z = 1} and requires independent
#' censoring are also implemented in \code{reReg};
#' these includes \code{method = "cox.LWYY"} for Lin et al. (2000),
#' \code{method = "cox.GL"} for Ghosh and Lin (2002),
#' and \code{method = "am.GL"} for Ghosh and Lin (2003).
#'
#' The available methods for variance estimation are:
#' \describe{
#'   \item{NULL}{variance estimation will not be performed. This is equivalent to setting \code{B = 0}.}
#'   \item{resampling}{performs the efficient resampling-based variance estimation.}
#'   \item{bootstrap}{performs nonparametric bootstrap.}
#' }
#'
#' The \code{control} list consists of the following parameters:
#' \describe{
#'   \item{tol}{absolute error tolerance.}
#'   \item{a0, b0}{initial guesses used for root search.}
#'   \item{solver}{the equation solver used for root search. The available options are \code{BB::BBsolve}, \code{BB::dfsane}, \code{BB:BBoptim}, and \code{optim}.}
#'   \item{baseSE}{an logical value indicating whether the 95\% confidence bounds for the baseline functions will be computed.}
#'   \item{parallel}{an logical value indicating whether parallel computation will be applied when \code{se = "bootstrap"} is called.}
#'   \item{parCl}{an integer value specifying the number of CPU cores to be used when \code{parallel = TRUE}. The default value is half the CPU cores on the current host.}
#' }
#' 
#' @param formula a formula object, with the response on the left of a "~" operator, and the predictors on the right.
#' The response must be a recurrent event survival object as returned by function \code{Recur}.
#' @param data  an optional data frame in which to interpret the variables occurring in the \code{"formula"}.
#' @param B a numeric value specifies the number of resampling for variance estimation.
#' When \code{B = 0}, variance estimation will not be performed.
#' @param method a character string specifying the underlying model. See \bold{Details}.
#' @param se a character string specifying the method for standard error estimation. See \bold{Details}.
#' @param control a list of control parameters.
#'
#' @export
#' @references Lin, D., Wei, L., Yang, I. and Ying, Z. (2000). Semiparametric Regression for the Mean and Rate Functions of Recurrent Events.
#' \emph{Journal of the Royal Statistical Society: Series B (Methodological)}, \bold{62}: 711--730.
#' @references Wang, M.-C., Qin, J., and Chiang, C.-T. (2001). Analyzing Recurrent Event Data with Informative Censoring.
#' \emph{Journal of the American Statistical Association}, \bold{96}(455): 1057--1065.
#' @references Ghosh, D. and Lin, D.Y. (2002). Marginal Regression Models for Recurrent and Terminal Events. \emph{Statistica Sinica}: 663--688.
#' @references Ghosh, D. and Lin, D.Y. (2003). Semiparametric Analysis of Recurrent Events Data in the Presence of Dependent Censoring.
#' \emph{Biometrics}, \bold{59}: 877--885.
#' @references Huang, C.-Y. and Wang, M.-C. (2004). Joint Modeling and Estimation for Recurrent Event Processes and Failure Time Data.
#' \emph{Journal of the American Statistical Association}, \bold{99}(468): 1153--1165.
#' @references Xu, G., Chiou, S.H., Huang, C.-Y., Wang, M.-C. and Yan, J. (2017). Joint Scale-change Models for Recurrent Events and Failure Time.
#' \emph{Journal of the American Statistical Association}, \bold{112}(518): 796--805.
#' @references Xu, G., Chiou, S.H.,Yan, J., Marr, K., and Huang, C.-Y. (2019). Generalized Scale-Change Models for Recurrent Event
#' Processes under Informative Censoring. \emph{Statistica Sinica}: pre-print.
#'
#' @importFrom stats approxfun optim
#' 
#' @seealso \code{\link{Recur}}, \code{\link{simSC}}
#'
#' @example inst/examples/ex_reReg.R
reReg <- function(formula, data, B = 200, 
                  method = "cox", se = c("resampling", "bootstrap", "NULL"),
                  control = list()) {
    se <- match.arg(se)
    Call <- match.call()
    if (missing(data)) obj <- eval(formula[[2]], parent.frame()) 
    if (!missing(data)) obj <- eval(formula[[2]], data) 
    if (!is.Recur(obj)) stop("Response must be a `Recur` object")
    formula[[2]] <- NULL
    if (formula == ~ 1) {
        DF <- as.data.frame(cbind(obj@.Data, zero = 0))
    } else {
        ## remove intercept
        if (!missing(data)) DF <- as.data.frame(cbind(obj@.Data, model.matrix(formula, data)))
        if (missing(data)) DF <- as.data.frame(cbind(obj@.Data, model.matrix(formula, parent.frame())))
        DF <- DF[,-which(colnames(DF) == "(Intercept)")]
    }
    DF <- DF[order(DF$id, DF$time2), ]
    allMethod <- apply(expand.grid(c("cox", "am", "sc", "ar"), c("cox", "am", "sc", "ar", ".")), 1,
                       paste, collapse = "|")
    allMethod <- c(allMethod, "cox.LWYY", "cox.GL", "cox.HW", "am.GL", "am.XCHWY", "sc.XCYH")
    method <- match.arg(method, c("cox", "am", "sc", "ar", allMethod))
    recType <- temType <- NULL
    if (grepl("|", method, fixed = TRUE)) {
        recType <- substring(method, 1, regexpr("[|]", method) - 1)
        temType <- substring(method, regexpr("[|]", method) + 1)
        method <- "general"
    }
    if (method %in% c("cox", "am", "sc", "ar")) {
        recType <- temType <- method
        method <- "general"
    }
    ## Special cases:
    if (method == "cox.HW") {
        recType <- temType <- "cox"
        method <- "general"
    }
    if (method == "am.XCHWY") {
        recType <- temType <- "am"
        method <- "general"
    }
    if (method == "sc.XCYH") {
        recType <- "sc"
        temType <- "."
        method <- "general"        
    }
    if (method == "cox.LWYY") {
        recType <- "cox.LWYY"
        temType <- "."
    }
    if (method == "cox.GL") {
        recType <- "cox.GL"
        temType <- "cox.GL"
    }
    if (method == "am.GL") {
        recType <- "am.GL"
        temType <- "am.GL"
    }    
    engine.control <- control[names(control) %in% names(attr(getClass(method), "slots"))]
    engine <- do.call("new", c(list(Class = method), engine.control))
    engine@recType <- recType
    engine@temType <- temType
    if (se == "NULL" || B == 0)
        stdErr <- NULL
    else {
        stdErr.control <- control[names(control) %in% names(attr(getClass(se), "slots"))]
        stdErr <- do.call("new", c(list(Class = se), stdErr.control))
        stdErr@B <- B
    }
    p <- ncol(DF) - ncol(obj@.Data)
    if (length(engine@a0) == 1 & any(grepl("sc", c(method, engine@recType), fixed = FALSE)))
        engine@a0 <- rep(engine@a0, 2 * p + 1)
    if (length(engine@a0) == 1 & any(grepl("cox", c(method, engine@recType), fixed = FALSE)))
        engine@a0 <- rep(engine@a0, p + 1)
    if (length(engine@a0) == 1 & any(grepl("ar|am", c(method, engine@recType), fixed = FALSE)))    
        engine@a0 <- rep(engine@a0, p)
    if (length(engine@b0) == 1) {
        if (any(grepl("sc", c(method, engine@temType), fixed = FALSE)))
            engine@b0 <- rep(engine@b0, 2 * p)
        else engine@b0 <- rep(engine@b0, p)
    }
    if (formula == ~1) {
        if (engine@baseSE) fit <- npFit(DF, B)
        else fit <- npFit(DF)
        fit$method <- "nonparametric"
    } else {
        fit <- regFit(DF = DF, engine = engine, stdErr = stdErr)
        if (method == "general" & engine@baseSE) {
            if (fit$recType == "sc") {
                a0 <- c(fit$alpha[1:p], fit$log.muZ, fit$alpha[1:p + p] - fit$alpha[1:p])
                fit <- c(fit, npFitSE(DF, fit$recType, fit$temType, a0, fit$beta, fit$zi, B))
            }
            if (fit$recType != "sc") 
                fit <- c(fit, npFitSE(DF, fit$recType, fit$temType, fit$alpha, fit$beta, fit$zi, B))
        }
        fit$method <- method
    }    
    class(fit) <- "reReg"
    fit$reTb <- obj@.Data
    fit$DF <- DF
    fit$call <- Call
    fit$varNames <- names(DF)[-(1:6)]
    fit$se <- se
    fit
}

#' Equation wrapper
#'
#' @noRd
#' @importFrom BB spg
#' @importFrom rootSolve uniroot.all
#' @keywords internal
eqSolve <- function(par, fn, solver, ...) {
    if (length(fn(par, ...)) == 1) {
        tmp <- uniroot.all(Vectorize(fn), interval = c(par - 10, par + 10))
        out <- NULL
        out$par <- tmp[which.min(abs(tmp - par))]
        out$convergence <- 0 ## 1 * (abs(fn(out$par, ...)) < 1e-5)
        return(out)
    }
    if (solver == "dfsane") {
        out <- dfsane(par = par, fn = function(z) fn(z, ...), 
                      alertConvergence = FALSE, quiet = TRUE, control = list(trace = FALSE))
        if (max(abs(out$par)) > 10) solver <- "BBsolve"
    }
    if (solver == "BBsolve")
        out <- BBsolve(par = par, fn = fn, ..., quiet = TRUE)
    if (solver == "BBoptim")
        out <- BBoptim(par = par, fn = function(z) sum(fn(z, ...)^2),
                       quiet = TRUE, control = list(trace = FALSE))
    if (solver == "optim")
        out <- optim(par = par, fn = function(z) sum(fn(z, ...)^2),
                     control = list(trace = FALSE))
    return(out)
}

##############################################################################
## Background functions...
## Probably need to clean these up someday
##############################################################################

#' R function for equation 8 of Ghosh & Lin (2002);
#' Marginal regression models for recurrent and terminal events.
#' 
#' @keywords internal
#' @noRd
coxGLeq <- function(beta, X, Y, T, cl, wgt) {
    p <- ncol(X)
    res <- vector("double", p)
    xb <- exp(X %*% beta)
    .C("coxGL", as.double(T), as.double(Y), as.double(X), as.double(xb), as.double(wgt),
       as.integer(length(T)), as.integer(cl), as.integer(c(0, cumsum(cl)[-length(cl)])),
       as.integer(nrow(X)), as.integer(p),        
       out = double(p), PACKAGE = "reReg")$out       
}

## varOut <- function(dat, na.rm = TRUE) {
##     dat[which(dat %in% boxplot(dat, plot = FALSE)$out)] <- NA
##     dat <- dat[complete.cases(dat),]
##     var(dat, na.rm = na.rm)
## }
