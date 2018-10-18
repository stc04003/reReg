globalVariables(c("event", "vcov", "wgt")) ## global variables for reReg

##############################################################################
## Functions for different methods
## stdErr is estimated with resampling if method = sc or am.xc,
##        bootstrap otherwise
##############################################################################

regFit.am.XCHWY <- function(DF, engine, stdErr) {
    id <- DF$id
    event <- DF$event
    status <- DF$status
    X <- as.matrix(DF[,-c(1:4)])
    n <- length(unique(id))
    p <- ncol(X)
    T <- DF$Time
    mt <- aggregate(event ~ id, data = DF, sum)$event
    Y <- rep(DF$Time[event == 0], mt + 1)
    cluster <- unlist(sapply(mt + 1, function(x) 1:x))
    outA <- eqSolve(engine@a0, alphaEq, engine@solver,
                    X = X, Y = Y, T = T, cluster = cluster, mt = mt, weights = NULL)
    alpha <- outA$par
    Ystar <- log(Y) + X %*% alpha
    Tstar <- log(T) + X %*% alpha
    lambda <- npMLE(Ystar[event == 0], Tstar, Ystar)
    zHat <- as.numeric(mt * npMLE(log(max(Y)), Tstar, Ystar) / lambda)
    if (any(zHat > 1e5)) zHat <- (mt * npMLE(log(max(Y)), Tstar, Ystar) + .01) / (lambda + .01)
    ## zHat <- ifelse(zHat > 1e5, (mt * npMLE(log(max(Y)), Tstar, Ystar) + .01) / (lambda + .01), zHat)
    zHat <- ifelse(is.na(zHat), 0, zHat)
    outB <- eqSolve(engine@b0, betaEq, engine@solver,
                    X = X, Y = Y, T = T, cluster = cluster, delta = status[event == 0],
                    mt = mt, alpha = outA$par, zHat = zHat, weights = NULL)
    list(alpha = outA$par, aconv = outA$convergence, beta = outB$par, bconv = outB$convergence,
         muZ = mean(zHat, na.rm = TRUE), zHat = zHat)
}

regFit.am.GL <- function(DF, engine, stdErr) {
    DF0 <- subset(DF, event == 0)
    p <- ncol(DF0) - 4
    alpha <- beta <- gamma <- rep(0, p)
    Y <- log(DF0$Time)
    X <- as.matrix(DF0[,-(1:4)])
    status <- DF0$status
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
        tij <- log(DF$Time) - as.matrix(DF[,-(1:4)]) %*% a
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
    list(alpha = fit.a$par, aconv = fit.a$convergence,
         beta = fit.b$par, bconv = fit.b$convergence, muZ = NA)
}

regFit.am.GL.resampling <- function(DF, engine, stdErr) {
    res <- regFit(DF, engine, NULL)
    DF0 <- subset(DF, event == 0)
    p <- ncol(DF0) - 4
    alpha <- beta <- gamma <- rep(0, p)
    Y <- log(DF0$Time)
    X <- as.matrix(DF0[,-(1:4)])
    status <- DF0$status
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
        tij <- log(DF$Time) - as.matrix(DF[,-(1:4)]) %*% a
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
    J1 <- coef(lm(sqrt(n) * lmfit1 ~ t(Z)))[-1,]
    J2 <- coef(lm(sqrt(n) * lmfit2 ~ t(Z)))[-1,]
    if (qr(J1)$rank == p) aVar <- solve(J1) %*% V1 %*% t(solve(J1))
    else aVar <- ginv(J1) %*% V1 %*% t(ginv(J1))
    if (qr(J2)$rank == p) bVar <- solve(J2) %*% V2 %*% t(solve(J2))
    else bVar <- ginv(J2) %*% V2 %*% t(ginv(J2))    
    aSE <- sqrt(diag(aVar))
    bSE <- sqrt(diag(bVar))
    c(res, list(alphaSE = aSE, betaSE = bSE, alphaVar = aVar, betaVar = bVar))
}

regFit.cox.HW <- function(DF, engine, stdErr) {
    id <- DF$id
    event <- DF$event
    status <- DF$status
    X <- as.matrix(DF[,-c(1:4)])    
    n <- length(unique(id))
    p <- ncol(X)
    T <- DF$Time
    mt <- aggregate(event ~ id, data = DF, sum)$event
    Y <- rep(DF$Time[event == 0], mt + 1)
    cluster <- unlist(sapply(mt + 1, function(x) 1:x))
    X <- cbind(1, X[event == 0,])
    engine@a0 <- c(0, engine@a0)
    outA <- eqSolve(par = engine@a0, fn = HWeq, solver = engine@solver,
                    X = X, Y = Y, T = T, cluster = cluster, mt = mt, weights = NULL)
    muZ <- outA$par[1]
    alpha <- outA$par <- outA$par[-1]
    lambda <- npMLE(Y[event == 0], T, Y)
    zHat <- as.numeric(mt / (lambda * exp(as.matrix(X[, -1]) %*% alpha)))
    zHat <- ifelse(zHat > 1e5, (mt + .01) / (lambda * exp(as.matrix(X[, -1]) %*% alpha) + .01), zHat)
    zHat <- ifelse(is.na(zHat), 0, zHat)
    outB <- eqSolve(par = engine@b0, fn = HWeq2, solver = engine@solver,
                    X = as.matrix(X[,-1]), Y = Y[event == 0], delta = status[event == 0], zHat = zHat/muZ, weights = NULL)
    list(alpha = outA$par, aconv = outA$convergence,
         beta = outB$par, bconv = outB$convergence,
         muZ = muZ, zHat = zHat / muZ)
}

#' @importFrom survival cluster
#' @importFrom stats vcov
regFit.cox.LWYY <- function(DF, engine, stdErr) {
    id <- DF$id
    event <- DF$event
    X <- as.matrix(DF[,-c(1:4)])
    p <- ncol(X)
    T <- DF$Time
    T0 <- unlist(lapply(split(T, id), function(x) c(0, x[-length(x)])))
    fit.coxph <- coxph(Surv(T0, T, event) ~ X + cluster(id))
    list(alpha = coef(fit.coxph), alphaSE = sqrt(diag(vcov(fit.coxph))),
         beta = rep(NA, p), betaSE = rep(NA, p), muZ = NA)
    ## my implementation of coxph estimating equation
    ##
    ## status <- DF$status
    ## n <- length(unique(id))
    ## p <- ncol(X)
    ## mt <- aggregate(event ~ id, data = DF, sum)$event
    ## Y <- rep(DF$Time[event == 0], mt + 1)
    ## cluster <- unlist(sapply(mt + 1, function(x) 1:x))   
    ## out <- dfsane(par = engine@a0, fn = LWYYeq, X = as.matrix(X[event == 0, ]),
    ##               Y = Y[event == 0], T = ifelse(T == Y, 1e5, T), cl = mt + 1,
    ##               alertConvergence = FALSE, quiet = TRUE,
    ##               control = list(NM = FALSE, M = 100, noimp = 50, trace = FALSE))
    ## list(alpha = out$par, beta = rep(0, p), muZ = NA)
}

#' This is also the ARF in Luo
#' @importFrom survival basehaz
#' @noRd
regFit.cox.GL <- function(DF, engine, stdErr) {
    id <- DF$id
    event <- DF$event
    X <- as.matrix(DF[,-c(1:4)])
    p <- ncol(X)
    T <- DF$Time
    mt <- aggregate(event ~ id, data = DF, sum)$event
    Y <- rep(DF$Time[event == 0], mt + 1)
    ## T0 <- unlist(lapply(split(T, id), function(x) c(0, x[-length(x)])))
    ## fit.coxph <- coxph(Surv(T0, T, event) ~ X + cluster(id))
    fit.coxph <- coxph(Surv(Time, status) ~ ., data = subset(DF, !event, select = -c(id, event)))
    cumHaz <- basehaz(fit.coxph)
    ## cumHaz$hazard <- cumHaz$hazard / max(cumHaz$hazard)
    X0 <- subset(X, !event)
    wgt <- sapply(exp(X0 %*% coef(fit.coxph)), function(x)
        with(cumHaz, approxfun(time, exp(-hazard * x), yleft = 1, yright = min(exp(-hazard * x)),
                               method = "constant"))(T))
    wgt <- 1 / wgt ## ifelse(wgt == 0, 1 / sort(c(wgt))[2], 1 / wgt)
    wgt <- ifelse(wgt > 1e5, 1e5, wgt)
    out <- dfsane(par = engine@a0, fn = coxGLeq, wgt = wgt, 
                  X = as.matrix(X[!event, ]),
                  Y = Y[!event], T = ifelse(T == Y, 1e5, T), cl = mt + 1,
                  alertConvergence = FALSE, quiet = TRUE, control = list(trace = FALSE))
    list(alpha = out$par, beta = coef(fit.coxph), betaSE = sqrt(diag(vcov(fit.coxph))), muZ = NA,
         wgt = wgt, haz0 = with(cumHaz, approxfun(time, hazard, yleft = 0, yright = max(hazard), method = "constant")))
}

regFit.sc.XCYH <- function(DF, engine, stdErr) {
    if (is.na(match(engine@solver, c("dfsane", "BBsolve", "optim", "BBoptim")))) {
        print("Warning: Unidentified solver; BB::dfsane is used.")
        engine@solver <- "dfsane"
    }
    DF0 <- subset(DF, event == 0)
    p <- ncol(DF) - 4
    X <- as.matrix(DF0[,-(1:4)])
    n <- nrow(DF0)
    m <- aggregate(event ~ id, data = DF, sum)[,2]
    index <- c(1, cumsum(m)[-n] + 1)
    U1 <- function(a) {
        yi <- log(DF0$Time) + X %*% a
        tij <- log(DF$Time) + as.matrix(DF[,-(1:4)]) %*% a
        tij <- tij[DF$event == 1]
        if (engine@eqType %in% c("Logrank", "logrank")) 
            return(.C("sc1Log", as.integer(n), as.integer(p), as.integer(index - 1),
                      as.integer(m), as.double(yi), as.double(tij), as.double(X),
                      as.double(rep(1, length(index))), 
                      result = double(p), PACKAGE = "reReg")$result)
        if (engine@eqType %in% c("Gehan", "gehan"))
            return(.C("sc1Gehan", as.integer(n), as.integer(p), as.integer(index - 1),
                      as.integer(m), as.double(yi), as.double(tij), as.double(X),
                      as.double(rep(1, length(index))), 
                      result = double(p), PACKAGE = "reReg")$result / n)
    }
    fit.a <- eqSolve(engine@a0, U1, engine@solver)
    ahat <- fit.a$par
    yi <- log(DF0$Time) + X %*% ahat
    tij <- log(DF$Time) + as.matrix(DF[,-(1:4)]) %*% ahat
    tij <- tij[DF$event == 1]
    t0.rate <- unique(sort(c(tij, yi)))
    rate <- .C("scRate", as.integer(n), as.integer(index - 1), as.integer(m),
               as.integer(length(t0.rate)), as.double(rep(1, n)), as.double(yi),
               as.double(tij), as.double(t0.rate), result = double(length(t0.rate)),
               PACKAGE = "reReg")$result
    Lam <- exp(-approx(t0.rate, rate, yi)$y)
    U2 <- function(b) {
        tmp <- m / Lam
        tmp <- ifelse(tmp > 1e5, (m + .01) / (Lam + .01), tmp)
        xb <- exp(cbind(1, X) %*% b)
        xb <- ifelse(xb == Inf, 1e5, xb)
        .C("sc2", as.integer(n), as.integer(p + 1), as.double(cbind(1, X)),
           as.double(xb), as.double(tmp), as.double(rep(1, n)),
           result = double(p + 1), PACKAGE = "reReg")$result / n
    }
    fit.b <- eqSolve(engine@b0, U2, engine@solver)
    list(alpha = fit.a$par, beta = fit.b$par[-1] + fit.a$par,
         aconv = fit.a$convergence, bconv = fit.b$convergence,
         log.muZ = fit.b$par[1])
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
    status <- DF$status
    X <- as.matrix(DF[,-c(1:4)])    
    n <- length(unique(id))
    p <- ncol(X)
    T <- DF$Time
    mt <- aggregate(event ~ id, data = DF, sum)$event
    clsz <- mt + 1
    Y <- rep(DF$Time[event == 0], mt + 1)
    cluster <- unlist(sapply(mt + 1, function(x) 1:x))
    B <- stdErr@B
    uID <- unique(DF$ID)
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
            with(regFit(DF2, engine, NULL), c(alpha, beta))            
        })
        stopCluster(cl)
        betaMatrix <- t(out)
        convergence <- apply(betaMatrix, 1, function(x)
            1 * (x %*% x > 1e3 * with(res, c(alpha, beta) %*% c(alpha, beta))))
    } else {
        betaMatrix <- matrix(0, B, p * 2)
        convergence <- rep(0, B)
            for (i in 1:B) {
            sampled.id <- sample(unique(id), n, TRUE)
            ind <- unlist(sapply(sampled.id, function(x) which(id == x)))
            DF2 <- DF[ind,]
            DF2$id <- rep(1:n, clsz[sampled.id])
            betaMatrix[i,] <- with(regFit(DF2, engine, NULL), c(alpha, beta))
            convergence[i] <- 1 * (betaMatrix[i,] %*% betaMatrix[i,] >
                                   1e3 * with(res, c(alpha, beta) %*% c(alpha, beta)))
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
                alphaVar = betaVar[1:p, 1:p], betaVar = betaVar[1:p + p, 1:p + p],
                SEmat = betaMatrix, B = length(converged)))
}

regFit.cox.HW.resampling <- function(DF, engine, stdErr) {
    res <- regFit(DF, engine, NULL)
    id <- DF$id
    event <- DF$event
    status <- DF$status
    X <- as.matrix(DF[,-c(1:4)])    
    n <- length(unique(id))
    p <- ncol(X)
    T <- DF$Time
    mt <- aggregate(event ~ id, data = DF, sum)$event
    Y <- rep(DF$Time[event == 0], mt + 1)
    cluster <- unlist(sapply(mt + 1, function(x) 1:x))
    X <- cbind(1, X[event == 0,])
    B <- stdErr@B
    E <- matrix(rexp(n * B), nrow = n)
    Z <- matrix(rnorm((p + 1) * B), nrow = p + 1)
    Sn <- function(g, b, w, r = "both") {
        s1 <- HWeq(g, X, Y, T, cluster, mt, rep(w, mt + 1))
        if (r == "s1") return(s1)            
        lam <- npMLE(Y[event == 0], T, Y, rep(w, mt + 1))
        zHat <- mt / (lam * exp(as.matrix(X[,-1]) %*% g[-1]))
        zHat <- ifelse(zHat > 1e5, (mt + .01) / (lam * exp(as.matrix(X[,-1]) %*% g[-1]) + .01), zHat)
        zHat <- ifelse(is.na(zHat), 0, zHat)
        ## s2 <- HWeq2(b, X[,-1], Y[event == 0], status[event == 0], zHat / g[1], rep(w, mt + 1))
        s2 <- HWeq2(b, X[,-1], Y[event == 0], status[event == 0], zHat, w) ## rep(w, mt + 1))
        if (r == "s2") return(s2)
        return(c(s1, s2))
    }
    g <- c(res$muZ, res$alpha)
    V <- var(t(apply(E, 2, function(x) Sn(g, res$beta, x))))
    V1 <- V[2:(p + 1), 2:(p + 1)]
    V2 <- V[1 + p + 1:p, 1 + p + 1:p]
    lmfit1 <- t(apply(Z, 2, function(x) Sn(g + c(0, x[-1]) / sqrt(n), res$beta, rep(1, n), "s1")[-1]))
    ## V1 <- V[1:(p + 1), 1:(p + 1)]
    ## V2 <- V[1 + p + 1:p, 1 + p + 1:p]
    ## lmfit1 <- t(apply(Z, 2, function(x) Sn(g + x / sqrt(n), res$beta, rep(1, n), "s1")))
    lmfit2 <- t(apply(Z, 2, function(x) Sn(g, res$beta + x[-1] / sqrt(n), rep(1, n), "s2")))
    ## J1 <- coef(lm(sqrt(n) * lmfit1 ~ t(Z)))[-1,]
    J1 <- coef(lm(sqrt(n) * lmfit1 ~ t(Z[-1,])))[-1,]
    J2 <- coef(lm(sqrt(n) * lmfit2 ~ t(Z[-1,])))[-1,]
    if (qr(J1)$rank == (p + 1)) aVar <- solve(J1) %*% V1 %*% t(solve(J1))
    else aVar <- ginv(J1) %*% V1 %*% t(ginv(J1))
    if (qr(J2)$rank == p) bVar <- solve(J2) %*% V2 %*% t(solve(J2))
    else bVar <- ginv(J2) %*% V2 %*% t(ginv(J2))
    aSE <- sqrt(diag(aVar))
    bSE <- sqrt(diag(bVar))
    ## c(res, list(alphaSE = aSE[-1], betaSE = bSE, alphaVar = aVar[-1, -1], betaVar = bVar))
    c(res, list(alphaSE = aSE, betaSE = bSE, alphaVar = aVar, betaVar = bVar))
}

regFit.am.XCHWY.resampling <- function(DF, engine, stdErr) {
    res <- regFit(DF, engine, NULL)
    id <- DF$id
    event <- DF$event
    status <- DF$status
    X <- as.matrix(DF[,-c(1:4)])    
    n <- length(unique(id))
    p <- ncol(X)
    T <- DF$Time
    mt <- aggregate(event ~ id, data = DF, sum)$event
    Y <- rep(DF$Time[event == 0], mt + 1)
    cluster <- unlist(sapply(mt + 1, function(x) 1:x))        
    B <- stdErr@B
    E <- matrix(rexp(n * B), nrow = n)
    Z <- matrix(rnorm(p * B), nrow = p)
    Sn <- function(a, b, w, r = "both") {
        Ystar <- log(Y) + X %*% a
        Tstar <- log(T) + X %*% a
        Lam <- npMLE(Ystar[event ==  0], Tstar, Ystar, rep(w, mt + 1))
        zHat <- mt / Lam
        if (any(zHat > 1e5)) zHat <- (mt + .01) / (Lam + .01)
        ## zHat <- ifelse(zHat > 1e5, (mt + .01) / (Lam + .01), zHat)
        zHat <- ifelse(is.na(zHat), 0, zHat)
        s1 <- .C("alphaEqC", as.double(X[event == 0,]), as.double(zHat),
                 as.integer(n), as.integer(p), as.double(w),
                 res = double(p), PACKAGE = "reReg")$res / n
        if (r == "s1") return(s1)
        Lam <- npMLE(Ystar[event == 0], Tstar, Ystar)
        zHat <- mt / Lam
        if (any(zHat > 1e5)) zHat <- (mt + .01) / (Lam + .01)
        ## zHat <- ifelse(zHat > 1e5, (mt + .01) / (Lam + .01), zHat)
        zHat <- ifelse(is.na(zHat), 0, zHat)
        Ystar <- (log(Y) + X %*% b)[event == 0]
        s2 <- .C("betaEst", as.double(Ystar), as.double(X[event == 0,]),
                 as.double(status[event == 0]), as.double(zHat),
                 as.double(w), as.integer(n), as.integer(p), res = double(p),
                 PACKAGE = "reReg")$res / n
        if (r == "s2") return(s2)
        return(c(s1, s2))
    }
    V <- var(t(apply(E, 2, function(x) Sn(res$alpha, res$beta, x))))
    V1 <- V[1:p, 1:p]
    V2 <- V[1:p + p, 1:p + p]
    lmfit1 <- t(apply(Z, 2, function(x) Sn(res$alpha + x / sqrt(n), res$beta, rep(1, n), "s1")))
    lmfit2 <- t(apply(Z, 2, function(x) Sn(res$alpha, res$beta + x / sqrt(n), rep(1, n), "s2")))
    J1 <- coef(lm(sqrt(n) * lmfit1 ~ t(Z)))[-1,]
    J2 <- coef(lm(sqrt(n) * lmfit2 ~ t(Z)))[-1,]
    if (qr(J1)$rank == p) aVar <- solve(J1) %*% V1 %*% t(solve(J1))
    else aVar <- ginv(J1) %*% V1 %*% t(ginv(J1))
    if (qr(J2)$rank == p) bVar <- solve(J2) %*% V2 %*% t(solve(J2))
    else bVar <- ginv(J2) %*% V2 %*% t(ginv(J2))    
    aSE <- sqrt(diag(aVar))
    bSE <- sqrt(diag(bVar))
    c(res, list(alphaSE = aSE, betaSE = bSE, alphaVar = aVar, betaVar = bVar))
    ## c(res, list(alphaSE = aSE, betaSE = bSE, da = da, va = va, db = db, vb = vb, B = stdErr@B))
}

regFit.sc.XCYH.resampling <- function(DF, engine, stdErr) {
    if (is.na(match(engine@solver, c("dfsane", "BBsolve", "optim", "BBoptim")))) {
        print("Warning: Unidentified solver; BB::dfsane is used.")
        engine@solver <- "dfsane"
    }
    res <- regFit(DF, engine, NULL)
    DF0 <- subset(DF, event == 0)
    p <- ncol(DF) - 4
    X <- as.matrix(DF0[,-(1:4)])
    n <- nrow(DF0)
    m <- aggregate(event ~ id, data = DF, sum)[,2]
    index <- c(1, cumsum(m)[-n] + 1)
    B <- stdErr@B
    E <- matrix(rexp(n * B), nrow = n)
    Z <- matrix(rnorm((2 * p  + 1) * B), nrow = 2 * p + 1)
    a <- res$alpha
    b <- c(res$log.muZ, res$beta - res$alpha)
    U12 <- function(a, b, w) {
        yi <- log(DF0$Time) + X %*% a
        tij <- log(DF$Time) + as.matrix(DF[,-(1:4)]) %*% a
        tij <- tij[DF$event == 1]
        if (engine@eqType %in% c("Logrank", "logrank")) 
            s1 <- .C("sc1Log", as.integer(n), as.integer(p), as.integer(index - 1),
                     as.integer(m), as.double(yi), as.double(tij), as.double(X), as.double(w), 
                     result = double(p), PACKAGE = "reReg")$result
        if (engine@eqType %in% c("Gehan", "gehan"))
            s1 <- .C("sc1Gehan", as.integer(n), as.integer(p), as.integer(index - 1),
                     as.integer(m), as.double(yi), as.double(tij), as.double(X), as.double(w), 
                     result = double(p), PACKAGE = "reReg")$result / n
        t0.rate <- unique(sort(c(tij, yi)))
        rate <- .C("scRate", as.integer(n), as.integer(index - 1), as.integer(m),
                   as.integer(length(t0.rate)), as.double(w), as.double(yi),
                   as.double(tij), as.double(t0.rate), result = double(length(t0.rate)),
                   PACKAGE = "reReg")$result
        Lam <- exp(-approx(t0.rate, rate, yi)$y)
        tmp <- m / Lam
        tmp <- ifelse(tmp > 1e5, (m + .01) / (Lam + .01), tmp)
        s2 <-  .C("sc2", as.integer(n), as.integer(p + 1), as.double(cbind(1, X)),
                  as.double(exp(as.matrix(cbind(1, X)) %*% b)),
                  as.double(tmp), as.double(w),
                  result = double(p + 1), PACKAGE = "reReg")$result / n
        return(c(s1, s2))
    }
    V <- var(t(apply(E, 2, function(x) U12(a, b, x))))
    tmp <- t(apply(Z, 2, function(x) n^.5 * U12(a + x[1:p] * n^-.5, b + x[-(1:p)] * n^-.5, rep(1, n))))
    J0 <- t(coef(lm(tmp[,1:p] ~ t(Z[1:p,])))[-1,])
    Jtmp <- t(coef(lm(tmp[,-(1:p)] ~ t(Z)))[-1,])
    J <- rbind(cbind(J0, matrix(0, ncol = p + 1, nrow = nrow(J0))), cbind(Jtmp))
    if (qr(J)$rank == nrow(J)) J <- solve(J) else J <- ginv(J)
    if (qr(J0)$rank == nrow(J0)) J0 <- solve(J0) else J0 <- ginv(J0)
    ase <- J %*% V %*% t(J)    
    list(alpha = a, beta = res$beta,
         alphaSE = sqrt(diag(J0 %*% V[1:p, 1:p] %*% t(J0))),
         betaSE = sqrt(diag(ase[1:p, 1:p] + ase[(p + 2):(2 * p + 1), (p + 2):(2 * p + 1)] +
                            2 * ase[1:p, (p + 2):(2 * p + 1)])),
         gamma = b,
         varMat = ase,
         log.muZ = res$log.muZ,
         values = c(res$a.value, res$b.value),
         alphaVar = ase[1:p, 1:p],
         betaVar = ase[1:p, 1:p] + ase[(p+2):(2*p+1), (p+2):(2*p+1)] + 2 * ase[1:p, (p+2):(2*p+1)])
}

##############################################################################
# Nonparametric (~1)
##############################################################################
npFit.sc.XCYH <- function(DF, alpha, beta, engine, stdErr) {
    DF0 <- subset(DF, event == 0)
    p <- ncol(DF) - 4
    X <- as.matrix(DF0[,-(1:4)])
    yi <- log(DF0$Time) + X %*% alpha
    n <- nrow(DF0)
    tij <- log(DF$Time) + as.matrix(DF[,-(1:4)]) %*% alpha
    tij <- tij[DF$event == 1]
    m <- aggregate(event ~ id, data = DF, sum)[,2]
    index <- c(1, cumsum(m)[-n] + 1)
    t0.rate <- unique(sort(c(tij, yi))) ## log scale
    rate <- .C("scRate", as.integer(n), as.integer(index - 1), as.integer(m),
               as.integer(length(t0.rate)), as.double(rep(1, n)), as.double(yi),
               as.double(tij), as.double(t0.rate), result = double(length(t0.rate)),
               PACKAGE = "reReg")$result
    if (engine@muZ == 0) 
        engine@muZ <- sum(m / approx(t0.rate, exp(-rate), yi, yleft = 0, yright = max(rate),
                                     method = "constant")$y, na.rm = TRUE) / n
    rate <- exp(-rate) * engine@muZ
    rate0 <- approxfun(exp(t0.rate), rate, yleft = 0, yright = max(rate), method = "constant")
    list(rate0 = rate0,  rate0.lower = NULL, rate0.upper = NULL, t0.rate = exp(t0.rate),
         haz0 = NULL, haz0.lower = NULL, haz0.upper = NULL, t0.haz = NULL)
}

npFit.SE.sc.XCYH <- function(DF, alpha, beta, engine, stdErr) {
    DF0 <- subset(DF, event == 0)
    p <- ncol(DF) - 4
    X <- as.matrix(DF0[,-(1:4)])
    yi <- log(DF0$Time) + X %*% alpha
    n <- nrow(DF0)
    tij <- log(DF$Time) + as.matrix(DF[,-(1:4)]) %*% alpha
    tij <- tij[DF$event == 1]
    m <- aggregate(event ~ id, data = DF, sum)[,2]
    index <- c(1, cumsum(m)[-n] + 1)
    t0.rate <- unique(sort(c(tij, yi))) ## log scale
    rate <- .C("scRate", as.integer(n), as.integer(index - 1), as.integer(m),
               as.integer(length(t0.rate)), as.double(rep(1, n)), as.double(yi),
               as.double(tij), as.double(t0.rate), result = double(length(t0.rate)),
               PACKAGE = "reReg")$result
    if (engine@muZ == 0)
        engine@muZ <- sum(m / approx(t0.rate, exp(-rate), yi, yleft = 0, yright = max(rate),
                                     method = "constant")$y, na.rm = TRUE) / n
    rate <- exp(-rate) * engine@muZ
    rate0 <- approxfun(exp(t0.rate), rate, yleft = 0, yright = max(rate), method = "constant")
    B <- stdErr@B
    rateMat <- matrix(NA, B, length(t0.rate))
    for (i in 1:B) {
        W <- rexp(n)
        rateMat[i,] <- .C("scRate", as.integer(n), as.integer(index - 1), as.integer(m),
                          as.integer(length(t0.rate)), as.double(W), as.double(yi),
                          as.double(tij), as.double(t0.rate), result = double(length(t0.rate)),
                          PACKAGE = "reReg")$result
        rateMat[i,] <- exp(-rateMat[i,]) * engine@muZ
    }
    rl <- apply(rateMat, 2, quantile, prob = .025)
    ru <- apply(rateMat, 2, quantile, prob = .975)
    list(rate0 = rate0,
         rate0.lower = approxfun(exp(t0.rate), rl, yleft = 0, yright = max(rl), method = "constant"),
         rate0.upper = approxfun(exp(t0.rate), ru, yleft = 0, yright = max(ru), method = "constant"),
         t0.rate = exp(t0.rate),
         haz0 = NULL, haz0.lower = NULL, haz0.upper = NULL, t0.haz = NULL)
}

npFit.cox.GL <- function(DF, alpha, beta, engine, stdErr) {
    id <- DF$id
    event <- DF$event
    X <- as.matrix(DF[,-c(1:4)])
    p <- ncol(X)
    T <- DF$Time
    mt <- aggregate(event ~ id, data = DF, sum)$event
    Y <- rep(DF$Time[!event], mt + 1)
    t0.rate <- unique(sort(T))
    xb <- exp(X[!event,] %*% beta) 
    cl <- mt + 1
    rate <- .C("glCoxRate", as.double(ifelse(T == Y, 1e5, T)), as.double(Y[!event]),
               as.double(xb), as.double(engine@wgt), as.double(t0.rate), as.integer(length(t0.rate)), 
               as.integer(length(T)), as.integer(cl), as.integer(c(0, cumsum(cl)[-length(cl)])),
               as.integer(sum(!event)), out = double(length(t0.rate)), PACKAGE = "reReg")$out
    rate0 <- approxfun(t0.rate, rate, yleft = 0, yright = max(rate), method = "constant")
    fit.coxph <- coxph(Surv(Time, status) ~ ., data = subset(DF, !event, select = -c(id, event)))
    cumHaz <- basehaz(fit.coxph)
    haz0 <- with(cumHaz, approxfun(time, hazard, yleft = 0, yright = max(hazard), method = "constant"))
    list(rate0 = rate0, rate0.lower = NULL, rate0.upper = NULL, t0.rate = t0.rate,
         haz0 = haz0, haz0.lower = NULL, haz0.upper = NULL, t0.haz = cumHaz$time)
}

npFit.SE.cox.GL <- function(DF, alpha, beta, engine, stdErr) {
    id <- DF$id
    event <- DF$event
    X <- as.matrix(DF[,-c(1:4)])
    p <- ncol(X)
    T <- DF$Time
    mt <- aggregate(event ~ id, data = DF, sum)$event
    Y <- rep(DF$Time[!event], mt + 1)
    t0.rate <- unique(sort(T))
    xb <- exp(X[!event,] %*% beta) 
    cl <- mt + 1
    rate <- .C("glCoxRate", as.double(ifelse(T == Y, 1e5, T)), as.double(Y[!event]),
               as.double(xb), as.double(engine@wgt), as.double(t0.rate), as.integer(length(t0.rate)), 
               as.integer(length(T)), as.integer(cl), as.integer(c(0, cumsum(cl)[-length(cl)])),
               as.integer(sum(!event)), out = double(length(t0.rate)), PACKAGE = "reReg")$out
    rate0 <- approxfun(t0.rate, rate, yleft = 0, yright = max(rate), method = "constant")
    fit.coxph <- coxph(Surv(Time, status) ~ ., data = subset(DF, !event, select = -c(id, event)))
    cumHaz <- basehaz(fit.coxph)
    t0.haz <- cumHaz$time
    haz0 <- with(cumHaz, approxfun(time, hazard, yleft = 0, yright = max(hazard), method = "constant"))
    B <- stdErr@B
    rateMat <- matrix(NA, B, length(t0.rate))
    hazMat <- matrix(NA, B, length(cumHaz$time))
    id0 <- unique(id)
    for (i in 1:B) {
        sampled.id <- sample(id0, length(id0), TRUE)
        ind <- unlist(sapply(sampled.id, function(x) which(DF$id == x)))
        DF2 <- DF[ind,]
        DF2$id <- rep(1:length(id0), table(DF$id)[sampled.id])
        rownames(DF2) <- NULL
        fit.coxphB <- coxph(Surv(Time, status) ~ ., data = subset(DF2, !event, select = -c(id, event)))
        cumHaz <- basehaz(fit.coxphB)
        ## cumHaz$hazard <- cumHaz$hazard / max(cumHaz$hazard)
        X0 <- as.matrix(DF2[!DF2$event, -(1:4)])
        wgt <- sapply(exp(X0 %*% coef(fit.coxphB)), function(x)
            with(cumHaz, approxfun(time, exp(-hazard * x), yleft = 1, yright = min(exp(-hazard * x)),
                                   method = "constant"))(DF2$Time))
        engine@wgt <- 1 / wgt ## ifelse(wgt == 0, 1 / sort(c(wgt))[2], 1 / wgt)
        engine@wgt <- ifelse(engine@wgt > 1e5, 1e5, engine@wgt)
        tmp <- npFit.cox.GL(DF2, alpha, beta, engine, NULL)
        rateMat[i,] <- tmp$rate0(t0.rate)
        hazMat[i,] <- tmp$haz0(t0.haz)
    }
    rl <- apply(rateMat, 2, quantile, prob = .025)
    ru <- apply(rateMat, 2, quantile, prob = .975)
    hl <- apply(hazMat, 2, quantile, prob = .025)
    hu <- apply(hazMat, 2, quantile, prob = .975)
    list(rate0 = rate0,
         rate0.lower = approxfun(t0.rate, rl, yleft = 0, yright = max(rl), method = "constant"),
         rate0.upper = approxfun(t0.rate, ru, yleft = 0, yright = max(ru), method = "constant"),
         t0.rate = t0.rate,
         haz0 = haz0,
         haz0.lower = approxfun(t0.haz, hl, yleft = 0, yright = max(hl), method = "constant"),
         haz0.upper = approxfun(t0.haz, hu, yleft = 0, yright = max(hu), method = "constant"),
         t0.haz = t0.haz)
}

npFit.am.GL <- function(DF, alpha, beta, engine, stdErr) {
    alpha <- -alpha
    beta <- -beta
    DF0 <- subset(DF, event == 0)
    p <- ncol(DF0) - 4
    X <- as.matrix(DF0[,-(1:4)])
    status <- DF0$status
    n <- nrow(DF0)
    ## d <- max(X %*% (alpha - beta), 0)
    d <- max(X %*% (alpha - beta))
    tij <- log(DF$Time) - as.matrix(DF[,-(1:4)]) %*% alpha
    tij <- tij[DF$event == 1]
    yi <- log(DF0$Time) - X %*% beta ## beta or alpha?
    m <- aggregate(event ~ id, data = DF, sum)[,2]
    index <- c(1, cumsum(m)[-n] + 1)
    t0.rate <- unique(sort(tij)) # log scale
    t0.haz <- unique(sort(yi))
    rate <- .C("glRate", as.integer(n), as.integer(index - 1), as.integer(m),
               as.integer(length(t0.rate)),
               as.double(yi - d), as.double(tij), as.double(t0.rate), 
               result = double(length(t0.rate)), 
               PACKAGE = "reReg")$result
    haz <- .C("glHaz", as.integer(n), as.integer(status), as.integer(length(t0.haz)),
              as.double(yi), as.double(t0.haz), result = double(length(t0.haz)), 
              PACKAGE = "reReg")$result
    rate0 <- approxfun(exp(t0.rate), rate, yleft = 0, yright = max(rate), method = "constant")
    haz0 <- approxfun(exp(t0.haz), haz, yleft = 0, yright = max(haz), method = "constant")
    list(rate0 = rate0, rate0.lower = NULL, rate0.upper = NULL, t0.rate = exp(t0.rate),
         haz0 = haz0, haz0.lower = NULL, haz0.upper = NULL, t0.haz = exp(t0.haz))
}

npFit.SE.am.GL <- function(DF, alpha, beta, engine, stdErr) {
    id <- subset(DF, event == 0)$id
    B <- stdErr@B
    PE <- npFit.am.GL(DF, alpha, beta, engine, NULL)
    rateMat <- matrix(NA, B, length(PE$t0.rate))
    hazMat <- matrix(NA, B, length(PE$t0.haz))
    for (i in 1:B) {
        sampled.id <- sample(id, length(id), TRUE)
        ind <- unlist(sapply(sampled.id, function(x) which(DF$id == x)))
        DF2 <- DF[ind,]
        DF2$id <- rep(1:length(id), table(DF$id)[sampled.id])
        tmp <- npFit.am.GL(DF2, alpha, beta, engine, NULL)
        rateMat[i,] <- tmp$rate0(PE$t0.rate)
        hazMat[i,] <- tmp$haz0(PE$t0.haz)
    }
    rl <- apply(rateMat, 2, quantile, prob = .025)
    ru <- apply(rateMat, 2, quantile, prob = .975)
    hl <- apply(hazMat, 2, quantile, prob = .025)
    hu <- apply(hazMat, 2, quantile, prob = .975)
    list(rate0 = PE$rate0,
         rate0.lower = approxfun(PE$t0.rate, rl, yleft = 0, yright = max(rl), method = "constant"),
         rate0.upper = approxfun(PE$t0.rate, ru, yleft = 0, yright = max(ru), method = "constant"),
         haz0 = PE$haz0,
         haz0.lower = approxfun(PE$t0.haz, hl, yleft = 0, yright = max(hl), method = "constant"),
         haz0.upper = approxfun(PE$t0.haz, hu, yleft = 0, yright = max(hu), method = "constant"))
}

npFit.am.XCHWY <- function(DF, alpha, beta, engine, stdErr) {
    ly <- hy <- lyU <- lyL <- hyU <- hyL <- NULL
    id <- DF$id
    event <- DF$event
    status <- DF$status
    X <- as.matrix(DF[,-c(1:4)])    
    n <- length(unique(id))
    p <- ncol(X)
    T <- DF$Time
    mt <- aggregate(event ~ id, data = DF, sum)$event
    Y <- rep(DF$Time[event == 0], mt + 1)
    ## cluster <- unlist(sapply(mt + 1, function(x) 1:x))
    ## t0 <- seq(0, max(Y), length.out = 5 * nrow(DF))
    t0 <- sort(unique(T, Y))
    ng <- length(t0)
    Ya <- log(Y) + X %*% alpha
    Ta <- log(T) + X %*% alpha
    ## Ya <- Y * exp(X %*% alpha)
    ## Ta <- T * exp(X %*% alpha)
    lambda <- npMLE(Ya[which(event == 0)], Ta, Ya)
    ly <- npMLE(t0, Ta, Ya)
    zHat <- as.numeric(mt * max(ly) / lambda)
    zHat <- ifelse(zHat > 1e5, (mt * max(ly) + .01) / (lambda + .01), zHat)
    zHat <- ifelse(is.na(zHat), 0, zHat)
    ## zHat <- ifelse(zHat %in% c("Inf", "NA", "NaN"), 0, zHat)
    ly <- ly / max(ly)
    win.ly <- max(ly)
    muZ <- mean(zHat)
    Yb <- log(Y) + X %*% beta
    Yb <- Yb[event == 0]
    hy <- baseHaz(log(t0), Yb, zHat / muZ, status[event == 0])
    win.hy <- max(hy)
    list(rate0 = approxfun(t0, ly * muZ, yleft = 0, yright = max(ly * muZ), method = "constant"),
         rate0.lower = NULL, rate0.upper = NULL, t0.rate = t0,
         haz0 = approxfun(t0, hy, yleft = 0, yright = max(hy), method = "constant"),
         haz0.lower = NULL, haz0.upper = NULL, t0.haz = t0)
}

npFit.cox.NA <- function(DF, alpha, beta, engine, stdErr) {
    ## t0 <- seq(0, max(DF$Time), length.out = 5 * nrow(DF))
    T <- DF$Time
    id <- DF$id
    event <- DF$event
    status <- DF$status
    mt <- aggregate(event ~ id, data = DF, sum)$event
    Y <- rep(DF$Time[event == 0], mt + 1)  
    t0 <- sort(unique(T, Y))
    ng <- length(t0)
    event <- DF$event
    status <- DF$status
    ly <- hy <- lyU <- lyL <- hyU <- hyL <- NULL
    id <- DF$id
    cluster <- unlist(lapply(split(id, id), function(z) 1:length(z)))
    X <- as.matrix(DF[,-(1:4)])
    DF$T0 <- with(DF, unlist(lapply(split(Time, id), function(x) c(0, x)[1:length(x)])))
    if (!all(X == 0)) {
        tmp <- coxph(as.formula(paste("Surv(T0, Time, event)~",
                                      paste(colnames(DF)[-c(1:4, ncol(DF))], collapse = "+"))),
                     data = DF)
        base <- data.frame(matrix(0, nrow = 1, ncol = ncol(X)))
        names(base) <- names(coef(tmp))
        tmp0 <- survfit(tmp, newdata = base)
        ly <- with(tmp0, approx(time, cumhaz, t0)$y)
        lyU <- -log(with(tmp0, approx(time, upper, t0)$y))
        lyL <- -log(with(tmp0, approx(time, lower, t0)$y))
        tmp <- coxph(as.formula(paste("Surv(Time, status)~",
                                      paste(colnames(DF)[-c(1:4, ncol(DF))], collapse = "+"))),
                     data = DF[event == 0,])
        ## tmp0 <- survfit(tmp, newdata = base)
        ## hy <- with(tmp0, approx(time, cumhaz, t0)$y)
        ## hyU <- -log(with(tmp0, approx(time, upper, t0)$y))
        ## hyL <- -log(with(tmp0, approx(time, lower, t0)$y))
    }
    if (all(X == 0)) {
        tmp <- coxph(Surv(T0, Time, event) ~ 1, data = DF)
        ly <- with(basehaz(tmp), approx(time, hazard, t0)$y)
        lyU <- -log(with(survfit(tmp), approx(time, upper, t0)$y))
        lyL <- -log(with(survfit(tmp), approx(time, lower, t0)$y))
        ## tmp <- coxph(Surv(Time, status) ~ 1, data = DF[event == 0,])
        ## hy <- with(basehaz(tmp), approx(time, hazard, t0)$y)
        ## hyU <- -log(with(survfit(tmp), approx(time, upper, t0)$y))
        ## hyL <- -log(with(survfit(tmp), approx(time, lower, t0)$y))
    }
    list(rate0 = approxfun(t0, ly, yleft = 0, yright = max(ly, na.rm = TRUE), method = "constant"),
         rate0.lower = approxfun(t0, lyL, yleft = 0, yright = max(lyL, na.rm = TRUE), method = "constant"),
         rate0.upper = approxfun(t0, lyU, yleft = 0, yright = max(lyU, na.rm = TRUE), method = "constant"),
         t0.rate = t0,
         haz0 = NULL, haz0.lower = NULL, haz0.upper = NULL,
         ## haz0 = approxfun(t0, hy, yleft = 0, yright = max(hy, na.rm = TRUE), method = "constant"),
         ## haz0.lower = approxfun(t0, hyL, yleft = 0, yright = max(hyL, na.rm = TRUE), method = "constant"),
         ## haz0.upper = approxfun(t0, hyU, yleft = 0, yright = max(hyU, na.rm = TRUE), method = "constant"),
         t0.haz = t0)
}

## npFit.SE.cox.NA <- function(DF, alpha, beta, engine, stdErr) {
##     ly <- hy <- lyU <- lyL <- hyU <- hyL <- NULL
##     id <- DF$id
##     B <- stdErr@B
##     cluster <- unlist(lapply(split(id, id), function(z) 1:length(z)))
##     clsz <- unlist(lapply(split(id, id), length))
##     Y <- rep(DF$Time[cumsum(clsz)], clsz)
##     T <- DF$Time[-cumsum(clsz)]
##     nT <- length(T)
##     nY <- length(clsz)
##     tmp <- basehaz(coxph(Surv(T, rep(1, nT)) ~ 1))
##     t0 <- seq(0, max(Y), length.out = 5 * nrow(DF))
##     ng <- length(t0)
##     ly <- with(tmp, approx(time, hazard, t0)$y)
##     tmp <- basehaz(coxph(Surv(Time, event) ~ 1, data = DF[cumsum(clsz),]))    
##     hy <- with(tmp, approx(time, hazard, t0)$y)
##     bootH <- bootL <- matrix(NA, length(t0), B) 
##     for (i in 1:B) {
##         bootL[,i] <- with(basehaz(coxph(Surv(T[sample(1:nT, nT, TRUE)], rep(1, nT)) ~ 1)),
##                           approx(time, hazard, t0, yleft = min(hazard), yright = max(hazard))$y)
##         bootH[,i] <- with(basehaz(coxph(Surv(Time, event) ~ 1,
##                                         data = DF[cumsum(clsz),][sample(1:nY, nY, TRUE),])),
##                           approx(time, hazard, t0, yleft = min(hazard), yright = max(hazard))$y)
##     }
##     lyU <- apply(bootL, 1, function(z) quantile(z, .975))
##     lyL <- apply(bootL, 1, function(z) quantile(z, .025))
##     hyU <- apply(bootH, 1, function(z) quantile(z, .975))
##     hyL <- apply(bootH, 1, function(z) quantile(z, .025))
##     list(t0 = t0, ly = ly, lyU = lyU, lyL = lyL, hy = hy, hyU = hyU, hyL = hyL)
## }

npFit.cox.HW <- function(DF, alpha, beta, engine, stdErr) {
    ly <- hy <- lyU <- lyL <- hyU <- hyL <- NULL
    id <- DF$id
    T <- DF$Time
    X <- as.matrix(DF[,-c(1:4)])
    event <- DF$event
    status <- DF$status
    mt <- aggregate(event ~ id, data = DF, sum)$event
    Y <- rep(DF$Time[event == 0], mt + 1)
    if (all(X == 0)) alpha <- beta <- 0
    delta <- DF$event
    ## t0 <- seq(0, max(Y), length.out = 5 * nrow(DF))
    t0 <- sort(unique(T, Y))
    ng <- length(t0)
    Ya <- ifelse(Y <= 0, 0, log(Y))
    Ta <- ifelse(T <= 0, 0, log(T))
    lambda <- npMLE(Ya[event == 0], Ta, Ya)
    ly <- npMLE(t0, exp(Ta), exp(Ya))
    zHat <- as.numeric(mt * max(ly) / (lambda * exp(as.matrix(X[event == 0,]) %*% alpha)))
    zHat <- ifelse(zHat > 1e5, (mt * max(ly) + .01) / (lambda * exp(as.matrix(X[event == 0,]) %*% alpha) + .01), zHat)
    zHat <- ifelse(is.na(zHat), 0, zHat)
    ## zHat <- ifelse(zHat %in% c("Inf", "NA", "NaN"), 0, zHat)
    muZ <- mean(zHat)
    ly <- ly / max(ly)
    win.ly <- max(ly)
    Yb <- log(Y) ## + X %*% beta
    Yb <- Yb[event == 0]
    ## hy <- sapply(t0, function(z) baseHaz(z, exp(Yb),
    ##                                      exp(as.matrix(X[event == 0,]) %*% beta) * zHat / muZ,
    ##                                      status[event == 0]))
    hy <- baseHaz(t0, exp(Yb), exp(as.matrix(X[event == 0,]) %*% beta) * zHat / muZ, status[event == 0])
    win.hy <- max(hy)
    list(rate0 = approxfun(t0, ly * muZ, yleft = 0, yright = max(ly * muZ, na.rm = TRUE), method = "constant"),
         rate0.lower = NULL, rate0.upper = NULL, t0.rate = t0,
         haz0 = approxfun(t0, hy, yleft = 0, yright = max(hy, na.rm = TRUE), method = "constant"),
         haz0.lower = NULL, haz0.upper = NULL, t0.haz = t0)
}

npFit.SE.am.XCHWY <- function(DF, alpha, beta, engine, stdErr) {
    B <- stdErr@B
    ly <- hy <- lyU <- lyL <- hyU <- hyL <- NULL
    id <- DF$id
    event <- DF$event
    status <- DF$status
    X <- as.matrix(DF[,-c(1:4)])    
    n <- length(unique(id))
    p <- ncol(X)
    T <- DF$Time
    mt <- aggregate(event ~ id, data = DF, sum)$event
    Y <- rep(DF$Time[event == 0], mt + 1)
    cluster <- unlist(sapply(mt + 1, function(x) 1:x))
    t0 <- sort(unique(T, Y))
    ng <- length(t0)
    ## Ya <- log(Y) + X %*% alpha
    ## Ta <- log(T) + X %*% alpha
    Ya <- Y * exp(X %*% alpha)
    Ta <- T * exp(X %*% alpha)
    lambda <- npMLE(Ya[event == 0], Ta, Ya)
    ly <- npMLE(t0, Ta, Ya)
    zHat <-  as.numeric(mt * max(ly) / lambda)
    zHat <- ifelse(zHat > 1e5, (mt * max(ly) + .01) / (lambda + .01), zHat)
    zHat <- ifelse(is.na(zHat), 0, zHat)
    ## zHat <- ifelse(zHat %in% c("Inf", "NA", "NaN"), 0, zHat)
    ly <- ly / max(ly)
    E <- matrix(rexp(length(Ya) * B), nrow = length(Ya))
    ## E <- matrix(rexp(length(t0) * B), nrow = length(t0))
    lytmp <- apply(E, 2, function(x) npMLE(t0, Ta, Ya, x))
    lytmp <- apply(lytmp, 2, function(z) z / max(z))
    lyU <- apply(lytmp, 1, function(z) quantile(z, 0.975))
    lyL <- apply(lytmp, 1, function(z) quantile(z, 0.025))
    Yb <- log(Y) + X %*% beta
    Yb <- Yb[event == 0]
    muZ <- mean(zHat)
    ## hy <- sapply(t0, function(z) baseHaz(z, exp(Yb), zHat / muZ, status[event == 0]))
    hy <- baseHaz(t0, exp(Yb), zHat / muZ, status[event == 0])
    E <- matrix(rexp(n * B), nrow = n)
    hytmp <- apply(E, 2, function(z) baseHaz(t0, exp(Yb), zHat / muZ, status[event == 0], z))
    hyU <- apply(hytmp, 1, function(z) quantile(z, 0.975))
    hyL <- apply(hytmp, 1, function(z) quantile(z, 0.025))
    list(rate0 = approxfun(t0, ly * muZ, yleft = 0, yright = max(ly * muZ, na.rm = TRUE), method = "constant"),
         rate0.lower = approxfun(t0, lyL * muZ,
                                 yleft = 0, yright = max(lyL * muZ, na.rm = TRUE), method = "constant"),
         rate0.upper = approxfun(t0, lyU * muZ,
                                 yleft = 0, yright = max(lyU * muZ, na.rm = TRUE), method = "constant"),
         t0.rate = t0,
         haz0 = approxfun(t0, hy, yleft = 0, yright = max(hy, na.rm = TRUE), method = "constant"),
         haz0.lower = approxfun(t0, hyL,
                                yleft = 0, yright = max(hyL, na.rm = TRUE), method = "constant"),
         haz0.upper = approxfun(t0, hyU,
                                yleft = 0, yright = max(hyU, na.rm = TRUE), method = "constant"),
         t0.haz = t0)
}

npFit.SE.cox.HW <- function(DF, alpha, beta, engine, stdErr) {
    B <- stdErr@B
    ly <- hy <- lyU <- lyL <- hyU <- hyL <- NULL
    id <- DF$id
    event <- DF$event
    status <- DF$status
    X <- as.matrix(DF[,-c(1:4)])    
    n <- length(unique(id))
    p <- ncol(X)
    T <- DF$Time
    mt <- aggregate(event ~ id, data = DF, sum)$event
    Y <- rep(DF$Time[event == 0], mt + 1)
    cluster <- unlist(sapply(mt + 1, function(x) 1:x))
    if (all(X == 0)) alpha <- beta <- 0
    delta <- DF$event
    t0 <- sort(unique(T, Y))
    ng <- length(t0)
    Ya <- log(Y)
    Ta <- log(T)
    lambda <- npMLE(Ya[event == 0], Ta, Ya)
    ly <- npMLE(t0, exp(Ta), exp(Ya))
    zHat <- as.numeric(mt * max(ly) / (lambda * exp(as.matrix(X[event == 0,]) %*% alpha)))
    zHat <- ifelse(zHat > 1e5, (mt * max(ly) + .01) / (lambda * exp(as.matrix(X[event == 0,]) %*% alpha) + .01), zHat)
    zHat <- ifelse(is.na(zHat), 0, zHat)
    ## zHat <- ifelse(zHat %in% c("Inf", "NA", "NaN"), 0, zHat)
    ly <- ly / max(ly)
    E <- matrix(rexp(ng * B), nrow = ng)
    lytmp <- apply(E, 2, function(x) npMLE(t0, exp(Ta), exp(Ya), x))
    lytmp <- apply(lytmp, 2, function(z) z / max(z))
    lyU <- apply(lytmp, 1, function(z) quantile(z, 0.975))
    lyL <- apply(lytmp, 1, function(z) quantile(z, 0.025))
    Yb <- log(Y) ## + X %*% beta
    Yb <- Yb[event == 0]
    muZ <- mean(zHat)
    hy <- baseHaz(t0, exp(Yb), exp(as.matrix(X[event == 0,]) %*% beta) * zHat / muZ, status[event == 0])
    E <- matrix(rexp(n * B), nrow = n)
    hytmp <- apply(E, 2, function(z) baseHaz(t0, exp(Yb), exp(as.matrix(X[event == 0,]) %*% beta) * zHat / muZ, status[event == 0], z))
    hyU <- apply(hytmp, 1, function(z) quantile(z, 0.975))
    hyL <- apply(hytmp, 1, function(z) quantile(z, 0.025))
    list(rate0 = approxfun(t0, ly * muZ, yleft = 0, yright = max(ly * muZ, na.rm = TRUE), method = "constant"),
         rate0.lower = approxfun(t0, lyL * muZ,
                                 yleft = 0, yright = max(lyL * muZ, na.rm = TRUE), method = "constant"),
         rate0.upper = approxfun(t0, lyU * muZ,
                                 yleft = 0, yright = max(lyU * muZ, na.rm = TRUE), method = "constant"),
         t0.rate = t0,
         haz0 = approxfun(t0, hy, yleft = 0, yright = max(hy, na.rm = TRUE), method = "constant"),
         haz0.lower = approxfun(t0, hyL, yleft = 0, yright = max(hyL, na.rm = TRUE), method = "constant"),
         haz0.upper = approxfun(t0, hyU, yleft = 0, yright = max(hyU, na.rm = TRUE), method = "constant"),
         t0.haz = t0)
}

##############################################################################
# Class Definition
##############################################################################

setClass("Engine",
         representation(tol = "numeric", a0 = "numeric", b0 = "numeric", solver = "character"),
         prototype(tol = 1e-7, a0 = 0, b0 = 0, solver = "dfsane"),
         contains = "VIRTUAL")
setClass("cox.LWYY", contains = "Engine")
setClass("cox.HW", contains = "Engine")
setClass("am.XCHWY", contains = "Engine")
setClass("am.GL", contains = "Engine")
setClass("sc.XCYH",
         representation(eqType = "character", muZ = "numeric"),
         prototype(eqType = "Logrank", muZ = 0),
         contains = "Engine")
setClass("cox.GL",
         representation(wgt = "matrix"), prototype(wgt = matrix(0)), contains = "Engine")

setClass("stdErr",
         representation(B = "numeric", parallel = "logical", parCL = "integer"),
         prototype(B = 100, parallel = FALSE, parCl = parallel::detectCores() / 2),
         contains = "VIRTUAL")

setClass("bootstrap", contains="stdErr")
setClass("resampling", contains="stdErr")


##############################################################################
# Method Dispatch
##############################################################################
setGeneric("regFit", function(DF, engine, stdErr) {standardGeneric("regFit")})

setMethod("regFit", signature(engine = "cox.LWYY", stdErr = "NULL"), regFit.cox.LWYY)
setMethod("regFit", signature(engine = "cox.LWYY", stdErr = "bootstrap"), regFit.cox.LWYY)
setMethod("regFit", signature(engine = "cox.GL", stdErr = "NULL"), regFit.cox.GL)
setMethod("regFit", signature(engine = "cox.HW", stdErr = "NULL"), regFit.cox.HW)
setMethod("regFit", signature(engine = "am.XCHWY", stdErr = "NULL"), regFit.am.XCHWY)
setMethod("regFit", signature(engine = "am.GL", stdErr = "NULL"), regFit.am.GL)
setMethod("regFit", signature(engine = "sc.XCYH", stdErr = "NULL"), regFit.sc.XCYH)
setMethod("regFit", signature(engine = "Engine", stdErr = "bootstrap"),
          regFit.Engine.Bootstrap)
setMethod("regFit", signature(engine = "cox.HW", stdErr = "resampling"),
          regFit.cox.HW.resampling)
setMethod("regFit", signature(engine = "am.GL", stdErr = "resampling"),
          regFit.am.GL.resampling)
setMethod("regFit", signature(engine = "am.XCHWY", stdErr = "resampling"),
          regFit.am.XCHWY.resampling)
setMethod("regFit", signature(engine = "sc.XCYH", stdErr = "resampling"),
          regFit.sc.XCYH.resampling)

setGeneric("npFit", function(DF, alpha, beta, engine, stdErr) {standardGeneric("npFit")})
setMethod("npFit", signature(engine = "cox.LWYY", stdErr = "NULL"), npFit.cox.NA)
setMethod("npFit", signature(engine = "cox.GL", stdErr = "NULL"), npFit.cox.GL)
setMethod("npFit", signature(engine = "cox.HW", stdErr = "NULL"), npFit.cox.HW)
setMethod("npFit", signature(engine = "am.XCHWY", stdErr = "NULL"), npFit.am.XCHWY)
setMethod("npFit", signature(engine = "cox.LWYY", stdErr = "bootstrap"), npFit.cox.NA)
setMethod("npFit", signature(engine = "cox.HW", stdErr = "bootstrap"), npFit.SE.cox.HW)
setMethod("npFit", signature(engine = "cox.GL", stdErr = "bootstrap"), npFit.SE.cox.GL)
setMethod("npFit", signature(engine = "cox.HW", stdErr = "resampling"), npFit.SE.cox.HW)
setMethod("npFit", signature(engine = "am.XCHWY", stdErr = "resampling"), npFit.SE.am.XCHWY)
setMethod("npFit", signature(engine = "am.XCHWY", stdErr = "bootstrap"), npFit.SE.am.XCHWY)
setMethod("npFit", signature(engine = "am.GL", stdErr = "bootstrap"), npFit.SE.am.GL)
setMethod("npFit", signature(engine = "am.GL", stdErr = "resampling"), npFit.SE.am.GL)
setMethod("npFit", signature(engine = "am.GL", stdErr = "NULL"), npFit.am.GL)
setMethod("npFit", signature(engine = "sc.XCYH", stdErr = "bootstrap"), npFit.SE.sc.XCYH)
setMethod("npFit", signature(engine = "sc.XCYH", stdErr = "resampling"), npFit.SE.sc.XCYH)
setMethod("npFit", signature(engine = "sc.XCYH", stdErr = "NULL"), npFit.sc.XCYH)

## #' When a joint model is fitted (e.g., \code{method = "cox.HW"} or \code{method = "am.XCHWY"}),
## #' the hazard function of the terminal event is either in a Cox model or an accelerated failure time model.

#' Fits Semiparametric Regression Models for Recurrent Event Data
#'
#' Fits a semiparametric regression model for the recurrent event data.
#' The rate function of the underlying process for the recurrent event process
#' can be specified as a Cox-type model, an accelerated mean model, or a generalized scale-change model.
#' See details for model specifications.
#'
#' Suppose the recurrent event process and the failure events are observed in the time interval \eqn{t\in[0,\tau]},
#' for some constant \eqn{\tau}.
#' We formulate the rate function, \eqn{\lambda(t)}, for the recurrent event process and
#' the hazard function, \eqn{h(t)}, for the censoring time
#' under the following model specifications:
#' \describe{
#'   \item{Cox-type model:}{
#' \deqn{\lambda(t) = Z \lambda_0(t) e^{X^\top\alpha}, h(t) = Z h_0(t)e^{X^\top\beta},}}
#'   \item{Accelerated mean model:}{
#' \deqn{\lambda(t) = Z \lambda_0(te^{X^\top\alpha})e^{X^\top\alpha}, h(t) = Z h_0(te^{X^\top\beta})e^{X^\top\beta},}}
#'   \item{Scale-change model:}{
#' \deqn{\lambda(t) = Z \lambda_0(te^{X^\top\alpha})e^{X^\top\beta},}}
#' }
#' where \eqn{\lambda_0(t)} is the baseline rate function, \eqn{h_0(t)} is the baseline hazard function,
#' \eqn{X} is a \eqn{n} by \eqn{p} covariate matrix and \eqn{\alpha},
#' \eqn{Z} is an unobserved shared frailty variable,
#' and \eqn{\beta} are unknown \eqn{p}-dimensional regression parameters.
#'
#' The \code{reReg} function fits models with the following available methods:
#' \describe{
#'   \item{\code{method = "cox.LWYY"}}{
#' assumes the Cox-type model with \code{Z = 1} and requires independent censoring. 
#' The returned result is equivalent to that from \code{coxph}. See reference Lin et al. (2000).}
#'   \item{\code{method = "cox.HW"}}{
#' assumes the Cox-type model with unspecified \code{Z}, thus accommodate informative censoring.
#' See the references See reference Wang, Qin and Chiang (2001) and Huang and Wang (2004).}
#'   \item{\code{method = "am.GL"}}{
#' assumes the accelerated mean model with \code{Z = 1} and requires independent censoring. 
#' See the reference Ghosh and Lin (2003).}
#'   \item{\code{method = "am.XCHWY"}}{
#' assumes the accelerated mean model with unspecified \code{Z}, thus accommodate informative censoring.
#' See the reference Xu et al. (2017).}
#'   \item{\code{method = "sc.XCYH"}}{
#' assumes the generalized scale-change model, and includes the methods \code{"cox.HW"} and \code{"am.XCHWY"} as special cases.
#' Informative censoring is accounted for through the unspecified frailty variable \code{Z}.
#' The methods also provide a hypothesis test of these submodels.}
#' }
#'
#' The available methods for variance estimation are:
#' \describe{
#'   \item{\code{NULL}}{variance estimation will not be performed. This is equivalent to setting \code{B = 0}.}
#'   \item{\code{"resampling"}}{performs the efficient resampling-based sandwich estimator that works with methods \code{"cox.HW"}, \code{"am.XCHWY"} and \code{"sc.XCYH"}.}
#'   \item{\code{"bootstrap"}}{works with all fitting methods.}
#' }
#'
#' The \code{control} list consists of the following parameters:
#' \describe{
#'   \item{\code{tol}}{absolute error tolerance.}
#'   \item{\code{a0, b0}}{initial guesses used for root search.}
#'   \item{\code{solver}}{the equation solver used for root search.
#' The available options are \code{BB::BBsolve}, \code{BB::dfsane}, \code{BB:BBoptim}, and \code{optim}.}
#'   \item{\code{parallel}}{an logical value indicating whether parallel computation will be applied when \code{se = "bootstrap"} is called.}
#'   \item{\code{parCl}}{an integer value specifying the number of CPU cores to be used when \code{parallel = TRUE}.
#' The default value is half the CPU cores on the current host.}
#' }
#' 
#' @param formula a formula object, with the response on the left of a "~" operator, and the predictors on the right.
#' The response must be a recurrent event survival object as returned by function \code{reSurv}.
#' @param data  an optional data frame in which to interpret the variables occurring in the \code{"formula"}.
#' @param B a numeric value specifies the number of resampling for variance estimation.
#' When \code{B = 0}, variance estimation will not be performed.
#' @param method a character string specifying the underlying model. See \bold{Details}.
#' @param se a character string specifying the method for standard error estimation. See \bold{Details}.
#' @param contrasts an optional list.
#' @param control a list of control parameters.
#'
#' @export
#' @references Xu, G., Chiou, S.H., Huang, C.-Y., Wang, M.-C. and Yan, J. (2017). Joint Scale-change Models for Recurrent Events and Failure Time.
#' \emph{Journal of the American Statistical Association}, \bold{112}(518): 796--805.
#' @references Lin, D., Wei, L., Yang, I. and Ying, Z. (2000). Semiparametric Regression for the Mean and Rate Functions of Recurrent Events.
#' \emph{Journal of the Royal Statistical Society: Series B (Methodological)}, \bold{62}: 711--730.
#' @references Wang, M.-C., Qin, J., and Chiang, C.-T. (2001). Analyzing Recurrent Event Data with Informative Censoring.
#' \emph{Journal of the American Statistical Association}, \bold{96}(455): 1057--1065.
#' @references Ghosh, D. and Lin, D.Y. (2003). Semiparametric Analysis of Recurrent Events Data in the Presence of Dependent Censoring.
#' \emph{Biometrics}, \bold{59}: 877--885.
#' @references Huang, C.-Y. and Wang, M.-C. (2004). Joint Modeling and Estimation for Recurrent Event Processes and Failure Time Data.
#' \emph{Journal of the American Statistical Association}, \bold{99}(468): 1153--1165.
#'
#' @importFrom stats approxfun optim
#' 
#' @seealso \code{\link{reSurv}}, \code{\link{simDat}}
#'
#' @examples
#' ## readmission data
#' data(readmission, package = "frailtypack")
#' set.seed(123)
#' ## Accelerated Mean Model
#' (fit <- reReg(reSurv(t.stop, id, event, death) ~ sex + chemo,
#'               data = subset(readmission, id < 50),
#'               method = "am.XCHWY", se = "resampling", B = 20))
#' summary(fit)
#'
#' ## Generalized Scale-Change Model
#' set.seed(123)
#' (fit <- reReg(reSurv(t.stop, id, event, death) ~ sex + chemo,
#'               data = subset(readmission, id < 50),
#'               method = "sc.XCYH", se = "resampling", B = 20))
#' summary(fit)
reReg <- function(formula, data, B = 200, 
                  method = c("cox.LWYY", "cox.GL", "cox.HW", "am.GL", "am.XCHWY", "sc.XCYH"),
                  se = c("NULL", "bootstrap", "resampling"), 
                  contrasts = NULL, control = list()) {
    method <- match.arg(method)
    se <- match.arg(se)
    Call <- match.call()
    if (missing(data)) obj <- eval(formula[[2]], parent.frame()) 
    if (!missing(data)) obj <- eval(formula[[2]], data) 
    if (!is.reSurv(obj)) stop("Response must be a reSurv object")
    formula[[2]] <- NULL
    if (formula == ~ 1) {
        DF <- cbind(obj$reDF[, -5], zero=0)
    } else {
        ## remove intercept
        if (!missing(data)) DF <- cbind(obj$reDF[,-5], model.matrix(formula, data))
        if (missing(data)) DF <- cbind(obj$reDF[,-5], model.matrix(formula, parent.frame()))
        DF <- DF[,-which(colnames(DF) == "(Intercept)")]
    }
    DF <- DF[order(DF$id, DF$Time), ]
    ## reset ID
    DF$id <- rep(1:length(unique(DF$id)), table(DF$id))
    engine.control <- control[names(control) %in% names(attr(getClass(method), "slots"))]
    engine <- do.call("new", c(list(Class = method), engine.control))
    if (se == "NULL")
        stdErr <- NULL
    else {
        stdErr.control <- control[names(control) %in% names(attr(getClass(se), "slots"))]
        stdErr <- do.call("new", c(list(Class = se), stdErr.control))
        stdErr@B <- B
    }
    p <- ncol(DF) - 4
    if (length(engine@a0) == 1) engine@a0 <- rep(engine@a0, p)
    if (length(engine@b0) == 1) engine@b0 <- rep(engine@b0, p)
    if (method %in% "sc.XCYH") {
        if (length(engine@b0) == 1) engine@b0 <- rep(engine@b0, p + 1)
        if (length(engine@b0) == p) engine@b0 <- c(0, engine@b0)
    }
    if (formula == ~1) {
        fit <- NULL
        fit$alpha <- fit$beta <- rep(NA, p)
        fit <- c(fit, npFit(DF = DF, alpha = 0, beta = 0, engine = engine, stdErr = stdErr))
    } else {
        fit <- regFit(DF = DF, engine = engine, stdErr = stdErr)
        if (method == "sc.XCYH") engine@muZ <- exp(fit$log.muZ)
        if (method == "cox.GL") engine@wgt <- fit$wgt
        fit <- c(fit, npFit(DF = DF, alpha = fit$alpha, beta = fit$beta, engine = engine, stdErr = stdErr))
    }
    class(fit) <- "reReg"
    fit$reTb <- obj$reTb
    fit$DF <- DF
    fit$call <- Call
    fit$varNames <- names(DF)[-(1:4)]
    fit$method <- method
    fit$se <- se
    if (!is.null(fit$alphaVar)) rownames(fit$alphaVar) <- colnames(fit$alphaVar) <- fit$varNames
    if (!is.null(fit$betaVar)) rownames(fit$betaVar) <- colnames(fit$betaVar) <- fit$varNames
    if (!is.null(fit$alphaSE)) names(fit$alphaSE) <- fit$varNames
    if (!is.null(fit$betaSE)) names(fit$betaSE) <- fit$varNames
    if (!is.null(fit$alpha)) names(fit$alpha) <- fit$varNames
    if (!is.null(fit$beta)) names(fit$beta) <- fit$varNames
    fit
}

#' equation wrapper
#'
#' @noRd
#' @keywords internal
eqSolve <- function(par, fn, solver, ...) {
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

npMLE <- function(t, tij, yi, weights = NULL) {
    if (is.null(weights)) weights <- rep(1, length(yi))
    ttmp <- tij[tij != yi]
    ord <- order(ttmp)
    sl <- unique(ttmp[ord])
    l <- ifelse(min(t) < max(sl), which(sl > min(t))[1], length(sl))
    tmp <- sl[l:length(sl)]
    tmp <- rev(tmp)
    tij <- rev(tij)
    yi <- rev(yi)
    ## print(length(weights))
    res <- .C("plLambda", as.double(tmp), as.double(tij), as.double(yi), as.double(weights), 
              as.integer(length(tmp)), as.integer(length(yi)),
              out = double(length(tmp)),
              PACKAGE = "reReg")$out
    out <- rev(res)[sapply(t, function(x) which(rev(tmp) >= x)[1])]
    out <- ifelse(is.na(out), 0, out)
    out <- exp(-out)
    return(out)
}

baseHaz <- function(t0, Y, zhat, delta, weights  = NULL) {
    if (is.null(weights)) weights <- rep(1, length(Y))
    .C("hwHaz", as.double(t0), as.double(Y), as.double(zhat), as.double(delta),
       as.double(weights), as.integer(length(Y)), as.integer(length(t0)), 
       out = double(length(t0)), PACKAGE = "reReg")$out
}

alphaEq <- function(alpha, X, Y, T, cluster, mt, weights = NULL) {
    n <- sum(cluster == 1)
    if (is.null(weights)) weights <- rep(1, n)
    p <- ncol(X)
    Ystar <- log(Y) + X %*% alpha
    Tstar <- log(T) + X %*% alpha
    Lambda <- npMLE(Ystar[which(cluster == 1)], Tstar, Ystar, weights = rep(weights, mt + 1))
    zHat <- mt / Lambda
    if (any(zHat > 1e5)) zHat <- (mt + .01) / (Lambda + .01)
    ## zHat <- ifelse(zHat > 1e5, (mt + .01) / (Lambda + .01), zHat)
    zHat <- ifelse(is.na(zHat), 0, zHat)
    res <- .C("alphaEqC", as.double(X[which(cluster == 1), ]), as.double(zHat),
              as.integer(n), as.integer(p), as.double(weights), 
              out = double(p), PACKAGE = "reReg")$out
    res / n
}

betaEq <- function(X, Y, T, cluster, mt, delta, zHat = NULL, alpha, beta, weights = NULL) {
    p <- ncol(X)
    n <- sum(cluster == 1)
    if (is.null(weights)) weights <- rep(1, n)
    if (is.null(zHat)) {
        Ystar <- Y * exp(X %*% alpha)
        Tstar <- T * exp(X %*% alpha)
        ## Ystar <- log(Y) + X %*% alpha
        ## Tstar <- log(T) + X %*% alpha
        lambda <- npMLE(Ystar[which(cluster == 1)], Tstar, Ystar, weights = rep(weights, mt + 1))
        ## weights = rep(weights, diff(c(which(cluster ==1), length(cluster)+1))))
        ## zHat <- as.numeric(weights * mt / lambda)
        zHat <- weights * mt / lambda
        if (any(zHat > 1e5)) zHat <- (mt + .01) / (lambda + .01)
        zHat <- ifelse(zHat > 1e5, (mt + .01) / (lambda + .01), zHat)
        zHat <- ifelse(is.na(zHat), 0, zHat)
        ## zHat <- ifelse(zHat %in% c("Inf", "NA", "NaN"), 0, zHat)
    }
    Y <- log(Y) + X %*% beta
    Y <- Y[which(cluster == 1)]
    X <- X[which(cluster == 1), ]
    res <- .C("betaEst", as.double(Y), as.double(X), as.double(delta), as.double(zHat),
              as.double(weights), as.integer(n), as.integer(p),
              out = double(p), PACKAGE = "reReg")$out
    res / n
}

LWYYeq <- function(beta, X, Y, T, cl) {
    p <- ncol(X)
    res <- vector("double", p)
    .C("lwyy", as.double(T), as.double(Y), as.double(X), as.double(wgt), as.integer(cl),
       as.integer(c(0, cumsum(cl)[-length(cl)])), as.integer(nrow(X)), as.integer(p),        
       out = double(p), PACKAGE = "reReg")$out       
}

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
                    
HWeq <-function(gamma, X, Y, T, cluster, mt, weights = NULL) {
    if (is.null(weights)) weights <- rep(1, length(T))
    n <- sum(cluster == 1)
    Lambda <- npMLE(Y[cluster == 1], T, Y, weights)
    zHat <- mt / Lambda
    zHat <- ifelse(zHat > 1e5, (mt + .01) / (Lambda + .01), zHat)
    zHat <- ifelse(is.na(zHat), 0, zHat) 
    p <- ncol(X)
    xr <- exp(X %*% gamma)
    .C("sarm1", as.double(X), as.double(weights), as.double(xr), 
       as.integer(zHat), as.integer(n), as.integer(p),
       res = double(p), PACKAGE = "reReg")$res / n
}

HWeq2 <-function(beta, X, Y, delta, zHat, weights = NULL) {
    if (is.null(weights)) weights <- rep(1, length(Y))        
    n <- nrow(X)
    p <- ncol(X)
    res <- .C("HWb", as.double(Y), as.double(X), as.double(delta), as.double(zHat),
              as.double(X %*% beta), as.double(weights), as.integer(n), as.integer(p), 
              out = double(p), PACKAGE = "reReg")$out
    res / n
}

varOut <- function(dat, na.rm = TRUE) {
    dat[which(dat %in% boxplot(dat, plot = FALSE)$out)] <- NA
    dat <- dat[complete.cases(dat),]
    var(dat, na.rm = na.rm)
}
