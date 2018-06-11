##############################################################################
## Functions for different methods
## stdErr is estimated with resampling if method = sc or am.xc,
##        bootstrap otherwise
##############################################################################

doREFit.am.XCHWY <- function(DF, engine, stdErr) {
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
    if (engine@solver %in% c("dfsane", "BBsolve")) 
        suppressWarnings(
            outA <- do.call(
                engine@solver, list(par = engine@a0, fn = alphaEq,
                                    X = X, Y = Y, T = T, cluster = cluster, mt = mt,
                                    weights = NULL, alertConvergence = FALSE, quiet = TRUE,
                                    control = list(trace = FALSE))))
    if (engine@solver == "BBoptim")
        suppressWarnings(
            outA <- do.call(
                engine@solver, list(par = engine@a0, fn = function(x)
                    sum(alphaEq(x, X = X, Y = Y, T = T, cluster = cluster, mt = mt, weights = NULL)^2),
                    quiet = TRUE, control = list(trace = FALSE))))
    if (engine@solver == "optim")
        suppressWarnings(
            outA <- do.call(
                engine@solver, list(par = engine@a0, fn = function(x)
                    sum(alphaEq(x, X = X, Y = Y, T = T, cluster = cluster, mt = mt, weights = NULL)^2),
                    control = list(trace = FALSE))))
    alpha <- outA$par
    Ystar <- log(Y) + X %*% alpha
    Tstar <- log(T) + X %*% alpha
    lambda <- npMLE(Ystar[event == 0], Tstar, Ystar)
    zHat <- as.numeric(mt * npMLE(log(max(Y)), Tstar, Ystar) / lambda)
    zHat <- ifelse(zHat %in% c("Inf", "NA", "NaN"), 0, zHat)
    if (engine@solver %in% c("dfsane", "BBsolve")) 
        suppressWarnings(
            outB <- do.call(
                engine@solver, list(par = engine@b0, fn = betaEq,
                                    X = X, Y = Y, T = T, cluster = cluster,                                    
                                    delta = status[event == 0], mt = mt,
                                    alpha = outA$par, zHat = zHat, weights = NULL,
                                    quiet = TRUE, control = list(trace = FALSE))))
    if (engine@solver == "BBoptim")
        suppressWarnings(
            outB <- do.call(
                engine@solver, list(par = engine@b0, fn = function(x)
                    sum(betaEq(x, X = X, Y = Y, T = T, cluster = cluster,
                               delta = status[event == 0], mt = mt,
                               alpha = outA$par, zHat = zHat, weights = NULL)^2),
                    quiet = TRUE, control = list(trace = FALSE))))
    if (engine@solver == "optim") 
        suppressWarnings(
            outB <- do.call(
                engine@solver, list(par = engine@b0, fn = function(x)
                    sum(betaEq(x, X = X, Y = Y, T = T, cluster = cluster,
                               delta = status[event == 0], mt = mt,
                               alpha = outA$par, zHat = zHat, weights = NULL)^2),
                    control = list(trace = FALSE))))
    list(alpha = outA$par, aconv = outA$convergence, beta = outB$par, bconv = outB$convergence,
         muZ = mean(zHat), zHat = zHat)
}

doREFit.am.GL <- function(DF, engine, stdErr) {
    DF0 <- subset(DF, event == 0)
    p <- ncol(DF0) - 4
    alpha <- beta <- gamma <- rep(0, p)
    aSE <- bSE <- da <- va <- db <- vb <- NA
    Y <- log(DF0$Time)
    X <- as.matrix(DF0[,-(1:4)])
    status <- DF0$status
    n <- nrow(DF0)
    ## Obtaining AFT estimator first
    log.est <- function(b) {
        .C("log_ns_est", as.double(b), as.double(Y), as.double(X), as.double(status),
           as.integer(rep(1, n)), as.integer(n), as.integer(p), as.integer(n),
           as.double(rep(1, n)), as.double(rep(1, n)), 
           double(p), PACKAGE = "reReg")[[11]]
    }
    if (engine@solver %in% c("dfsane", "BBsolve")) {    
        suppressWarnings(
            fit.b <- do.call(
                engine@solver, list(par = engine@b0, fn = log.est, quiet = TRUE,
                                    control = list(trace = FALSE))))
    }
    if (engine@solver == "BBoptim") {
        suppressWarnings(
            fit.b <- do.call(
                engine@solver, list(par = engine@b0, fn = function(x)
                    sum(log.est(x)^2), quiet = TRUE, control = list(trace = FALSE))))
    }
    if (engine@solver == "optim") {
        suppressWarnings(
            fit.b <- do.call(
                engine@solver, list(par = engine@b0, fn = function(x)
                    sum(log.est(x)^2), control = list(trace = FALSE))))
    }
    bhat <- fit.b$par
    m <- aggregate(event ~ id, data = DF, sum)[,2]
    index <- c(1, cumsum(m)[-n] + 1)
    ghoshU2 <- function(a) {
        d <- max(X %*% (a - bhat), 0)
        tij <- log(DF$Time) - as.matrix(DF[,-(1:4)]) %*% a
        tij <- tij[DF$event == 1]
        yi <- Y - X %*% bhat - d
        .C("glU2", as.integer(n), as.integer(p), as.integer(index - 1), as.integer(m),
           as.double(yi), as.double(tij), as.double(X), result = double(p),
           PACKAGE = "reReg")$result
    }
    if (engine@solver %in% c("dfsane", "BBsolve")) {
        suppressWarnings(
            fit.a <- do.call(
                engine@solver, list(par = engine@a0, fn = ghoshU2, quiet = TRUE,
                                    control = list(trace = FALSE))))
    }
    if (engine@solver == "BBoptim") {
        suppressWarnings(
            fit.a <- do.call(
                engine@solver, list(par = engine@a0, fn = function(x)
                    sum(log.est(x)^2), quiet = TRUE, control = list(trace = FALSE))))
    }
    if (engine@solver == "optim"){
        suppressWarnings(
            fit.a <- do.call(
                engine@solver, list(par = engine@a0, fn = function(x)
                    sum(log.est(x)^2), control = list(trace = FALSE))))
    }
    fit.b$par <- -fit.b$par
    fit.a$par <- -fit.a$par
    list(alpha = fit.a$par, aconv = fit.a$convergence,
         beta = fit.b$par, bconv = fit.b$convergence, muZ = NA)
}

doREFit.cox.HW <- function(DF, engine, stdErr) {
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
    aSE <- bSE <- da <- va <- db <- vb <- NA
    X <- cbind(1, X[event == 0,])
    ## outA <- dfsane(gamma, HWeq, X = X, Y = Y, T = T, cluster = cluster, mt = mt,
    ##                alertConvergence = FALSE,
    ##                quiet = TRUE, control = list(NM = FALSE, M = 100, noimp = 50, trace = FALSE))
    outA <- BBsolve(engine@b0, HWeq, X = X, Y = Y, T = T, cluster = cluster, mt = mt, quiet = TRUE)
    alpha <- outA$par <- outA$par[-1]
    muZ <- outA$par[1]
    lambda <- npMLE(Y[event == 0], T, Y)
    ## zHat <- as.numeric(mt * npMLE(max(Y), T, Y) / (lambda * exp(X[, -1] %*% alpha)))
    zHat <- as.numeric(mt / (lambda * exp(as.matrix(X[, -1]) %*% alpha)))
    zHat <- ifelse(zHat %in% c("Inf", "NA", "NaN"), 0, zHat)
    ## muZ <- mean(zHat)
    outB <- dfsane(engine@a0, HWeq2, X = as.matrix(X[,-1]), Y = Y[event == 0],
                   delta = status[event == 0], zHat = zHat/muZ,
                   alertConvergence = FALSE, quiet = TRUE,
                   control = list(NM = FALSE, M = 100, noimp = 50, trace = FALSE))
    list(alpha = outA$par, aconv = outA$convergence,
         beta = outB$par, bconv = outB$convergence, muZ = muZ)
}

doREFit.cox.LWYY <- function(DF, engine, stdErr) {
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
    out <- dfsane(par = engine@a0, fn = LWYYeq, X = as.matrix(X[event == 0, ]),
                  Y = Y[event == 0], T = ifelse(T == Y, 1e5, T), cl = mt + 1,
                  ## cl = unlist(lapply(split(id, id), length)), 
                  alertConvergence = FALSE, quiet = TRUE,
                  control = list(NM = FALSE, M = 100, noimp = 50, trace = FALSE))
    list(alpha = out$par, beta = rep(0, p), muZ = NA)
}

doREFit.sc.XCYH <- function(DF, engine, stdErr) {
    id <- DF$id
    clsz <- unlist(lapply(split(id, id), length))
    ind <- cumsum(clsz)
    df <- DF
    m <- rep(clsz, clsz) - 1
    y <- rep(df$Time[ind], clsz)
    rmv <- sort(unique(rep(ind, clsz) * (m > 0)))[-1]
    df$Time <- df$Time * (m > 0)
    df <- df[-rmv, ]
    m <- m[-rmv]
    y <- y[-rmv]
    if (is.na(match(engine@solver, c("dfsane", "BBsolve", "optim", "BBoptim")))) {
        print("Warning: Unidentified solver; BB::dfsane is used.")
        engine@solver <- "dfsane"
    }
    out <- with(df, sarmRV(id, Time, y, df[,-(1:4)], m, seq(0, max(Time), length.out = 100), engine))
    list(alpha = out$ahat, beta = out$bhat, log.muZ = out$LamTau, lam0 = out$lamY,
         values = c(out$a.value, out$g.value))
}

##############################################################################
# Variance estimation 
##############################################################################
doREFit.Engine.Bootstrap <- function(DF, engine, stdErr) {
    res <- doREFit(DF, engine, NULL)
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
    betaMatrix <- matrix(0, B, p * 2)
    convergence <- rep(0, B)
    uID <- unique(DF$ID)
    for (i in 1:B) {
        sampled.id <- sample(unique(id), n, TRUE)
        ind <- unlist(sapply(sampled.id, function(x) which(id == x)))
        DF2 <- DF[ind,]
        DF2$id <- rep(1:n, clsz[sampled.id])
        betaMatrix[i,] <- with(doREFit(DF2, engine, NULL), c(alpha, beta))
        convergence[i] <- 1 * (betaMatrix[i,] %*% betaMatrix[i,] >
                               1e3 * with(res, c(alpha, beta) %*% c(alpha, beta)))
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

sdOut <- function(dat) {
    ol <- boxplot(dat, plot = FALSE)$out
    if (length(ol) > 0)
        dat <- dat[-which(dat %in% ol)]
    sd(dat, na.rm = TRUE)
}

doREFit.am.XCHWY.resampling <- function(DF, engine, stdErr) {
    res <- doREFit(DF, engine, NULL)
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
    aSE <- bSE <- da <- va <- db <- vb <- NA
    E <- matrix(rexp(n * B), nrow = n)
    Z <- matrix(rnorm(p * B), nrow = p)
    ua <- matrix(apply(Z, 2, function(x) alphaEq(res$alpha + n ^ (-0.5) * x, X, Y, T, cluster, mt)),
                 nrow = p)
    da <- t(apply(ua, 1, function(x) lm(n ^ (0.5) * x ~ t(Z))$coef[-1]))
    ua2 <- apply(E, 2, function(x) alphaEq(res$alpha, X, Y, T, cluster, mt, weights = x))
    va <- var(t(matrix(ua2, nrow = p)))
    if (qr(da)$rank == p)
        aVar <- solve(da) %*% va %*% t(solve(da))
    if (qr(da)$rank != p)
        aVar <- ginv(da) %*% va %*% t(ginv(da))
    aSE <- sqrt(diag(aVar))
    ub <- matrix(apply(Z, 2, function(x)
        betaEq(X, Y, T, cluster, mt, status[event == 0],
               res$zHat, res$alpha, res$beta + n ^ (-0.5) * x)), nrow = p)
    db <- t(apply(ub, 1, function(x) lm(n ^ (0.5) * x ~ t(Z))$coef[-1]))
    ub2 <- apply(E, 2, function(x) betaEq(X, Y, T, cluster, mt, status[event == 0],
                                          NULL, res$alpha, res$beta, weights = x))
    vb <- var(t(matrix(ub2, nrow = p)))
    if (qr(db)$rank == p)
        bVar <- solve(db) %*% vb %*% t(solve(db))
    if (qr(db)$rank != p)
        bVar <- ginv(db) %*% vb %*% t(ginv(db))
    bSE <- sqrt(diag(bVar))
    c(res, list(alphaSE = aSE, betaSE = bSE, da = da, va = va, db = db, vb = vb, B = stdErr@B))
}

doREFit.sc.XCYH.resampling <- function(DF, engine, stdErr) {
    id <- DF$id
    clsz <- unlist(lapply(split(id, id), length))
    ind <- cumsum(clsz)
    df <- DF
    m <- rep(clsz, clsz) - 1
    y <- rep(df$Time[ind], clsz)
    rmv <- sort(unique(rep(ind, clsz) * (m > 0)))[-1]
    df$Time <- df$Time * (m > 0)
    df <- df[-rmv, ]
    m <- as.numeric(m[-rmv])
    y <- as.numeric(y[-rmv])
    p <- ncol(df[,-(1:4)])
    if (is.na(match(engine@solver, c("dfsane", "BBsolve", "optim", "BBoptim")))) {
        print("Warning: Unidentified solver; BB::dfsane is used.")
        engine@solver <- "dfsane"
    }
    out <- with(df, sarmRV(id, Time, y, df[,-(1:4)], m, seq(0, max(Time), length.out = 100), engine))
    outSE <- with(df, sarmRV.sand(id, Time, y, df[,-(1:4)], m,
                                  a = out$ahat, b = out$ghat, Bootstrap = stdErr@B, engine))
    list(alpha = out$ahat, beta = out$bhat, alphaSE = outSE$alphaSE, betaSE = outSE$betaSE,
         gamma = out$ghat, gammaSE = sqrt(diag(outSE$ase[(p + 2):(2 * p + 1), (p + 2):(2 * p + 1)])),
         varMat = outSE$ase, log.muZ = out$LamTau, lam0 = out$lamY,
         values = c(out$a.value, out$g.value))
}

##############################################################################
# Nonparametric (~1)
##############################################################################

doNonpara.am.GL <- function(DF, alpha, beta, engine, stdErr) {
    DF0 <- subset(DF, event == 0)
    p <- ncol(DF0) - 4
    Y <- log(DF0$Time)
    X <- as.matrix(DF0[,-(1:4)])
    status <- DF0$status
    n <- nrow(DF0)
    d <- max(X %*% (alpha - beta), 0)
    tij <- log(DF$Time) - as.matrix(DF[,-(1:4)]) %*% alpha
    tij <- tij[DF$event == 1]
    yi <- Y - X %*% beta
    m <- aggregate(event ~ id, data = DF, sum)[,2]
    index <- c(1, cumsum(m)[-n] + 1)
    t0.rate <- unique(sort(tij)) # log scale
    t0.haz <- unique(sort(yi))
    rate <- .C("glRate", as.integer(n), as.integer(p), as.integer(index - 1), as.integer(m),
               as.integer(length(t0.rate)),
               as.double(yi - d), as.double(tij), as.double(X), as.double(t0.rate), 
               result = double(length(t0.rate)), 
               PACKAGE = "reReg")$result
    haz <- .C("glHaz", as.integer(n), as.integer(status), as.integer(length(t0.haz)),
              as.double(yi), as.double(t0.haz), result = double(length(t0.haz)), 
              PACKAGE = "reReg")$result
    rate0 <- approxfun(exp(t0.rate), rate, yleft = 0, yright = max(rate), method = "constant")
    haz0 <- approxfun(exp(t0.haz), haz, yleft = 0, yright = max(rate), method = "constant")
    list(rate0 = rate0, rate0.lower = NULL, rate0.upper = NULL, t0.rate = exp(t0.rate),
         haz0 = haz0, haz0.lower = NULL, haz0.upper = NULL, t0.haz = exp(t0.haz))
}

doNonpara.SE.am.GL <- function(DF, alpha, beta, engine, stdErr) {
    id <- subset(DF, event == 0)$id
    B <- stdErr@B
    PE <- doNonpara.am.GL(DF, alpha, beta, engine, NULL)
    rateMat <- matrix(NA, B, length(PE$t0.rate))
    hazMat <- matrix(NA, B, length(PE$t0.haz))
    for (i in 1:B) {
        sampled.id <- sample(id, length(id), TRUE)
        ind <- unlist(sapply(sampled.id, function(x) which(DF$id == x)))
        DF2 <- DF[ind,]
        DF2$id <- rep(1:length(id), table(DF$id)[sampled.id])
        tmp <- doNonpara.am.GL(DF2, alpha, beta, engine, NULL)
        rateMat[i,] <- tmp$rate0(PE$t0.rate)
        hazMat[i,] <- tmp$haz0(PE$t0.rate)
    }
    rlower <- apply(rateMat, 2, quantile, prob = .025)
    rupper <- apply(rateMat, 2, quantile, prob = .975)
    hlower <- apply(hazMat, 2, quantile, prob = .025)
    hupper <- apply(hazMat, 2, quantile, prob = .975)
    list(rate0 = PE$rate0,
         rate0.lower = approxfun(PE$t0.rate, rlower, yleft = 0, yright = max(rlower), method = "constant"),
         rate0.upper = approxfun(PE$t0.rate, rupper, yleft = 0, yright = max(rupper), method = "constant"),
         haz0 = PE$haz0,
         haz0.lower = approxfun(PE$t0.haz, hlower, yleft = 0, yright = max(hlower), method = "constant"),
         haz0.upper = approxfun(PE$t0.haz, hupper, yleft = 0, yright = max(hupper), method = "constant"))         
}

doNonpara.am.XCHWY <- function(DF, alpha, beta, engine, stdErr) {
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
    ## t0 <- seq(0, max(Y), length.out = 5 * nrow(DF))
    t0 <- sort(unique(T, Y))
    ng <- length(t0)
    ## Ya <- log(Y) + X %*% alpha
    ## Ta <- log(T) + X %*% alpha
    Ya <- Y * exp(X %*% alpha)
    Ta <- T * exp(X %*% alpha)
    lambda <- npMLE(Ya[which(event == 0)], Ta, Ya)
    ly <- npMLE(t0, Ta, Ya)
    zHat <- as.numeric(mt * max(ly) / lambda)
    ly <- ly / max(ly)
    win.ly <- max(ly)
    muZ <- mean(zHat)
    zHat <- ifelse(zHat %in% c("Inf", "NA", "NaN"), 0, zHat)
    Yb <- log(Y) + X %*% beta
    Yb <- Yb[which(cluster == 1)]
    hy <- sapply(t0, function(z) baseHaz(z, exp(Yb), zHat / muZ, status[event == 0]))
    win.hy <- max(hy)
    list(t0 = t0, lam = ly * muZ, lamU = rep(NA, ng), lamL = rep(NA, ng), 
         haz = hy, hazU = rep(NA, ng), hazL = rep(NA, ng))
}

doNonpara.cox.NA <- function(DF, alpha, beta, engine, stdErr) {
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
        tmp0 <- survfit(tmp, newdata = base)
        hy <- with(tmp0, approx(time, cumhaz, t0)$y)
        hyU <- -log(with(tmp0, approx(time, upper, t0)$y))
        hyL <- -log(with(tmp0, approx(time, lower, t0)$y))
    }
    if (all(X == 0)) {
        tmp <- coxph(Surv(T0, Time, event) ~ 1, data = DF)
        ly <- with(basehaz(tmp), approx(time, hazard, t0)$y)
        lyU <- -log(with(survfit(tmp), approx(time, upper, t0)$y))
        lyL <- -log(with(survfit(tmp), approx(time, lower, t0)$y))
        tmp <- coxph(Surv(Time, status) ~ 1, data = DF[event == 0,])
        hy <- with(basehaz(tmp), approx(time, hazard, t0)$y)
        hyU <- -log(with(survfit(tmp), approx(time, upper, t0)$y))
        hyL <- -log(with(survfit(tmp), approx(time, lower, t0)$y))
    }
    list(t0 = t0, lam = ly, lamU = lyU, lamL = lyL,
         haz = hy, hazU = hyU, hazL = hyL)
}

## doNonpara.SE.cox.NA <- function(DF, alpha, beta, engine, stdErr) {
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

doNonpara.cox.HW <- function(DF, alpha, beta, engine, stdErr) {
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
    ly <- ly / max(ly)
    win.ly <- max(ly)
    zHat <- ifelse(zHat %in% c("Inf", "NA", "NaN"), 0, zHat)
    muZ <- mean(zHat)
    Yb <- log(Y) ## + X %*% beta
    Yb <- Yb[event == 0]
    hy <- sapply(t0, function(z) baseHaz(z, exp(Yb),
                                         exp(as.matrix(X[event == 0,]) %*% beta) * zHat / muZ,
                                         status[event == 0]))
    win.hy <- max(hy)
    list(t0 = t0, lam = ly * muZ, lamU = rep(NA, ng), lamL = rep(NA, ng), 
         haz = hy, hazU = rep(NA, ng), hazL = rep(NA, ng))
}

doNonpara.SE.am.XCHWY <- function(DF, alpha, beta, engine, stdErr) {
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
    ## t0 <- seq(0, max(Y), length.out = 5 * nrow(DF))
    ng <- length(t0)
    ## Ya <- log(Y) + X %*% alpha
    ## Ta <- log(T) + X %*% alpha
    Ya <- Y * exp(X %*% alpha)
    Ta <- T * exp(X %*% alpha)
    lambda <- npMLE(Ya[event == 0], Ta, Ya)
    ly <- npMLE(t0, Ta, Ya)
    zHat <-  as.numeric(mt * max(ly) / lambda)
    ly <- ly / max(ly)
    E <- matrix(rexp(length(t0) * B), nrow = length(t0))
    lytmp <- apply(E, 2, function(x) npMLE(t0, Ta, Ya, x))
    lytmp <- apply(lytmp, 2, function(z) z / max(z))
    lyU <- apply(lytmp, 1, function(z) quantile(z, 0.975))
    lyL <- apply(lytmp, 1, function(z) quantile(z, 0.025))
    zHat <- ifelse(zHat %in% c("Inf", "NA", "NaN"), 0, zHat)
    Yb <- log(Y) + X %*% beta
    Yb <- Yb[event == 0]
    muZ <- mean(zHat)
    hy <- sapply(t0, function(z) baseHaz(z, exp(Yb), zHat / muZ, status[event == 0]))
    E <- matrix(rexp(n * B), nrow = n)
    hytmp <- apply(E, 2, function(z)
        sapply(t0, function(y)
            baseHaz(y, exp(Yb), zHat / muZ, status[event == 0], z)))
    hyU <- apply(hytmp, 1, function(z) quantile(z, 0.975))
    hyL <- apply(hytmp, 1, function(z) quantile(z, 0.025))
    list(t0 = t0, lam = ly * muZ, lamU = lyU * muZ, lamL = lyL * muZ, haz = hy, hazU = hyU, hazL = hyL)
}

doNonpara.SE.cox.HW <- function(DF, alpha, beta, engine, stdErr) {
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
    ## t0 <- seq(0, max(Y), length.out = 5 * nrow(DF))
    t0 <- sort(unique(T, Y))
    ng <- length(t0)
    Ya <- log(Y)
    Ta <- log(T)
    lambda <- npMLE(Ya[event == 0], Ta, Ya)
    ly <- npMLE(t0, exp(Ta), exp(Ya))
    zHat <- as.numeric(mt * max(ly) / (lambda * exp(as.matrix(X[event == 0,]) %*% alpha)))
    ly <- ly / max(ly)
    E <- matrix(rexp(length(t0) * B), nrow = length(t0))
    lytmp <- apply(E, 2, function(x) npMLE(t0, exp(Ta), exp(Ya), x))
    lytmp <- apply(lytmp, 2, function(z) z / max(z))
    lyU <- apply(lytmp, 1, function(z) quantile(z, 0.975))
    lyL <- apply(lytmp, 1, function(z) quantile(z, 0.025))
    zHat <- ifelse(zHat %in% c("Inf", "NA", "NaN"), 0, zHat)
    Yb <- log(Y) ## + X %*% beta
    Yb <- Yb[event == 0]
    muZ <- mean(zHat)
    hy <- sapply(t0, function(z) baseHaz(z, exp(Yb),
                                         exp(as.matrix(X[event == 0,]) %*% beta) * zHat / muZ,
                                         status[event == 0]))
    E <- matrix(rexp(n * B), nrow = n)
    hytmp <- apply(E, 2, function(z)
        sapply(t0, function(y)
            baseHaz(y, exp(Yb), exp(as.matrix(X[event == 0,]) %*% beta) * zHat / muZ,
                    status[event == 0], z)))
    hyU <- apply(hytmp, 1, function(z) quantile(z, 0.975))
    hyL <- apply(hytmp, 1, function(z) quantile(z, 0.025))
    list(t0 = t0, lam = ly * muZ, lamU = lyU * muZ, lamL = lyL * muZ,
         haz = hy, hazU = hyU, hazL = hyL)
}

##############################################################################
# Class Definition
##############################################################################

setClass("Engine",
         representation(tol = "numeric", a0 = "numeric", b0 = "numeric", solver = "character"),
         prototype(tol = 1e-7, a0 = 0, b0 = 0, solver = "dfsane"),
         contains="VIRTUAL")

setClass("cox.LWYY", contains = "Engine")
setClass("cox.HW", contains = "Engine")
setClass("am.XCHWY", contains = "Engine")
setClass("am.GL", contains = "Engine")
setClass("sc.XCYH",
         representation(eqType = "character"),
         prototype(eqType = "Logrank"), contains = "Engine")

setClass("stdErr")
setClass("bootstrap", representation(B = "numeric"),
         prototype(B = 100), contains="stdErr")
setClass("resampling", representation(B = "numeric"),
         prototype(B = 100), contains="stdErr")


##############################################################################
# Method Dispatch
##############################################################################
setGeneric("doREFit", function(DF, engine, stdErr) {standardGeneric("doREFit")})

setMethod("doREFit", signature(engine = "cox.LWYY", stdErr = "NULL"), doREFit.cox.LWYY)
setMethod("doREFit", signature(engine = "cox.HW", stdErr = "NULL"), doREFit.cox.HW)
setMethod("doREFit", signature(engine = "am.XCHWY", stdErr = "NULL"), doREFit.am.XCHWY)
setMethod("doREFit", signature(engine = "am.GL", stdErr = "NULL"), doREFit.am.GL)
setMethod("doREFit", signature(engine = "sc.XCYH", stdErr = "NULL"), doREFit.sc.XCYH)

setMethod("doREFit", signature(engine="Engine", stdErr="bootstrap"),
          doREFit.Engine.Bootstrap)

setMethod("doREFit", signature(engine="am.XCHWY", stdErr="resampling"),
          doREFit.am.XCHWY.resampling)
setMethod("doREFit", signature(engine="sc.XCYH", stdErr="resampling"),
          doREFit.sc.XCYH.resampling)

## --------------------------------------------------------------------------------------------------------
## Non-parametric 
## --------------------------------------------------------------------------------------------------------
setGeneric("doNonpara", function(DF, alpha, beta, engine, stdErr) {standardGeneric("doNonpara")})
setMethod("doNonpara", signature(engine = "cox.LWYY", stdErr = "NULL"), doNonpara.cox.NA)
setMethod("doNonpara", signature(engine = "cox.HW", stdErr = "NULL"), doNonpara.cox.HW)
setMethod("doNonpara", signature(engine = "am.XCHWY", stdErr = "NULL"), doNonpara.am.XCHWY)

setMethod("doNonpara", signature(engine = "cox.LWYY", stdErr = "bootstrap"), doNonpara.cox.NA)
setMethod("doNonpara", signature(engine = "cox.HW", stdErr = "bootstrap"), doNonpara.SE.cox.HW)
setMethod("doNonpara", signature(engine = "am.XCHWY", stdErr = "resampling"), doNonpara.SE.am.XCHWY)
setMethod("doNonpara", signature(engine = "am.XCHWY", stdErr = "bootstrap"), doNonpara.SE.am.XCHWY)

## GL method?
setMethod("doNonpara", signature(engine = "am.GL", stdErr = "bootstrap"), doNonpara.SE.am.GL)
setMethod("doNonpara", signature(engine = "am.GL", stdErr = "NULL"), doNonpara.am.GL)

## general model
setMethod("doNonpara", signature(engine = "sc.XCYH", stdErr = "bootstrap"), doNonpara.SE.am.XCHWY)
setMethod("doNonpara", signature(engine = "sc.XCYH", stdErr = "NULL"), doNonpara.SE.am.XCHWY)

##############################################################################
## User's Main Function
##############################################################################

#' Fits Semiparametric Regression Models for Recurrent Events and Failure Time
#'
#' Fits a survival model for the recurrent event data.
#' The rate function of the underlying process for the recurrent event process
#' can be specified as a Cox-type model, an accelerated mean model, an accelerated rate model, or a generalized scale-change model.
#' When a joint model is fitted (e.g., \code{method = "cox.HW"} or \code{method = "am.XCHWY"}),
#' the hazard function of the terminal event is either in a Cox model or an accelerated failure time model.
#'
#' The underlying models and assumptions are different for different methods.
#' The available methods are:
#' \describe{
#'   \item{\code{method == "cox.LWYY"}}{models the rate function  for the recurrent event through a Cox-type model.
#' This function returns results equivalent to that of \code{coxph}. See reference Lin et al. (2000).}
#'   \item{\code{method == "cox.HW"}}{jointly model the recurrent events and failure time.
#' This method assumes a Cox-type model for both the intensity function of the recurrent event process and the hazard function of the failure time.
#' Informative censoring is accounted for via a shared frailty variable.
#' See the references See reference Wang, Qin and Chiang(2001) and Huang and Wang (2004).}
#'   \item{\code{method == "am.GL"}}{jointly model the recurrent events and failure time.
#' This method assumes an accelerated mean model for both the recurrent event process and the failure time.
#' This method uses artificial censoring to allow for an unspecified association between the two types of outcomes.
#' Informative censoring is not allowed. See the reference Ghosh and Lin (2003).}
#'   \item{\code{method == "am.XCHWY"}}{jointly model the recurrent events and failure time.
#' This method assumes an accelerated mean model for both the recurrent event process and the failure time.
#' Informative censoring is accounted for via a shared frailty variable.
#' See the reference Xu et al. (2017).}
#'   \item{\code{method == "sc.XCYH"}}{models the rate function for the recurrent events via a generalized scale-change model that includes
#' Cox-type models, accelerated mean models, and accelerated rate models as special cases.
#' The methods also provide a hypothesis test of these submodels.
#' Informative censoring is accounted for via a shared frailty variable.}
#' }
#'
#' For available method for standard errors estimation are:
#'
#' \describe{
#'   \item{NULL}{means do not calculate standard errors.}
#'   \item{\code{se == "resampling"}}{performs the efficient resampling-based sandwich estimator that works with \code{method == "am.XCHWY"} and \code{method == "sc.XCYH"}.}
#'   \item{\code{"bootstrap"}}{works with all fitting methods.}
#' }
#' 
#' @param formula a formula object, with the response on the left of a "~" operator, and the terms on the right.
#' The response must be a recurrent event survival object as returned by function \code{reSurv}.
#' @param data  an optional data frame in which to interpret the variables occurring in the "formula".
#' @param B a numeric value specifies the number of resampling for variance estimation.
#' When B = 0, only the point estimates will be displayed.
#' @param method a character string specifying the underlying model. See \code{Details}.
#' @param se a character string specifying the method for standard error estimation. See \code{Details}.
#' @param plot.ci a logical value indicating whether the 95% confidence interval for the estimated cumulative rate function
#' and/or the estimated cumulative hazard function should be computed. Default is "FALSE".
#' @param contrasts an optional list.
#' @param control a list of control parameters.
#'
#' @export
#' @references Xu, G., Chiou, S.H., Huang, C.-Y., Wang, M.-C. and Yan, J. (2017). Joint Scale-change Models for Recurrent Events and Failure Time.
#' \emph{Journal of the American Statistical Association} 112.518: 796-805.
#' @references Lin, D., Wei, L., Yang, I. and Ying, Z. (2000). Semiparametric Regression for the Mean and Rate Functions of Recurrent Events.
#' \emph{Journal of the Royal Statistical Society: Series B (Methodological)}, \bold{62}: 711 -- 730.
#' @references Wang, M.C., Qin, J., and Chiang, C.T. (2001). Analyzing Recurrent Event Data with Informative Censoring.
#' \emph{Journal of the American Statistical Association} \bold{96}{455}: 1057--1065.
#' @references Ghosh, D. and D.Y. Lin (2003). Semiparametric Analysis of Recurrent Events Data in the Presence of Dependent Censoring.
#' \emph{Biometrics}, \bold{59}: 877 -- 885.
#' @references Huang, C.Y. and Wang, M.C. (2004). Joint Modeling and Estimation for Recurrent Event Processes and Failure Time Data.
#' \emph{Journal of the American Statistical Association} \bold{99}(468), 1153--1165.
#'
#' @seealso \code{\link{reSurv}}
#'
#' @examples
#' ## readmission data
#' data(readmission, package = "frailtypack")
#' set.seed(123)
#' ## Acceralted Mean Model
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
#'
#' \dontrun{
#' ## simulation data
#' simDat <- function(n, a, b, latent = FALSE) {
#'     ## setting rate function
#'     Lam.f <- function(t, z, x, w) .5 * z * exp(-x + w) * log(1 + t * exp(x))
#'     Lam.f0 <- function(t) .5 * log(1 + t)
#'     invLam.f  <- function(t, z, x, w) (exp((2 / z) * exp(x - w) * t )- 1) / exp(x)
#'     ## setting hazard funciton
#'     ## Haz.f0 <- function(t) .5 * log(1 + t) # assume constant hazard for now
#'     dat <- NULL
#'     for (id in 1:n) {
#'         z <- ifelse(latent, rgamma(1, 4, 4), 1)
#'         x1 <- rnorm(1)
#'         x2 <- rnorm(1)
#'         x <- c(x1, x2)
#'         cen <- rexp(1, z * exp(x %*% b) / 60) ## this gives constant hazard of 1/60
#'         y <- min(cen, 60)
#'         D <- 1 * (cen == y)
#'         tmpt <- NULL
#'         while(sum(tmpt) < Lam.f(y, z, c(x %*% a), c(x %*% b))) {
#'             tmpt <- c(tmpt, rexp(1))
#'         }
#'         m <- length(tmpt) - 1
#'         if (m > 0) {
#'             tt <- invLam.f(cumsum(tmpt[1:m]), z, c(x %*% a), c(x %*% b))
#'             dat <- rbind(dat, cbind(ID = id, Time = c(tt[order(tt)], y),
#'                                     event = c(rep(1, m), 0), status = c(rep(0, m), D),
#'                                     Z = z, M = m, X1 = x1, X2 = x2))
#'         } else {
#'             dat <- rbind(dat, cbind(ID = id, Time = y, event = 0, status = D,
#'                                     Z = z, M = m, X1 = x1, X2 = x2))
#'         }
#'     }
#'     return(data.frame(dat))
#' }
#' set.seed(2017)
#' dat <- simDat(200, c(0, 0), c(0, 0), TRUE) ## generate data under informative censoring
#' fm <- reSurv(Time, ID, event, status) ~ X1 + X2
#' fit.HW <- reReg(fm, data = dat, method = "cox.HW", B = 5)
#' }
reReg <- function(formula, data, B = 200, 
                  method = c("cox.LWYY", "cox.HW", "am.GL", "am.XCHWY", "sc.XCYH"),
                  se = c("NULL", "bootstrap", "resampling"), plot.ci = FALSE,
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
    engine <- do.call("new", c(list(Class=method), engine.control))
    if (se == "NULL")
        stdErr <- NULL
    else {
        stdErr.control <- control[names(control) %in% names(attr(getClass(se), "slots"))]
        stdErr <- do.call("new", c(list(Class=se), stdErr.control))
        stdErr@B <- B
    }
    p <- ncol(DF) - 4
    if (length(engine@a0) == 1) engine@a0 <- rep(engine@a0, p)
    if (length(engine@b0) == 1) engine@b0 <- rep(engine@b0, p)
    if (method %in% c("cox.HW", "sc.XCYH")) {
        if (length(engine@b0) == 1) engine@b0 <- rep(engine@b0, p + 1)
        if (length(engine@b0) == p) engine@b0 <- c(0, engine@b0)
    }
    if (formula == ~1) {
        fit <- NULL
        fit$alpha <- fit$beta <- rep(NA, p)
        fit$muZ <- NA
        if (plot.ci) {
            stdErr.np.control <- control[names(control) %in% names(attr(getClass("bootstrap"), "slots"))]
            stdErr.np <- do.call("new", c(list(Class = "bootstrap"), stdErr.np.control))
            stdErr.np@B <- B
            fit <- c(fit, doNonpara(DF = DF, alpha = fit$alpha, beta = fit$beta,
                                    engine = engine, stdErr = stdErr.np))
        } else {
            fit <- c(fit, doNonpara(DF = DF, alpha = 0, beta = 0,
                                    engine = engine, stdErr = stdErr))
        }
    } else {
        fit <- doREFit(DF = DF, engine = engine, stdErr = stdErr)
        if (method != "sc.XCYH") {
            if (plot.ci) {
                stdErr.np.control <- control[names(control) %in% names(attr(getClass("bootstrap"), "slots"))]
                stdErr.np <- do.call("new", c(list(Class = "bootstrap"), stdErr.np.control))
                stdErr.np@B <- B
            fit <- c(fit, doNonpara(DF = DF, alpha = fit$alpha, beta = fit$beta,
                                    engine = engine, stdErr = stdErr.np))
            } else {
                fit <- c(fit, doNonpara(DF = DF, alpha = fit$alpha, beta = fit$beta, 
                                        engine = engine, stdErr = NULL))
            }
        }
            
    }
    class(fit) <- "reReg"
    fit$DF <- DF
    fit$call <- Call
    fit$varNames <- names(DF)[-(1:4)]
    fit$method <- method
    fit$se <- se
    fit
}

##############################################################################
##############################################################################

npMLE <- function(t, tij, yi, weights = NULL) {
    if (is.null(weights))
        weights <- rep(1, length(yi))
    ttmp <- tij[tij != yi]
    ord <- order(ttmp)
    sl <- unique(ttmp[ord])
    l <- ifelse(min(t) < max(sl), which(sl > min(t))[1], length(sl))
    ## res <- vector("double", 1)
    ## tmp <- sl[l:length(sl)]
    tmp <- sl[l:length(sl)]
    tmp <- rev(tmp)
    tij <- rev(tij)
    yi <- rev(yi)
    ## yi <- ifelse(is.infinite(yi), max(yi[!is.infinite(yi)]), yi)
    ## tij <- ifelse(is.infinite(tij), max(tij[!is.infinite(tij)]), tij)
    res <- vector("double", length(tmp)) + 1
    res <- .C("plLambda", as.double(tmp), as.double(tij), as.double(yi), as.double(weights), 
              as.integer(length(tmp)), as.integer(length(yi)),
              out = as.double(res), PACKAGE = "reReg")$out
    out <- rev(res)[sapply(t, function(x) which(rev(tmp) >= x)[1])]
    out <- ifelse(is.na(out), 0, out)
    out <- exp(-out)
}


baseHaz <- function(t, Y, zHat, delta, weights  = NULL) {
    if (is.null(weights)) 
        weights <- rep(1, length(Y))
    ind <- which(delta == 1 & Y <= t)
    temp2 <- tmp <- weights[order(Y)]
    ## temp2 <- c(tmp[1], diff(cumsum(tmp)))
    ## temp2[order(Y)] <- temp2
    temp2[order(Y)] <- tmp
    if (length(ind) > 0) {
        out <- sapply(ind, function(x) temp2[x] / sum(zHat * weights * (Y >= Y[x])))
    }
    if (length(ind) == 0)
        out <- 0
    sum(out)
}

alphaEq <- function(alpha, X, Y, T, cluster, mt, weights = NULL) {
    n <- sum(cluster == 1)
    if (is.null(weights))
        weights <- rep(1, n)
    p <- ncol(X)
    Ystar <- Y * exp(X %*% alpha)
    Tstar <- T * exp(X %*% alpha)
    ## Ystar <- log(Y) + X %*% alpha
    ## Tstar <- log(T) + X %*% alpha
    Lambda <- npMLE(Ystar[which(cluster == 1)], Tstar, Ystar,
                    weights = rep(weights, diff(c(which(cluster ==1), length(cluster)+1))))
    res <- vector("double", p * length(weights) %/% n)
    res <- .C("alphaEqC", as.double(X[which(cluster == 1), ]), as.double(Lambda),
              as.integer(mt), as.integer(n), as.integer(p),
              out = as.double(res), PACKAGE = "reReg")$out
    res / rep(n * unlist(lapply(split(weights, rep(1:(length(weights) %/% n), each = n)), sum)), each = p)
}

betaEq <- function(X, Y, T, cluster, mt, delta, zHat = NULL, alpha, beta, weights = NULL) {
    p <- ncol(X)
    n <- sum(cluster == 1)
    if (is.null(weights))
        weights <- rep(1, n)
    if (is.null(zHat)) {
        Ystar <- Y * exp(X %*% alpha)
        Tstar <- T * exp(X %*% alpha)
        ## Ystar <- log(Y) + X %*% alpha
        ## Tstar <- log(T) + X %*% alpha
        lambda <- npMLE(Ystar[which(cluster == 1)], Tstar, Ystar,
                        weights = rep(weights, diff(c(which(cluster ==1), length(cluster)+1))))
        zHat <- as.numeric(weights * mt / lambda)
        zHat <- ifelse(zHat %in% c("Inf", "NA", "NaN"), 0, zHat)
    }
    Y <- log(Y) + X %*% beta
    Y <- Y[which(cluster == 1)]
    X <- X[which(cluster == 1), ]
    ## delta <- delta[which(cluster == 1)]
    res <- vector("double", p * length(weights) %/% n)
    res <- .C("betaEst", as.double(Y), as.double(X), as.double(delta), as.double(zHat),
              as.double(weights), as.integer(n), as.integer(p),
              as.integer(length(weights) %/% n), 
              out = as.double(res), PACKAGE = "reReg")$out
    res / n
}

## ghoshU2 <- function(alpha, beta, T, Y, X, cl) {
##     ## dim(X) = n by p, dim(Y) = n by 1, dim(T) > n by 1
##     d <- max(X %*% (alpha - beta), 0)
##     TT <- log(T) - rep(X %*% alpha, cl)
##     TY <- log(Y) - X %*% beta - d
##     p <- ncol(X)
##     .C("ghosh", as.double(TT), as.double(TY), as.double(X), as.integer(cl),
##        as.integer(c(0, cumsum(cl)[-length(cl)])),
##        as.integer(nrow(X)), as.integer(p), 
##        out = as.double(double(p)), PACKAGE = "reReg")$out
## }

LWYYeq <- function(beta, X, Y, T, cl) {
    p <- ncol(X)
    res <- vector("double", p)
    wgt <- exp(X %*% beta)
    .C("lwyy", as.double(T), as.double(Y), as.double(X), as.double(wgt), as.integer(cl),
       as.integer(c(0, cumsum(cl)[-length(cl)])), as.integer(nrow(X)), as.integer(p),        
       out = as.double(double(p)), PACKAGE = "reReg")$out       
}


##########################################################################################
## Paper 2: More general models
##########################################################################################

sarm <- function(X, Y, T, id, cluster, method, B = 200) {
    n <- sum(cluster == 1)
    mt <- unlist(lapply(split(cluster, id), length)) - 1
    p <- ncol(X)
    alpha <- beta <- gamma <- rep(0, p)
    muZ <- NULL
    if (method == "M1") {
        gamma <- c(0, gamma)
        out <- BBsolve(gamma, coefEq, alpha = alpha, X = X, Y = Y, T = T,
                       cluster = cluster, mt = mt, weights = NULL,
                       quiet = TRUE, control = list(M = c(1, 10)))
        muZ <- out$par[1]
        alpha <- out$par[-1]
    }
    if (method == "M3") {
        alpha <- BBsolve(alpha, M1eq, X = X, Y = Y, T = T, cluster = cluster, weights = NULL,
                         quiet = TRUE, control = list(M = c(1, 10)))$par
        gamma <- c(0, gamma)
        if (alpha %*% alpha > 100) {
            beta <- c(0, alpha)
        }
        else {
            out <- BBsolve(gamma, coefEq, alpha = alpha, X = X, Y = Y, T = T,
                           cluster = cluster, mt = mt, weights = NULL,
                           quiet = TRUE, control = list(M = c(1, 10)))
            muZ <- out$par[1]
            beta <- out$par[-1] + alpha
        }
    }
    if (method == "M2") {
        ## gamma <- c(0, gamma)
        ## out <- BBsolve(gamma, coefEq, alpha = NULL, X = X, Y = Y, T = T,
        ##                cluster = cluster, mt = mt, weights = NULL,
        ##                quiet = TRUE, control = list(M = c(1, 10)))
        ## muZ <- out$par[1]
        ## alpha <- out$par[-1]
        ## out <- dfsane(alpha, M1eq, X = X, Y = Y, T = T,
        ##                 cluster = cluster, weights = NULL,
        ##                 control = list(NM = TRUE, M = 1, noimp = 50, trace = FALSE))
        out <- BBsolve(alpha, M1eq, X = X, Y = Y, T = T, cluster = cluster, weights = NULL,
                       quiet = TRUE, control = list(M = c(1, 10)))
        alpha <- out$par
    }
    list(alpha = alpha, beta = beta, muZ = muZ)
}

HWeq <-function(gamma, X, Y, T, cluster, mt) {
    n <- sum(cluster == 1)
    Lambda <- npMLE(Y[cluster == 1], T, Y)
    res <- vector("double", length(gamma))
    p <- ncol(X)
    .C("sarm1", as.double(X), as.double(Lambda), as.double(rep(1, n)),
       as.double(gamma), as.integer(mt), as.integer(n), as.integer(p), as.integer(1),
       out = as.double(rep(0, p)), PACKAGE = "reReg")$out                            
}


HWeq2 <-function(beta, X, Y, delta, zHat) {
    n <- nrow(X)
    p <- ncol(X)
    res <- vector("double", p)
    res <- .C("HWb", as.double(Y), as.double(X), as.double(delta), as.double(zHat),
              as.double(X %*% beta), as.integer(n), as.integer(p), as.integer(1),
              out = as.double(res), PACKAGE = "reReg")$out
    res / n
}
    
coefEq <- function(alpha, gamma, X, Y, T, cluster, mt, weights = NULL) {
    n <- sum(cluster == 1)
    if (is.null(weights))
        weights <- rep(1, n)
    if (is.null(alpha))
        alpha <- -1 * gamma[-1]
    Ytmp <- log(Y) + X %*% alpha
    Ttmp <- log(T) + X %*% alpha    
    Lambda <- npMLE(Ytmp[cluster == 1], Ttmp, Ytmp,
                    weights = rep(weights, diff(c(which(cluster ==1), length(cluster)+1))))
    Lambda <- Lambda / max(Lambda)
    X <- X[cluster == 1,]
    p <- ncol(X)
    res <- vector("double", (p + 1) * length(weights) %/% n)
    res <- .C("sarm1", as.double(cbind(1, X)), as.double(Lambda), as.double(weights),
              as.double(gamma), as.integer(mt), as.integer(n), as.integer(p+1),
              as.integer(length(weights) %/% n),
              out = as.double(res), PACKAGE = "reReg")$out
    res / n    
}


## Martingal approach
M1eq <- function(alpha, X, Y, T, cluster, weights = NULL) {
    n <- sum(cluster == 1)
    ## mi <- unlist(lapply(split(cluster, id), length)) - 1
    if (is.null(weights))
        weights <- rep(1, n)
    p <- ncol(X)
    Ytmp <- log(Y) + X %*% alpha
    Ttmp <- log(T) + X %*% alpha
    ## Ytmp <- Y * exp(X %*% alpha)
    ## Ttmp <- T * exp(X %*% alpha)
    ind <- which(T != Y)
    Ttmp <- Ttmp[ind]
    Ytmp <- Ytmp[ind]
    X <- X[ind,]
    res <- vector("double", p * length(weights) %/% n)
    res <- .C("sarm2", as.double(X), as.double(Ttmp), as.double(Ytmp), as.double(weights), 
              as.integer(length(Ttmp)), as.integer(p), as.integer(length(weights) %/% n),
              out = as.double(res), PACKAGE = "reReg")$out
    res
}

## Method 1: Z\lambda(t)e^xa  CY's 2004 JASA
## Method 2: Z\lambda(te^xa)
## Method 3: Z\lambda(te^xa)e^xb
## Method 4: Z\lambda(te^xa)e^xa

#########################################################
## General modes in R, need to move this to C sometimes
#########################################################

sarmRV <- function(id, Tij, Yi, X, M, lamEva = NULL, engine) {
    n <- length(unique(id))
    X <- as.matrix(X)
    p <- ncol(X)
    mm <- matrix((M > 0), nrow(X), ncol(X))
    ## U1RV <- function(a) {
    ##     tx <- as.vector(log(Tij) + X %*% a)
    ##     yx <- as.vector(log(Yi) + X %*% a)
    ##     tx2 <-  outer(tx, tx,">=")
    ##     txy <-  outer(tx, yx, "<=")
    ##     A <- (tx2 * txy) %*% (X * mm)
    ##     B <- (tx2 * txy) %*% mm
    ##     B[B == 0] <- 1e10
    ##     colSums((X - A / B) * mm) 
    ## }
    index <- c(1, cumsum(tabulate(id))[-n] + 1)
    U1RV <- function(a) {
        tx <- as.vector(log(Tij) + X %*% a)
        yx <- as.vector(log(Yi) + X %*% a)
        tx <- ifelse(tx == -Inf, -1e10, tx)
        yx <- ifelse(yx == -Inf, -1e10, yx)
        if (engine@eqType %in% c("Logrank", "logrank")) 
            return(.C("scaleChangeLog", as.integer(n), as.integer(p), as.integer(index - 1),
                      as.integer(M[index]), as.double(yx), as.double(tx), as.double(X[index,]),
                      as.double(rep(1, length(index))), 
                      result = double(p), PACKAGE = "reReg")$result)
        if (engine@eqType %in% c("Gehan", "gehan"))
            return(.C("scaleChangeGehan", as.integer(n), as.integer(p), as.integer(index - 1),
                      as.integer(M[index]), as.double(yx), as.double(tx), as.double(X[index,]),
                      as.double(rep(1, length(index))), 
                      result = double(p), PACKAGE = "reReg")$result / n)
    }
    if (engine@solver %in% c("dfsane", "BBsolve")) {
        suppressWarnings(
            fit.a <- do.call(engine@solver, list(par = engine@a0, fn = U1RV,
                                                 quiet = TRUE, control = list(trace = FALSE))))
        a.value <- fit.a$residual
    }
    if (engine@solver == "BBoptim") {
        suppressWarnings(fit.a <- do.call(engine@solver, list(par = engine@a0, fn = function(x) sum(U1RV(x)^2), quiet = TRUE, control = list(trace = FALSE))))
        a.value <- fit.a$value
    }
    if (engine@solver == "optim") {
        suppressWarnings(fit.a <- do.call(engine@solver, list(par = engine@a0, fn = function(x) sum(U1RV(x)^2), control = list(trace = FALSE))))
        a.value <- fit.a$value
    }
    ahat <- fit.a$par
    ## tx <- as.vector(Tij * exp(X %*% ahat))
    ## yx <- as.vector(Yi * exp(X %*% ahat))
    tx <- as.vector(log(Tij) + X %*% ahat)
    yx <- as.vector(log(Yi) + X %*% ahat)
    tx <- ifelse(tx == -Inf, 0, tx)
    yx <- ifelse(yx == -Inf, 0, yx)
    tx2 <- outer(tx, tx,">=")
    txy <- outer(tx, yx, "<=")
    vv <- matrix((M > 0), nrow(X), n)
    yx0 <- as.numeric(unlist(lapply(split(yx, id), unique)))
    txy0 <- outer(tx, yx0, ">=")
    Rn <- (tx2 * txy) %*% (M > 0)
    Rn[Rn == 0] <- 1e10
    ## Lam <- apply(1 - (txy0 * vv) / matrix(Rn, nrow(X), n), 2, prod)
    Lam <- exp(-colSums((txy0 * vv) / matrix(Rn, nrow(X), n))) ## gives \Lambda(Y_i)
    Lam[Lam == 0] <- 1e10
    ind <- cumsum(unlist(lapply(split(id, id), length)))
    U2RV <- function(b) {
        tmp <- ifelse(M[ind] / Lam > 1e5, (M[ind] + .01) / (Lam + .01), M[ind] / Lam)
        as.numeric(t(cbind(1, X[ind,])) %*% (tmp - exp(cbind(1, X[ind,]) %*% b))) / n
        ## as.numeric(t(cbind(1, X[ind,])) %*% (M[ind] / Lam - exp(cbind(1, X[ind,]) %*% b))) / n
    }
    if (engine@solver %in% c("dfsane", "BBsolve")) {
        suppressWarnings(fit.g <- do.call(engine@solver, list(par = engine@b0, fn = U2RV, quiet = TRUE, control = list(trace = FALSE))))
        g.value <- fit.g$residual
    }
    if (engine@solver == "BBoptim") {
        suppressWarnings(fit.g <- do.call(engine@solver, list(par = engine@b0, fn = function(x) sum(U2RV(x)^2), quiet = TRUE, control = list(trace = FALSE))))
        g.value <- fit.g$value
    }
    if (engine@solver == "optim") {
        suppressWarnings(fit.g <- do.call(engine@solver, list(par = engine@b0, fn = function(x) sum(U2RV(x)^2), control = list(trace = FALSE))))
        g.value <- fit.g$value
    }
    ghat <- fit.g$par
    y0 <- exp(yx[ind])
    list(ahat = ahat, bhat = ghat[-1] + ahat, ghat = ghat, LamTau = ghat[1],
         lamY = stepfun(y0[order(y0)], c(0, Lam[order(y0)])),
         a.convergence = fit.a$convergence, a.value = a.value, g.convergence = fit.g$convergence, g.value = g.value)
}

varOut <- function(dat, na.rm = TRUE) {
    dat[which(dat %in% boxplot(dat, plot = FALSE)$out)] <- NA
    dat <- dat[complete.cases(dat),]
    var(dat, na.rm = na.rm)
}


sarmRV.sand <- function(id, Tij, Yi, X, M, a = NULL, b = NULL, Bootstrap = 200, engine) {
  n <- length(unique(id))
  X <- as.matrix(X)
  p <- ncol(X)
  mm <- matrix((M > 0), nrow(X), ncol(X))
  clut <- as.numeric(unlist(lapply(split(id, id), length)))
  tmpE <- matrix(rexp(n * Bootstrap), ncol = n)
  tmpN <- matrix(rnorm((2 * p + 1) * Bootstrap), ncol = 2 * p + 1)
  ## Sn <- function(a, b, e) {
  ##     ## Part S1
  ##     e1 <- rep(e, clut)
  ##     tx <- as.vector(Tij * exp(X %*% a))
  ##     yx <- as.vector(Yi * exp(X %*% a))
  ##     ee1 <- matrix(e1, nrow(X), ncol(X))
  ##     tx2 <-  outer(tx, tx,">=")
  ##     txy <-  outer(tx, yx, "<=")
  ##     A <- (tx2 * txy) %*% (X * mm * ee1)
  ##     B <- (tx2 * txy) %*% (mm * ee1)
  ##     B[B == 0] <- 1e15
  ##     if (engine@eqType %in% c("Logrank", "logrank"))
  ##         s1 <- colSums((X - A / B) * (mm * ee1)) / n
  ##     if (engine@eqType %in% c("Gehan", "gehan"))
  ##         s1 <- colSums((X * B- A) * (mm * ee1)) / n^2
  ##     ## lambda
  ##     vv <- matrix((M > 0), nrow(X), n)
  ##     yx0 <- as.numeric(unlist(lapply(split(yx, id), unique)))
  ##     txy0 <- outer(tx, yx0, ">=")
  ##     Rn <- (tx2 * txy) %*% (e1 * (M > 0))
  ##     Rn[Rn == 0] <- 1e15
  ##     ## Lam <- apply(1 - (txy0 * vv * e1) / matrix(Rn, nrow(X), n), 2, prod)
  ##     Lam <- exp(-colSums((txy0 * vv * e1) / matrix(Rn, nrow(X), n)))
  ##     Lam[Lam == 0] <- 1e15
  ##     ind <- cumsum(unlist(lapply(split(id, id), length)))
  ##     ee2 <- matrix(e, nrow(X[ind,]), ncol(X) + 1)
  ##     s2 <- as.numeric(t(cbind(1, X[ind,]) * ee2) %*%
  ##                      (M[ind] / Lam - exp(cbind(1, X[ind,]) %*% b))) / n
  ##     ## list(s1 = s1, s2 = s2)
  ##     return(c(s1, s2))
  ## }
  index <- c(1, cumsum(tabulate(id))[-n] + 1)
  Sn <- function(a, b, e) {
      ## Part S1
      e1 <- rep(e, clut)
      tx <- as.vector(log(Tij) + X %*% a)
      yx <- as.vector(log(Yi) + X %*% a)
      tx <- ifelse(tx == -Inf, -1e10, tx)
      yx <- ifelse(yx == -Inf, -1e10, yx)
      ## tx <- as.vector(Tij * exp(X %*% a))
      ## yx <- as.vector(Yi * exp(X %*% a))
      ## ee1 <- matrix(e1, nrow(X), ncol(X))
      tx2 <-  outer(tx, tx,">=")
      txy <-  outer(tx, yx, "<=")
      ## A <- (tx2 * txy) %*% (X * mm * ee1)
      ## B <- (tx2 * txy) %*% (mm * ee1)
      ## B[B == 0] <- 1e15
      if (engine@eqType %in% c("Logrank", "logrank")) {
          s1 <- .C("scaleChangeLog", as.integer(n), as.integer(p), as.integer(index - 1),
                   as.integer(M[index]), as.double(yx), as.double(tx), as.double(X[index,]), as.double(e),
                   result = double(p), PACKAGE = "reReg")$result / n}
      if (engine@eqType %in% c("Gehan", "gehan")) {
          s1 <- .C("scaleChangeGehan", as.integer(n), as.integer(p), as.integer(index - 1),
                   as.integer(M[index]), as.double(yx), as.double(tx), as.double(X[index,]), as.double(e),
                   result = double(p), PACKAGE = "reReg")$result / n^2}
      vv <- matrix((M > 0), nrow(X), n)
      yx0 <- as.numeric(unlist(lapply(split(yx, id), unique)))
      txy0 <- outer(tx, yx0, ">=")
      Rn <- (tx2 * txy) %*% (e1 * (M > 0))
      Rn[Rn == 0] <- 1e15
      ## Lam <- apply(1 - (txy0 * vv * e1) / matrix(Rn, nrow(X), n), 2, prod)
      Lam <- exp(-colSums((txy0 * vv * e1) / matrix(Rn, nrow(X), n)))
      Lam[Lam == 0] <- 1e15
      ind <- cumsum(unlist(lapply(split(id, id), length)))
      ee2 <- matrix(e, nrow(X[ind,]), ncol(X) + 1)
      s2 <- as.numeric(t(cbind(1, X[ind,]) * ee2) %*%
                       (M[ind] / Lam - exp(cbind(1, X[ind,]) %*% b))) / n
      ## list(s1 = s1, s2 = s2)
      return(c(s1, s2))
  }
  ## V <- varOut(t(apply(tmpE, 1, function(x) Sn(a, b, x)))) ## / sqrt(n)
  V <- var(t(apply(tmpE, 1, function(x) Sn(a, b, x)))) ## / sqrt(n)
  tmp <- t(apply(tmpN, 1, function(x)
                 sqrt(n) * Sn(a + x[1:p] / sqrt(n), b + x[-(1:p)] / sqrt(n), rep(1, n))))
  ## J0 <- t(coef(lm(tmp[,1:p] ~ tmpN[,1:p] - 1)))
  ## Jtmp <- t(coef(lm(tmp[,-c(1:p)] ~ tmpN - 1)))
  J0 <- t(coef(lm(tmp[,1:p] ~ tmpN[,1:p]))[-1,])
  Jtmp <- t(coef(lm(tmp[,-c(1:p)] ~ tmpN))[-1,])
  J <- rbind(cbind(J0, matrix(0, ncol = p + 1, nrow = nrow(J0))), cbind(Jtmp))
  if (qr(J)$rank == nrow(J)) J <- solve(J) else J <- ginv(J)
  if (qr(J0)$rank == nrow(J0)) J0 <- solve(J0) else J0 <- ginv(J0)
  ase <- J %*% V %*% t(J)
  list(ase = ase, J = J, V = V,
       alphaSE = sqrt(diag(J0 %*% V[1:p, 1:p] %*% t(J0))),## sqrt(diag(ase)[1:p]),
       betaSE = sqrt(diag(ase[1:p, 1:p] + ase[(p + 2):(2 * p + 1), (p + 2):(2 * p + 1)] +
         2 * ase[1:p, (p + 2):(2 * p + 1)])))
}
