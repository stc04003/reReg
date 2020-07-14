npFit.am.GL <- function(DF, alpha, beta, engine, stdErr) {
    alpha <- -alpha
    beta <- -beta
    DF0 <- subset(DF, event == 0)
    X <- as.matrix(DF0[,-(1:6)])
    p <- ncol(X)
    status <- DF0$terminal
    n <- nrow(DF0)
    ## d <- max(X %*% (alpha - beta), 0)
    d <- max(X %*% (alpha - beta))
    tij <- log(DF$time2) - as.matrix(DF[,-(1:6)]) %*% alpha
    tij <- tij[DF$event == 1]
    yi <- log(DF0$time2) - X %*% beta ## beta or alpha?
    m <- aggregate(event ~ id, data = DF, sum)[,2]
    index <- c(1, cumsum(m)[-n] + 1)
    t0.rate <- sort(unique(tij)) # log scale
    t0.haz <- sort(unique(yi))
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
    list(Lam0 = rate0, Lam0.lower = NULL, Lam0.upper = NULL,
         Haz0 = haz0, Haz0.lower = NULL, Haz0.upper = NULL)
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
    list(Lam0 = PE$rate0,
         Lam0.lower = approxfun(PE$t0.rate, rl, yleft = 0, yright = max(rl), method = "constant"),
         Lam0.upper = approxfun(PE$t0.rate, ru, yleft = 0, yright = max(ru), method = "constant"),
         Haz0 = PE$haz0,
         Haz0.lower = approxfun(PE$t0.haz, hl, yleft = 0, yright = max(hl), method = "constant"),
         Haz0.upper = approxfun(PE$t0.haz, hu, yleft = 0, yright = max(hu), method = "constant"))
}

npFit.SE.cox.GL <- function(DF, alpha, beta, engine, stdErr) {
    id <- DF$id
    event <- DF$event
    X <- as.matrix(DF[,-c(1:6)])
    p <- ncol(X)
    T <- DF$time2
    mt <- aggregate(event ~ id, data = DF, sum)$event
    Y <- rep(DF$time2[!event], mt + 1)
    t0.rate <- sort(unique(T))
    xb <- exp(X[!event,] %*% beta) 
    cl <- mt + 1
    rate <- .C("glCoxRate", as.double(ifelse(T == Y, 1e5, T)), as.double(Y[!event]),
               as.double(xb), as.double(engine@wgt), as.double(t0.rate), as.integer(length(t0.rate)), 
               as.integer(length(T)), as.integer(cl), as.integer(c(0, cumsum(cl)[-length(cl)])),
               as.integer(sum(!event)), out = double(length(t0.rate)), PACKAGE = "reReg")$out
    rate0 <- approxfun(t0.rate, rate, yleft = 0, yright = max(rate), method = "constant")
    fit.coxph <- coxph(Surv(time2, terminal) ~ .,
                       data = subset(DF, !event, select = -c(id, event, time1, origin)))
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
        fit.coxphB <- coxph(Surv(time2, terminal) ~ .,
                            data = subset(DF2, !event, select = -c(id, event, time1, origin)))
        cumHaz <- basehaz(fit.coxphB)
        ## cumHaz$hazard <- cumHaz$hazard / max(cumHaz$hazard)
        X0 <- as.matrix(DF2[!DF2$event, -(1:6)])
        wgt <- sapply(exp(X0 %*% coef(fit.coxphB)), function(x)
            with(cumHaz, approxfun(time, exp(-hazard * x), yleft = 1, yright = min(exp(-hazard * x)),
                                   method = "constant"))(DF2$time2))
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
    list(Lam0 = rate0,
         Lam0.lower = approxfun(t0.rate, rl, yleft = 0, yright = max(rl), method = "constant"),
         Lam0.upper = approxfun(t0.rate, ru, yleft = 0, yright = max(ru), method = "constant"),
         Haz0 = haz0,
         Haz0.lower = approxfun(t0.haz, hl, yleft = 0, yright = max(hl), method = "constant"),
         Haz0.upper = approxfun(t0.haz, hu, yleft = 0, yright = max(hu), method = "constant"))
}
