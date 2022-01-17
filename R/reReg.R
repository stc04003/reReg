## globalVariables("DF") ## global variables for reReg

##############################################################################
## Functions for different models
## stdErr is estimated with resampling if method = gsc or am.xc,
##        bootstrap otherwise
##############################################################################

regFit.am.GL <- function(DF, engine, stdErr) {
  DF0 <- DF[DF$event == 0,]
  p <- ncol(DF0) - 6
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
  fit.b <- eqSolve(engine@par2, log.est, engine@solver, engine@trace)
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
  fit.a <- eqSolve(engine@par1, ghoshU2, engine@solver)
  fit.b$par <- -fit.b$par
  fit.a$par <- -fit.a$par
  out <- list(par1 = fit.a$par, par1.conv = fit.a$convergence,
              par3 = fit.b$par, par3.conv = fit.b$convergence, muZ = NA)
  out$typeRec <- engine@typeRec
  out$typeTem <- engine@typeTem
  return(out)
}

regFit.am.GL.sand <- function(DF, engine, stdErr) {
  res <- regFit(DF, engine, NULL)
  DF0 <- DF[DF$event == 0,]
  p <- ncol(DF0) - 6
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
  V <- var(t(apply(E, 2, function(x) Sn(-res$par1, -res$par3, x))))
  V1 <- V[1:p, 1:p]
  V2 <- V[1:p + p, 1:p + p]
  lmfit1 <- t(apply(Z, 2, function(x) Sn(-res$par1 + x / sqrt(n), -res$par3, rep(1, n), "s1")))
  lmfit2 <- t(apply(Z, 2, function(x) Sn(-res$par1, -res$par3 + x / sqrt(n), rep(1, n), "s2")))
  if (p == 1) {
    J1 <- coef(lm(sqrt(n) * c(lmfit1) ~ c(Z)))[-1]
    J2 <- coef(lm(sqrt(n) * c(lmfit2) ~ c(Z)))[-1]
  } else {        
    J1 <- coef(lm(sqrt(n) * lmfit1 ~ t(Z)))[-1,]
    J2 <- coef(lm(sqrt(n) * lmfit2 ~ t(Z)))[-1,]
  }
  if (qr(J1)$rank == p) aVar <- AiBAi(J1, V1)
  ## aVar <- solve(J1) %*% V1 %*% t(solve(J1))
  else aVar <- ginv(J1) %*% V1 %*% t(ginv(J1))
  if (qr(J2)$rank == p) bVar <- AiBAi(J2, V2)
  ## bVar <- solve(J2) %*% V2 %*% t(solve(J2))
  else bVar <- ginv(J2) %*% V2 %*% t(ginv(J2)) 
  aSE <- sqrt(diag(aVar))
  bSE <- sqrt(diag(bVar))
  out <- c(res, list(par1.se = aSE, par3.se = bSE, par1.vcov = aVar, par3.vcov = bVar))
  out$typeRec <- engine@typeRec
  out$typeTem <- engine@typeTem
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
  out <- list(par1 = coef(fit.coxph), par1.se = sqrt(diag(vcov(fit.coxph))))
  out$typeRec <- engine@typeRec
  out$typeTem <- engine@typeTem
  out$log.muZ <- 0
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
    approxfun(cumHaz$time, exp(-cumHaz$hazard * x), yleft = 1,
              yright = min(exp(-cumHaz$hazard * x)),
              method = "constant")(T))
  wgt <- 1 / wgt ## ifelse(wgt == 0, 1 / sort(c(wgt))[2], 1 / wgt)
  wgt <- ifelse(wgt > 1e5, 1e5, wgt)
  out <- dfsane(par = engine@par1, fn = coxGLeq, wgt = wgt, 
                X = as.matrix(X[!event, ]),
                Y = Y[!event], T = ifelse(T == Y, 1e5, T), cl = mt + 1,
                alertConvergence = FALSE, quiet = TRUE,
                control = list(trace = engine@trace))
  out <- list(par1 = out$par,
              par3 = coef(fit.coxph),
              par3.se = sqrt(diag(vcov(fit.coxph))))
  out$typeRec <- engine@typeRec
  out$typeTem <- engine@typeTem
  return(out)
}

#' This implements the CPPL from M-Y
#'
#' Need to clean this function
#' 1. Xid is an array; need to convert this to a list.
#' 2. bootstrap weight is set at rep(1 / n, n) for now.
#' 3. Can it handle more than 2 weights?
#' 4. Have weights in a list, and re-organize the weights. 
#' 
#' @noRd
regFit.cox.HH <- function(DF, engine, stdErr) {
  DF0 <- subset(DF, event == 0)
  DF1 <- subset(DF, event == 1)
  tk <- sort(unique(DF1$time2))
  uid <- sort(unique(DF$id))
  n <- length(uid)
  dNit <- do.call(rbind, lapply(split(DF, DF$id), function(d) {
    tmp <- match(d$time2[d$event > 0], tk)
    out <- rep(0, length(tk))
    out[unique(tmp)] <- table(tmp)
    return(out)
  }))
  Yit <- do.call(rbind, lapply(split(DF, DF$id), function(d) {  
    out <- rep(0, length(tk))
    out[1:findInterval(max(d$time2), tk)] <- 1
    return(out)
  }))
  Xit.ind <- do.call(rbind, lapply(split(DF, DF$id), function(d) 
    pmax(1, cut(tk, c(unique(d$time2), Inf), labels = FALSE), na.rm = TRUE)))
  xMat <- as.matrix(DF[,-(1:6)])
  ## list version
  Xit.l <- apply(xMat, 2, function(xx)
    t(sapply(uid, function(i) (xx[DF$id == i])[Xit.ind[i,]])), simplify = FALSE)
  ## array version
  Xit <- array(NA, list(nrow(Xit.ind), ncol(xMat), ncol(Xit.ind)))
  for (i in 1:ncol(xMat)) Xit[,i,] <- Xit.l[[i]]
  attr(Xit, "dimnames") <- list(NULL, names(Xit.l), NULL)
  ## Define w1 and w2
  w1t <- w2t <- wfun1 <- wfun2 <- NULL
  wfun1 <- engine@cppl.wfun[[1]]
  if (length(engine@cppl.wfun) > 1) wfun2 <- engine@cppl.wfun[[2]]
  bootstrap.weights <- rep(1 / n, n)
  ## Default EL
  ini <- CPPL.EL(dNit, Yit, Xit, NULL, NULL,
                 initial = engine@par2,
                 maxit1 = engine@maxit1, maxit2 = engine@maxit2,
                 tol = engine@tol, trace = engine@trace)$coefficient   
  S0t <- colSums(exp(apply(aperm(Xit, c(2, 1, 3)) * ini, c(2, 3), sum)) * 
                 as.vector(Yit) * bootstrap.weights)
  dLhat <- colSums(dNit * bootstrap.weights) * (S0t != 0) / (S0t + (S0t == 0))
  Lhat <- cumsum(dLhat)
  if (is.function(wfun1)) {
    w1t <- wfun1(tk)
  } else {
    if (is.character(wfun1) && wfun1 == "cumbase") w1t <- Lhat
    else if (is.character(wfun1) && wfun1 == "Gehan") w1t <- S0t
  }
  if (is.function(wfun2)) {
    w2t <- wfun2(tk)
  } else {if (is.character(wfun2) && wfun2 == "cumbase") w2t <- Lhat
          else if (is.character(wfun2) && wfun2 == "Gehan") w2t <- S0t
  }
  if (engine@cppl == "EL")
    out <- CPPL.EL(dNit, Yit, Xit, w1t, w2t,
                   initial = engine@par2,
                   maxit1 = engine@maxit1, maxit2 = engine@maxit2,
                   tol = engine@tol, trace = engine@trace)
  if (engine@cppl == "GMM")
    out <- CPPL.GMM(dNit, Yit, Xit, w1t, w2t,
                    initial = engine@par2,
                    maxit1 = engine@maxit1, maxit2 = engine@maxit2,
                    tol = engine@tol, trace = engine@trace)
  out$Lam0 <- function(x)
    approx(x = tk, y = cumsum(out$diff.baseline), xout = x, 
           yleft = min(out$diff.baseline), yright = sum(out$diff.baseline))$y
  out$par1 <- as.numeric(out$coefficient)
  out$par1.vcov <- out$covariance
  out$par1.se <- sqrt(diag(out$par1.vcov))
  dimnames(out$par1.vcov) <- out$coefficient <- out$covariance <- NULL
  out$typeRec <- engine@typeRec
  out$typeTem <- engine@typeTem
  out$log.muZ <- 0
  return(out)
}

## #' @importFrom rlang is_empty
regFit.general <- function(DF, engine, stdErr) {
  eqList <- c("dfsane", "BBsolve", "optim", "BBoptim", "optimr", "hjk", "mads", "nleqslv")
  if (is.na(match(engine@solver, eqList))) {
    warning("Unidentified solver; BB::dfsane is used.")
    engine@solver <- "dfsane"
  }
  out <- s1(engine@typeRec, DF, engine@eqType, engine@solver, engine@par1, engine@par2,
            trace = engine@trace)
  if (engine@typeTem != ".") 
    out <- c(out, s2(engine@typeTem, DF, engine@eqType, engine@solver,
                     engine@par3, engine@par4, out$zi, trace = engine@trace))
  out$typeRec <- engine@typeRec
  out$typeTem <- engine@typeTem
  return(out)
}

s1 <- function(type, DF, eqType, solver, par1, par2,
               Lam0 = NULL, w1 = NULL, w2 = NULL, trace = FALSE) {
  if (type == "gsc") return(reSC(DF, eqType, solver, par1, par2, Lam0, w1, w2, trace))
  if (type == "cox") return(reCox(DF, eqType, solver, par1, Lam0, w1, trace))
  if (type == "am") return(reAM(DF, eqType, solver, par1, Lam0, w1, trace))
  if (type == "ar") return(reAR(DF, eqType, solver, par1, Lam0, w1, trace))
  return(NULL)
}

s2 <- function(type, DF, eqType, solver, par3, par4, zi, wgt = NULL, trace = FALSE) {
  if (type == "gsc") return(temSC(DF, eqType, solver, par3, par4, zi, wgt, trace))
  if (type == "cox") return(temCox(DF, eqType, solver, par3, zi, wgt, trace))
  if (type == "am") return(temAM(DF, eqType, solver, par3, zi, wgt, trace))
  if (type == "ar") return(temAR(DF, eqType, solver, par3, zi, wgt, trace))
  return(NULL)
}

regFit.general.sand <- function(DF, engine, stdErr) {
  if (is.na(match(engine@solver, c("dfsane", "BBsolve", "optim", "BBoptim")))) {
    warning("Unidentified solver; BB::dfsane is used.")
    engine@solver <- "dfsane"
  }
  res <- regFit(DF, engine, NULL)
  n <- length(unique(DF$id))
  B <- stdErr@B
  p <- ncol(DF) - 6
  tmpV <- replicate(B,
                    c(s1(engine@typeRec, DF, engine@eqType, NULL, res$par1, res$par2,
                         res$Lam0, rexp(n), rexp(n))$value,
                      s2(engine@typeTem, DF, engine@eqType, NULL, res$par3, res$par4,
                         rexp(n) * res$zi, rexp(n))))    
  V <- var(t(tmpV))
  Z <- matrix(rnorm(ncol(V) * B), B)
  len1 <- length(res$par1)
  len2 <- length(res$par2)
  len3 <- length(res$par3)
  len4 <- length(res$par4)
  na <- len1 + len2
  nb <- len3 + len4
  L <- apply(Z, 1, function(zz) {
    c(s1(engine@typeRec, DF, engine@eqType, NULL,
         res$par1 + zz[1:len1] / sqrt(n), res$par2 + zz[1:len2 + len1] / sqrt(n), res$Lam0)$value,
      s2(engine@typeTem, DF, engine@eqType, NULL,
         res$par3 + zz[1:len3 + len1 + len2] / sqrt(n),
         res$par4 + zz[1:len4 + len1 + len2 + len3] / sqrt(n), res$zi))})
  L <- t(L)
  ## if (engine@typeRec == "gsc") {
  ##     J <- rbind(cbind(t(Axb(Z[,1:len1], sqrt(n) * L[,1:len1])), matrix(0, len1, len2)),
  ##                t(Axb(Z, sqrt(n) * L[,1:len2 + len1])))        
  ## } else
  J <- t(Axb(Z, sqrt(n) * L))
  recVar <- AiBAi(J[1:na, 1:na], V[1:na, 1:na])
  par1.vcov <- recVar[1:len1, 1:len1, drop = FALSE]
  par1.se <- sqrt(diag(par1.vcov))
  res <- c(res, list(par1.vcov = par1.vcov, par1.se = par1.se))
  if (len2 > 0) {
    A <- cbind(diag(len1), 0, diag(len2 - 1))
    par2.vcov <- A %*% recVar %*% t(A)
    par2.se <- sqrt(diag(par2.vcov))
    res <- c(res, list(par2.vcov = par2.vcov, par2.se = par2.se))
    res$vcovRec12 <- recVar[1:len1, 2:len2 + len1]
  }
  if (nb > 0) {
    if (engine@typeTem == "gsc") {
      bCoef <- matrix(0, B, nb)
      convergence <- rep(0, B)
      id <- DF$id
      uID <- unique(id)
      mt <- aggregate(event ~ id, data = DF, sum)$event
      clsz <- mt + 1
      for (i in 1:B) {
        sampled.id <- unlist(lapply(split(sort(uID), mt > 0), function(x) {
          if (length(x) == 1) return(x)
          else return(sample(x, replace = TRUE))
        }))
        names(sampled.id) <- NULL
        ind <- unlist(sapply(sampled.id, function(x) which(id == x)))
        DF2 <- DF[ind,]
        DF2$id <- rep(1:n, clsz[sampled.id])
        tmp <- temSC(DF2, engine@eqType, engine@solver, res$par3, res$par4, res$zi)
        bCoef[i,] <- c(tmp$par3, tmp$par4)
      }
      tmp <- apply(bCoef, 1, crossprod)
      ps <- c(length(res$par3), length(res$par4))
      ps <- ps[ps > 0]
      bound <- c(res$par3, res$par4)
      convergence <- rowSums(sapply(1:length(ps), function(x) {
        p <- rev(seq(sum(ps[1:x])))[1:ps[x]]
        apply(bCoef[,p, drop = FALSE], 1, crossprod) > 10 * drop(crossprod(bound[p]))
      }))
      converged <- which(convergence == 0)
      if (sum(convergence) > 0) {
        message("Some bootstrap samples failed to converge")
        rm <- unique(c(which(convergence > 0),
                       which(tmp %in% boxplot(tmp, plot = FALSE)$out)))
        converged <- (1:B)[-rm]
      }
      temVar <- var(bCoef[converged, ], na.rm = TRUE)
      res$par3.vcov <- temVar[1:len3, 1:len3, drop = FALSE]
      res$par3.se <- sqrt(diag(res$par3.vcov))            
      res$par4.vcov <- temVar[1:len4 + len3, 1:len4 + len3, drop = FALSE]
      res$par4.se <- sqrt(diag(res$par4.vcov))
      res$vcovTem12 <- temVar[1:len3, 1:len4 + len3, drop = FALSE]
    } else { 
      ind2 <- tail(1:nrow(J), nb)
      J2 <- Axb(Z[,ind2], sqrt(n) * L[,ind2])
      temVar <- AiBAi(J2, V[ind2, ind2])
      res$par3.vcov <- temVar[1:len3, 1:len3, drop = FALSE]
      res$par3.se <- sqrt(diag(res$par3.vcov))
      if (!is.null(res$par4)) {
        res$par4.vcov <- temVar[1:len4 + len3, 1:len4 + len3, drop = FALSE]
        res$par4.se <- sqrt(diag(res$par4.vcov))
      }
      res$vcovTem12 <- temVar[1:len3, 2:len4 + len3, drop = FALSE]
    }
  }
  return(res)
}

##############################################################################
                                        # Variance estimation 
##############################################################################
regFit.Engine.boot <- function(DF, engine, stdErr) {
  res <- regFit(DF, engine, NULL)
  id <- DF$id
  uID <- unique(id)
  event <- DF$event
  status <- DF$terminal
  X <- as.matrix(DF[,-c(1:6)])    
  n <- length(unique(id))
  T <- DF$time2
  mt <- aggregate(event ~ id, data = DF, sum)$event
  clsz <- mt + 1
  Y <- rep(DF$time2[event == 0], mt + 1)
  cluster <- unlist(sapply(mt + 1, function(x) 1:x))
  B <- stdErr@B
  bound <- c(res$par1, res$par2, res$par3, res$par4)
  engine2 <- engine
  engine2@par1 <- res$par1
  if (!is.null(res$par2)) engine2@par2 <- res$par2
  if (stdErr@parallel) {
    cl <- makeCluster(stdErr@parCl)
    clusterExport(cl = cl,
                  varlist = c("DF", "engine2"),
                  envir = environment())
    out <- parSapply(cl, 1:B, function(x) {
      sampled.id <- sample(unique(id), n, TRUE)
      ind <- unlist(sapply(sampled.id, function(x) which(id == x)))
      DF2 <- DF[ind,]
      DF2$id <- rep(1:n, clsz[sampled.id])
      tmp <- regFit(DF2, engine2, NULL)
      return(c(tmp$par1, tmp$par2, tmp$par3, tmp$par4))
    })
    stopCluster(cl)
    bCoef <- t(out)
  } else {
    bCoef <- matrix(0, B, length(bound))
    convergence <- rep(0, B)
    for (i in 1:B) {
      ## sampled.id <- sample(unique(id), n, TRUE)
      sampled.id <- unlist(lapply(split(sort(uID), mt > 0), function(x) {
        if (length(x) == 1) return(x)
        else return(sample(x, replace = TRUE))
      }))
      names(sampled.id) <- NULL
      ind <- unlist(sapply(sampled.id, function(x) which(id == x)))
      DF2 <- DF[ind,]
      DF2$id <- rep(1:n, clsz[sampled.id])
      tmp <- regFit(DF2, engine2, NULL)
      bCoef[i,] <- c(tmp$par1, tmp$par2, tmp$par3, tmp$par4)
    }
  }
  tmp <- apply(bCoef, 1, crossprod)
  ps <- c(length(res$par1), length(res$par2), length(res$par3), length(res$par4))
  ps <- ps[ps > 0]
  convergence <- rowSums(sapply(1:length(ps), function(x) {
    p <- rev(seq(sum(ps[1:x])))[1:ps[x]]
    apply(bCoef[,p,drop = FALSE], 1, crossprod) > 10 * drop(crossprod(bound[p]))
  }))
  converged <- which(convergence == 0)
  ## res$bCoef <- bCoef[converged,]
  if (sum(convergence) > 0) {
    message("Some bootstrap samples failed to converge")
    rm <- unique(c(which(convergence > 0),
                   which(tmp %in% boxplot(tmp, plot = FALSE)$out)))
    converged <- (1:B)[-rm]
  }
  ## if (all(convergence != 0) || sum(convergence == 0) == 1) {
  ##     message("Some bootstrap samples failed to converge")
  ##     converged <- (1:B)[- which(tmp %in% boxplot(tmp, plot = FALSE)$out)]
  ## }
  bVar <- var(bCoef[converged, ], na.rm = TRUE)
  bSE <- sqrt(diag(as.matrix(bVar)))
  len1 <- length(res$par1)
  len2 <- length(res$par2)
  len3 <- length(res$par3)
  len4 <- length(res$par4)
  res <- c(res, list(par1.vcov = bVar[1:len1, 1:len1, drop = FALSE],
                     par1.se = bSE[1:len1], B = length(converged)))
  if (len2 > 0)
    res <- c(res, list(par2.vcov = bVar[1:len2 + len1, 1:len2 + len1, drop = FALSE],
                       par2.se = bSE[1:len2 + len1]))
  if (len3 > 0)
    res <- c(res, list(par3.vcov = bVar[1:len3 + len1 + len2, 1:len3 + len1 + len2, drop = FALSE],
                       par3.se = bSE[1:len3 + len1 + len2]))
  if (len4 > 0)
    res <- c(res, list(par4.vcov = bVar[1:len4 + len1 + len2 + len3,
                                        1:len4 + len1 + len2 + len3, drop = FALSE],
                       par4.se = bSE[1:len4 + len1 + len2 + len3]))
  if (engine@typeRec == "gsc") {        
    res$par2.vcov <- var(bCoef[converged, 1:len1, drop = FALSE] +
                         bCoef[converged, 2:len2 + len1, drop = FALSE])
    res$par2.se <- sqrt(diag(res$par2.vcov))  
    res$vcovRec12 <- bVar[1:len1, 2:len2 + len1, drop = FALSE]
  }
  if (engine@typeTem == "gsc")
    res$vcovTem12 <- bVar[1:len3 + len1 + len2, 1:len4 + len1 + len2 + len3, drop = FALSE]
  return(res)
}

##############################################################################
                                        # Class Definition
##############################################################################

setClass("Engine",
         representation(tol = "numeric",
                        par1 = "numeric", par2 = "numeric",
                        par3 = "numeric", par4 = "numeric",
                        baseSE = "logical", 
                        solver = "character", eqType = "character", 
                        typeRec = "character", typeTem = "character",
                        maxit1 = "numeric", maxit2 = "numeric",
                        trace = "logical"),
         prototype(eqType = "logrank", tol = 1e-7,
                   par1 = 0, par2 = 0, par3 = 0, par4 = 0, 
                   baseSE = FALSE, solver = "dfsane",
                   maxit1 = 100, maxit2 = 10, trace = FALSE),
         contains = "VIRTUAL")
setClass("general", contains = "Engine")
setClass("cox.LWYY", contains = "Engine")
setClass("cox.HH",
         representation(cppl.wfun = "list", cppl = "character"),
         prototype(cppl.wfun =  list(NULL, NULL), cppl = "EL"),
         contains = "Engine")
setClass("am.GL", contains = "Engine")
## setClass("gsc.XCYH", representation(muZ = "numeric"), prototype(muZ = 0), contains = "Engine")
setClass("cox.GL",
         representation(wgt = "matrix"), prototype(wgt = matrix(0)), contains = "Engine")

setClass("stdErr",
         representation(B = "numeric", parallel = "logical", parCl = "numeric"),
         prototype(B = 100, parallel = FALSE, parCl = parallel::detectCores() / 2L),
         contains = "VIRTUAL")

setClass("boot", contains = "stdErr")
setClass("sand", contains = "stdErr")


##############################################################################
                                        # Method Dispatch
##############################################################################
setGeneric("regFit", function(DF, engine, stdErr) {standardGeneric("regFit")})

setMethod("regFit", signature(engine = "general", stdErr = "NULL"), regFit.general)
setMethod("regFit", signature(engine = "general", stdErr = "sand"), regFit.general.sand)
setMethod("regFit", signature(engine = "cox.LWYY", stdErr = "NULL"), regFit.cox.LWYY)
setMethod("regFit", signature(engine = "cox.LWYY", stdErr = "boot"), regFit.cox.LWYY)
setMethod("regFit", signature(engine = "cox.LWYY", stdErr = "sand"), regFit.cox.LWYY)
setMethod("regFit", signature(engine = "cox.HH", stdErr = "NULL"), regFit.cox.HH)
setMethod("regFit", signature(engine = "cox.HH", stdErr = "boot"), regFit.cox.HH)
setMethod("regFit", signature(engine = "cox.HH", stdErr = "sand"), regFit.cox.HH)
setMethod("regFit", signature(engine = "cox.GL", stdErr = "NULL"), regFit.cox.GL)
setMethod("regFit", signature(engine = "cox.GL", stdErr = "sand"), regFit.cox.GL)
setMethod("regFit", signature(engine = "am.GL", stdErr = "NULL"), regFit.am.GL)
setMethod("regFit", signature(engine = "Engine", stdErr = "boot"),
          regFit.Engine.boot)
setMethod("regFit", signature(engine = "am.GL", stdErr = "sand"),
          regFit.am.GL.sand)


#' Fits Semiparametric Regression Models for Recurrent Event Data
#'
#' Fits a general (joint) semiparametric regression model for the recurrent event data,
#' where the rate function of the underlying recurrent event process and
#' the hazard function of the terminal event can be specified as a Cox-type model,
#' an accelerated mean model, an accelerated rate model, or a generalized scale-change model.
#' See details for model specifications.
#'
#' @details
#'
#' \bold{Model specification:}
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
#' \eqn{(\alpha, \eta)} and \eqn{(\beta, \theta)} correspond to the shape and size parameters,
#' respectively.
#' The model includes several popular semiparametric models as special cases,
#' which can be specified via the \code{model} argument with the rate function
#' and the hazard function separated by "\code{|}".
#' For examples,
#' Wang, Qin and Chiang (2001) (\eqn{\alpha = \eta = \theta = 0})
#' can be called with \code{model = "cox"};
#' Huang and Wang (2004) (\eqn{\alpha = \eta = 0})
#' can be called with \code{model = "cox|cox"};
#' Xu et al. (2017) (\eqn{\alpha = \beta} and \eqn{\eta = \theta})
#' can be called with \code{model = "am|am"};
#' Xu et al. (2019) (\eqn{\eta = \theta = 0}) can be called with \code{model = "gsc"}.
#' Users can mix the models depending on the application. For example,
#' \code{model = "cox|ar"} postulate a Cox proportional model for the
#' recurrent event rate function and an accelerated rate model for
#' the terminal event hazard function (\eqn{\alpha = \theta = 0}).
#' If only one model is specified without an "\code{|}",
#' it is used for both the rate function and the hazard function.
#' For example, specifying \code{model = "cox"} is equivalent to \code{model = "cox|cox"}.
#' Some models that assumes \code{Z = 1} and requires independent
#' censoring are also implemented in \code{reReg};
#' these includes \code{model = "cox.LWYY"} for Lin et al. (2000),
#' \code{model = "cox.GL"} for Ghosh and Lin (2002),
#' and \code{model = "am.GL"} for Ghosh and Lin (2003).
#' Additionally, an improved estimation of the proportional rate model
#' (Huang and Huang 2022) can be called by \code{model = "cox.HH"} with
#' additional \code{control} options to specify the underlying procedure.
#' See \href{https://www.sychiou.com/reReg/articles/reReg-reg.html}{online vignette}
#' for a detailed discussion of the implemented regression models.
#' 
#' \bold{Variance estimation:}
#' 
#' The available methods for variance estimation are:
#' \describe{
#'   \item{boot}{performs nonparametric bootstrap.}
#'   \item{sand}{performs the efficient resampling-based variance estimation.}
#' }
#'
#' \bold{Improving proportional rate model:}
#' A common semiparametric regression model for recurrent event process
#' under the noninformative censoring assumption is the Cox-type proportional rate model
#' (available in \code{reReg()} via \code{model = "cox.LWYY"}).
#' However, the construction of the pseudo-partial score function ignores the
#' dependency among recurrent events and thus could be inefficient. 
#' To improve upon this popular method, Huang and Huang (2022) proposed to combine
#' a system of weighted pseudo-partial score equations via the generalized method of moments (GMM)
#' and empirical likelihood (EL) estimation.
#' The proposed GMM and EL procedures are available in \code{reReg} via \code{model = "cox.HH"}
#' with additional control specifications.
#' See \href{https://www.sychiou.com/reReg/articles/reReg-cppl.html}{online vignette}
#' for an illustration of this feature.
#' 
#' \bold{Control options:}
#' 
#' The \code{control} list consists of the following parameters:
#' \describe{
#'   \item{tol}{absolute error tolerance.}
#'   \item{init}{a list contains initial guesses used for root search.}
#'   \item{solver}{the equation solver used for root search.
#' The available options are \code{BB::BBsolve}, \code{BB::dfsane}, \code{BB::BBoptim}, 
#' \code{optimx::optimr}, \code{dfoptim::hjk}, \code{dfoptim::mads}, \code{optim},
#' and \code{nleqslv::nleqslv}.}
#'   \item{eqType}{a character string indicating whether the log-rank type estimating equation or
#' the Gehan-type estimating equation (when available) will be used. }
#'   \item{boot.parallel}{an logical value indicating whether parallel computation
#' will be applied when \code{se = "boot"} is called.}
#'   \item{boot.parCl}{an integer value specifying the number of CPU cores to be used when
#' \code{parallel = TRUE}. The default value is half the CPU cores on the current host.}
#' \item{cppl}{A character string indicating either to improve the proportional rate model via
#' the generalized method of moments (\code{cppl = "GMM"}) or empirical likelihood estimation (\code{cppl = "EL"}).
#' This option is only used when \code{model = "cox.HH"}.}
#' \item{cppl.wfun}{A list of (up to two) weight functions to be combined with the weighted pseudo-partial likelihood scores.
#' Avaialble options are \code{"Gehan"} and \code{"cumbase"},
#' which correspond to the Gehan's weight and the cumulative baseline hazard function, respectively.
#' Alternatively, the weight functions can be specified with function formulas.
#' This option is only used when \code{model = "cox.HH"}.}
#' \item{trace}{A logical variable denoting whether some of the
#' intermediate results of iterations should be displayed to the user.  Default is \code{FALSE}.}
#' }
#' 
#' @param formula a formula object, with the response on the left of a "~" operator,
#' and the predictors on the right.
#' The response must be a recurrent event survival object as returned by function \code{Recur}.
#' @param data  an optional data frame in which to interpret the variables occurring
#' in the \code{"formula"}.
#' @param subset an optional logical vector specifying a subset of observations to be used
#' in the fitting process.
#' @param B a numeric value specifies the number of bootstraps for variance estimation.
#' When \code{B = 0}, variance estimation will not be performed.
#' @param model a character string specifying the underlying model.
#' The available functional form for the rate function and the hazard function include a Cox-type model,
#' an accelerated mean model, an accelerated rate model, or a generalized scale-change model,
#' and can be specified via "cox", "am", "ar", or "gsc", respectively.
#' The rate function and hazard function separated by "\code{|}".
#' See \bold{Details}.
#' @param se a character string specifying the method for the variance estimation. See \bold{Details}.
#' \describe{
#'    \item{\code{boot}}{ nonparametric bootstrap approach}
#'    \item{\code{sand}}{ resampling-based sandwich estimator}
#' }
#' @param control a list of control parameters. See \code{\link{reReg.control}} for default values.
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
#' Processes under Informative Censoring. \emph{Statistica Sinica}, \bold{30}: 1773--1795.
#' @references Huang, M.-Y. and Huang, C.-Y. (2022). Improved semiparametric estimation of the proportional rate model with recurrent event data.
#' \emph{In revision}.
#'
#' @importFrom stats approxfun optim model.response 
#' @importFrom stats .getXlevels 
#' @seealso \code{\link{Recur}}, \code{\link{simGSC}}
#'
#' @example inst/examples/ex_reReg.R

reReg <- function(formula, data, subset,
                  model = "cox",
                  B = 0, se = c("boot", "sand"),
                  control = list()) {
  ## se = c("resampling", "bootstrap", "NULL"),
  ## se <- ifelse(is.null(se), "NULL", se)
  se <- match.arg(se)
  Call <- match.call()
  if (missing(formula)) stop("Argument 'formula' is required.")
  if (missing(data)) 
    data <- environment(formula)
  if (!missing(subset)) {
    sSubset <- substitute(subset)
    subIdx <- eval(sSubset, data, parent.frame())
    if (!is.logical(subIdx)) 
      stop("'subset' must be logical")
    subIdx <- subIdx & !is.na(subIdx)
    data <- data[subIdx, ]
  }    
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$data <- data
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mm <- stats::model.matrix(formula, data = mf)
  obj <- stats::model.extract(mf, "response")
  DF <- cbind(obj, mm)
  DF <- as.data.frame(DF)
  DF <- DF[,colnames(DF) != "(Intercept)"]
  if (!is.Recur(obj)) stop("Response must be a `Recur` object")
  formula[[2]] <- NULL
  if (formula == ~ 1) DF$zero = 0 
  ctrl <- reReg.control()
  if (!is.null(control$init)) {
    if(is.null(control$init$alpha)) control$par1 <- 0
    else control$par1 <- control$init$alpha
    if(is.null(control$init$beta)) control$par2 <- 0
    else control$par2 <- control$init$beta
    if(is.null(control$init$eta)) control$par3 <- 0
    else control$par3 <- control$init$eta
    if(is.null(control$init$theta)) control$par4 <- 0
    else control$par4 <- control$init$theta
    control$init <- NULL
  }
  namc <- names(control)
  if (!all(namc %in% names(ctrl))) 
    stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
  ctrl[namc] <- control
  if (!is.null(ctrl$cppl)) model <- "cox.HH"
  DF <- DF[order(DF$id, DF$time2), ]
  allModel <- apply(expand.grid(c("cox", "am", "gsc", "ar"),
                                c("cox", "am", "gsc", "ar", ".")), 1, paste, collapse = "|")
  allModel <- c(allModel, "cox.LWYY", "cox.HH", "cox.GL", "cox.HW", "am.GL", "am.XCHWY", "gsc.XCYH")
  model <- match.arg(model, c("cox", "am", "gsc", "ar", allModel))
  typeRec <- typeTem <- NULL
  if (grepl("|", model, fixed = TRUE)) {
    typeRec <- substring(model, 1, regexpr("[|]", model) - 1)
    typeTem <- substring(model, regexpr("[|]", model) + 1)
    model <- "general"
  }
  if (model %in% c("cox", "am", "gsc", "ar")) {
    typeRec <- model
    typeTem <- "."
    model <- "general"
  }
  ## Special cases:
  if (model == "cox.HW") {
    typeRec <- typeTem <- "cox"
    model <- "general"
  }
  if (model == "am.XCHWY") {
    typeRec <- typeTem <- "am"
    model <- "general"
  }
  if (model == "gsc.XCYH") {
    typeRec <- "gsc"
    typeTem <- "."
    model <- "general"        
  }
  if (model == "cox.LWYY") {
    typeRec <- "cox.LWYY"
    typeTem <- "."
  }
  if (model == "cox.HH") {
    typeRec <- "cox.HH"
    typeTem <- "."
  }
  if (model == "cox.GL") typeRec <- typeTem <- "cox.GL"
  if (model == "am.GL") typeRec <- typeTem <- "am.GL"
  if (length(unique(DF$time2[DF$event == 0])) == 1 & typeTem != ".") {
    typeTem <- "."
    message("Only one unique censoring time is detected, terminal event model is not fitted.")
  }
  ## Temporary fix 
  if (typeRec != "gsc")  se <- "boot"
  engine.ctrl <- ctrl[names(ctrl) %in% names(attr(getClass(model), "slots"))]
  engine <- do.call("new", c(list(Class = model), engine.ctrl))
  engine@typeRec <- typeRec
  engine@typeTem <- typeTem
  if (se == "NULL" || B == 0)
    stdErr <- NULL
  else {
    stdErr.ctrl <- ctrl[names(ctrl) %in% names(attr(getClass(se), "slots"))]
    stdErr <- do.call("new", c(list(Class = se), stdErr.ctrl))
    stdErr@B <- B
  }
  ## initial values
  p <- ncol(DF) - ncol(mf[[1]])
  if (model %in% c("cox.GL", "am.GL")) {
    if (length(engine@par1) == 1) engine@par1 <- rep(engine@par1, p)
    if (length(engine@par1) != p)
      stop("The length of initial value does not match with the number of covariates.")
  }
  if (model == "general") {
    if (typeRec == "cox") {
      if (length(engine@par1) == 1) engine@par1 <- rep(engine@par1, p + 1)
      if (length(engine@par1) == p) engine@par1 <- c(0, engine@par1)
      if (length(engine@par1) != (p + 1))
        stop("The length of initial value does not match with the number of covariates.")
      if (typeTem != ".") {
        engine@par3 <- engine@par2
        if (length(engine@par3) == 1) engine@par3 <- rep(engine@par3, p)
        if (length(engine@par3) != p)
          stop("The length of initial value does not match with the number of covariates.")
        engine@par4 <- engine@par2
        if (length(engine@par4) == 1) engine@par4 <- rep(engine@par4, p)
        if (length(engine@par4) != p)
          stop("The length of initial value does not match with the number of covariates.")
      }
    }
    if (typeRec == "gsc") {
      if (length(engine@par1) == 1) engine@par1 <- rep(engine@par1, p)
      if (length(engine@par2) == 1) engine@par2 <- rep(engine@par2, p + 1)
      if (length(engine@par2) == p) engine@par2 <- c(0, engine@par2)
      if (length(engine@par1) != p | length(engine@par2) != (p + 1))
        stop("The length of initial value does not match with the number of covariates.")
      if (typeTem != ".") {
        if (length(engine@par3) == 1) engine@par3 <- rep(engine@par3, p)
        if (length(engine@par4) == 1) engine@par4 <- rep(engine@par4, p)
        if (length(engine@par3) != p | length(engine@par4) != p)
          stop("The length of initial value does not match with the number of covariates.")
      }
    }
    if (typeRec %in% c("ar", "am")) {
      if (length(engine@par1) == 1) engine@par1 <- rep(engine@par1, p)
      if (length(engine@par1) != p)
        stop("The length of initial value does not match with the number of covariates.")
      if (typeTem != ".") {
        engine@par3 <- engine@par2
        if (length(engine@par3) == 1) engine@par3 <- rep(engine@par3, p)
        if (length(engine@par3) != p)
          stop("The length of initial value does not match with the number of covariates.")
        if (length(engine@par4) == 1) engine@par4 <- rep(engine@par4, p)
        if (length(engine@par4) != p)
          stop("The length of initial value does not match with the number of covariates.")
      }
    }
  }
  engine@baseSE <- B > 0
  if (formula == ~1) {
    if (engine@baseSE) fit <- npFit(DF, B, typeTem)
    else fit <- npFit(DF, 0, typeTem)
    fit$typeRec <- "nonparametric"
    fit$typeTem <- typeTem
  } else {
    fit <- regFit(DF = DF, engine = engine, stdErr = stdErr)
  }    
  fit$DF <- DF
  fit$call <- Call
  fit$varNames <- names(DF)[-(1:6)]
  fit$se <- se
  fit$xlevels <- .getXlevels(attr(mf, "terms"), mf)
  if (engine@typeRec == "cox") fit$par1 <- fit$par1[-1]
  if (engine@typeRec == "gsc") fit$par2 <- fit$par1 + fit$par2[-1]
  ## if (se != "NULL" & se != "boot" & engine@typeRec == "gsc") {
  ##     fit$par2.se <- fit$par2.se[-1]
  ##     fit$par2.vcov <- fit$par2.vcov[-1, -1, drop = FALSE]
  ## }
  if (se != "NULL" & engine@typeRec == "cox") fit$par1.se <- fit$par1.se[-1]
  fit <- fit[order(names(fit))]
  class(fit) <- "reReg"
  return(fit)
}

#' Equation wrapper
#'
#' @noRd
#' @importFrom BB spg
#' @importFrom optimx optimr
#' @importFrom dfoptim hjk mads nmk
#' @importFrom rootSolve uniroot.all
#' @keywords internal
eqSolve <- function(par, fn, solver, trace, ...) {
  if (length(fn(par, ...)) == 1) {
    tmp <- uniroot.all(Vectorize(fn), interval = c(par - 10, par + 10))
    out <- NULL
    out$par <- tmp[which.min(abs(tmp - par))]
    out$convergence <- 0 ## 1 * (abs(fn(out$par, ...)) < 1e-5)
    return(out)
  }
  if (solver == "dfsane") {
    out <- dfsane(par = par, fn = function(z) fn(z, ...), 
                  alertConvergence = FALSE, quiet = TRUE,
                  control = list(trace = trace))
    if (max(abs(out$par)) > 10) solver <- "BBsolve"
  }
  if (solver == "BBsolve")
    out <- BBsolve(par = par, fn = fn, ..., quiet = TRUE)
  if (solver == "nleqslv") {
    out <- nleqslv(x = par, fn = fn, ...)
    out$par <- out$x
  }
  if (solver == "BBoptim")
    out <- BBoptim(par = par, fn = function(z) sum(fn(z, ...)^2),
                   quiet = TRUE, control = list(trace = trace))
  if (solver == "optim")
    out <- optim(par = par, fn = function(z) sum(fn(z, ...)^2),
                 control = list(trace = trace))
  if (solver == "optimr")
    out <- optimr(par = par, fn = function(z) sum(fn(z, ...)^2))
  if (solver == "hjk")
    out <- hjk(par = par, fn = function(z) sum(fn(z, ...)^2))
  if (solver == "mads")
    out <- mads(par = par, fn = function(z) sum(fn(z, ...)^2),
                control = list(trace = trace))
  if (solver == "nmk")
    out <- nmk(par = par, fn = function(z) sum(fn(z, ...)^2))    
  return(out)
}

#' Package options for reReg
#'
#' This function provides the fitting options for the \code{reReg()} function. 
#'
#' @param eqType a character string indicating whether the log-rank type estimating equation
#' or the Gehan-type estimating equation (when available) will be used.
#' @param solver a character string specifying the equation solver to be used for root search.
#' @param tol a numerical value specifying the absolute error tolerance in root search.
#' @param init a list contains the initial guesses used for root search.
#' @param boot.parallel an logical value indicating whether parallel computation will be
#' applied when \code{se = "boot"} is specified in \code{reReg()}.
#' @param boot.parCl an integer value specifying the number of CPU cores to be used when
#' \code{parallel = TRUE}. The default value is half the CPU cores on the current host.
#' @param cppl a character string indicating either to improve the proportional rate model via
#' the generalized method of moments (\code{cppl = "GMM"}) or empirical likelihood estimation (\code{cppl = "EL"}).
#' This option is only used when \code{model = "cox.HH"}.
#' @param cppl.wfun A list of (up to two) weight functions to be combined with the weighted pseudo-partial likelihood scores.
#' Avaialble options are \code{"Gehan"} and \code{"cumbase"},
#' which correspond to the Gehan's weight and the cumulative baseline hazard function, respectively.
#' Alternatively, the weight functions can be specified with function formulas.
#' This option is only used when \code{model = "cox.HH"}.
#' @param maxit1,maxit2 max number of iteration used when \code{model = "cox.HH"}.
#' @param trace A logical variable denoting whether some of the
#' intermediate results of iterations should be displayed to the user.  Default is \code{FALSE}.
#' 
#' @seealso \code{\link{reReg}}
#' @export
reReg.control <- function(eqType = c("logrank", "gehan", "gehan_s"),
                          solver = c("BB::dfsane", "BB::BBsolve", "BB::BBoptim", "optimx::optimr",
                                     "dfoptim::hjk", "dfoptim::mads", "optim", "nleqslv::nleqslv"),
                          tol = 1e-7, cppl = NULL, cppl.wfun = list(NULL, NULL),
                          init = list(alpha = 0, beta = 0, eta = 0, theta = 0),
                          boot.parallel = FALSE, boot.parCl = NULL,
                          maxit1 = 100, maxit2 = 10, trace = FALSE) {
  if (is.null(boot.parCl)) boot.parCl <- parallel::detectCores() / 2L
  solver <- match.arg(solver)
  if (solver == "nleqslv::nleqslv") solve <- "nleqslv"
  if (solver == "BB::dfsane") solver <- "dfsane"
  if (solver == "BB::BBsolve") solver <- "BBsolve"
  if (solver == "BB::BBoptim") solver <- "BBoptim"
  if (solver == "optimx::optimr") solver <- "optimr"
  if (solver == "dfoptim::hjk") solver <- "hjk"
  if (solver == "dfoptim::mads") solver <- "mads"
  eqType <- match.arg(eqType)
  if (!is.null(cppl) && !(cppl %in% c("GMM", "EL"))) stop("Invalid 'cppl' method.")
  list(tol = tol, eqType = eqType, solver = solver, cppl = cppl, cppl.wfun = cppl.wfun,
       par1 = init$alpha, par2 = init$beta, par3 = init$eta, par4 = init$theta,
       parallel = boot.parallel, parCl = boot.parCl, maxit1 = maxit1, maxit2 = maxit2, trace = trace)
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
