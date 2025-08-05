#' Functions used for nonparametric estiamtion
#' This function gives a point estimates assuming one type of event
#' @noRd
npFit0 <- function(DF, typeTem = ".", numAdj = 1e-3) {
  df0 <- DF[DF$event == 0,]
  df1 <- DF[DF$event > 0,]
  rownames(df0) <- rownames(df1) <- NULL
  m <- aggregate(event > 0~ id, data = DF, sum)[, 2]
  yi <- df0$time2
  ti <- df1$time2
  zi <- wi <- rep(1, length(ti))
  di <- df0$terminal
  xi <- as.matrix(df0[,-c(1:6)])
  p <- ncol(xi)
  yi2 <- sort(unique(yi))
  rate <- c(reRate(ti, rep(yi, m), wi, yi))
  Lam <- exp(-rate)
  keep <- !duplicated(yi)
  Lam0 <- approxfun(yi[keep], Lam[keep],
                    yleft = min(Lam), yright = max(Lam))
  if (numAdj > min(Lam)) {
    numAdj <- numAdj * min(Lam)
  }
  zi <- (m + numAdj) / (Lam + numAdj)
  if (typeTem != ".") {
    Haz <- c(temHaz(rep(0, p), rep(0, p), xi, yi, zi, di, wi, yi2))
    Haz0 <- approxfun(yi2, Haz, yleft = min(Haz), yright = max(Haz))
    return(list(Lam0 = Lam0, Haz0 = Haz0, log.muZ = log(mean(zi))))
  } else {
    return(list(Lam0 = Lam0, log.muZ = log(mean(zi))))
  }    
}

npFit <- function(DF, B = 0, typeTem = ".", model) {
  if (model == "cox.GLIPCW") out <- npFit.coxGLIPCW(DF, typeTem = typeTem)
  if (model == "cox.GLIPSW") out <- npFit.coxGLIPSW(DF, typeTem = typeTem)
  if (!(model %in% c("cox.GLIPCW", "cox.GLIPSW"))) out <- npFit0(DF, typeTem = typeTem)
  df0 <- DF[DF$event == 0,]
  mt <- aggregate(event > 0 ~ id, data = DF, sum)$event
  clsz <- mt + 1
  t0 <- sort(unique(DF$time2))
  LamB <- HazB <- matrix(NA, length(t0), B)
  if (B > 0) {
    for (i in 1:B) {
      sampled.id <- sample(df0$id, nrow(df0), TRUE)
      ind <- unlist(sapply(sampled.id, function(x) which(DF$id == x)))
      DF2 <- DF[ind,]
      DF2$id <- rep(1:nrow(df0), clsz[sampled.id])
      if (model == "cox.GLIPCW") tmp <- npFit.coxGLIPCW(DF2)
      if (model == "cox.GLIPSW") tmp <- npFit.coxGLIPSW(DF2)
      if (!(model %in% c("cox.GLIPCW", "cox.GLIPSW"))) tmp <- npFit0(DF2)
      LamB[,i] <- tmp$Lam0(t0) * exp(tmp$log.muZ)
      if (typeTem != ".") HazB[,i] <- tmp$Haz0(t0) * exp(tmp$log.muZ)
    }
    Lam.sd <- apply(LamB, 1, sd)
    if (typeTem != ".") Haz.sd <- apply(HazB, 1, sd)
    LamB.lower <- out$Lam0(t0) * exp(out$log.muZ) - 1.96 * Lam.sd
    LamB.upper <- out$Lam0(t0) * exp(out$log.muZ) + 1.96 * Lam.sd
    if (typeTem != ".") HazB.lower <- out$Haz0(t0) * exp(out$log.muZ) - 1.96 * Haz.sd
    if (typeTem != ".") HazB.upper <- out$Haz0(t0) * exp(out$log.muZ) + 1.96 * Haz.sd
    out$Lam0.lower <- approxfun(t0, LamB.lower,
                                yleft = min(LamB.lower), yright = max(LamB.lower))
    out$Lam0.upper <- approxfun(t0, LamB.upper,
                                yleft = min(LamB.upper), yright = max(LamB.upper))
    if (typeTem != ".")
      out$Haz0.lower <-
        approxfun(t0, HazB.lower, yleft = min(HazB.lower), yright = max(HazB.lower))
    if (typeTem != ".")
      out$Haz0.upper <-
        approxfun(t0, HazB.upper, yleft = min(HazB.upper), yright = max(HazB.upper))
  }
  out$Lam0 <- approxfun(t0, out$Lam0(t0) * exp(out$log.muZ),
                        yleft = min(out$Lam0(t0) * exp(out$log.muZ)),
                        yright = max(out$Lam0(t0) * exp(out$log.muZ)))
  if (typeTem != ".") 
    out$Haz0 <- approxfun(t0, out$Haz0(t0) * exp(out$log.muZ),
                          yleft = min(out$Haz0(t0) * exp(out$log.muZ)),
                          yright = max(out$Haz0(t0) * exp(out$log.muZ)))
  return(out)
}

#' @importFrom survival coxph.control
npFit.coxGLIPCW <- function(DF, typeTem = ".", numAdj = 1e-3) {
  df0 <- DF[DF$event == 0,]
  df1 <- DF[DF$event > 0,]
  rownames(df0) <- rownames(df1) <- NULL  
  yi <- df0$time2
  ti <- df1$time2
  di <- df0$terminal
  id <- df1$id
  tt <- sort(unique(ti))
  ii <- di | outer(yi, tt, ">=")
  sc <- survfit(Surv(yi, 1 - di) ~ 1)
  t1 <- summary(sc, times = matrix(tt, length(yi), length(tt), byrow = TRUE))$surv 
  t2 <- summary(sc, times = outer(yi, tt, pmin))$surv
  wij1 <- ii * matrix(t1 / t2, length(yi))
  wij1[is.na(wij1)] <- 0
  ipcw <- reRateGhosh(ti, yi, id, tt, wij1) |> as.numeric()
  if (typeTem != ".") {
    st <- coxph(Surv(yi, di) ~ 1)
    Haz0 <- basehaz(st, centered = FALSE)
    return(
      list(Lam0 = approxfun(tt, ipcw, yleft = min(ipcw), yright = max(ipcw)),
           Haz0 = approxfun(Haz0$time, Haz0$hazard,
                            yleft = min(Haz0$hazard), yright = max(Haz0$hazard)),
           log.muZ = 0)
    )
  }
  return(list(Lam0 = approxfun(tt, ipcw, yleft = min(ipcw), yright = max(ipcw)),
              log.muZ = 0))
}

npFit.coxGLIPSW <- function(DF, typeTem = ".", numAdj = 1e-3) {
  df0 <- DF[DF$event == 0,]
  df1 <- DF[DF$event > 0,]
  rownames(df0) <- rownames(df1) <- NULL  
  yi <- df0$time2
  ti <- df1$time2
  di <- df0$terminal
  id <- df1$id
  tt <- sort(unique(ti))
  st <- coxph(Surv(yi, di) ~ 1)
  wij2 <- t(outer(tt, yi, "<=") / summary(survfit(st), times = tt)$surv)
  wij2[is.na(wij2)] <- 0
  ipsw <- reRateGhosh(ti, yi, id, tt, wij2) |> as.numeric()
  if (typeTem != ".") {  
    Haz0 <- basehaz(st, centered = FALSE)
    return(
      list(Lam0 = approxfun(tt, ipsw, yleft = min(ipsw), yright = max(ipsw)),
           Haz0 = approxfun(Haz0$time, Haz0$hazard,
                            yleft = min(Haz0$hazard), yright = max(Haz0$hazard)),
           log.muZ = 0)
    )
  }
  list(Lam0 = approxfun(tt, ipsw, yleft = min(ipsw), yright = max(ipsw)), log.muZ = 0)
}

## npFitSE <- function(DF, typeRec, typeTem, par1, par2, par3, par4, zi, B) {
##     n <- length(unique(DF$id))
##     E1 <- matrix(rexp(n * B), n)
##     E2 <- matrix(rexp(n * B), n)
##     c(s1(typeRec, DF, NULL, NULL, par1, par2, E1),
##       s2(typeTem, DF, NULL, NULL, par3, par4, zi, E2))
## }


## ~1
## ## Using perturbation
## npFit <- function(DF, B = 0) {
##     df0 <- DF[DF$event == 0,]
##     df1 <- DF[DF$event == 1,]
##     rownames(df0) <- rownames(df1) <- NULL
##     m <- aggregate(event ~ id, data = DF, sum)[, 2]
##     yi <- df0$time2
##     ti <- df1$time2
##     zi <- wi <- rep(1, length(ti))
##     di <- df0$terminal
##     xi <- as.matrix(df0[,-c(1:6)])
##     p <- ncol(xi)
##     yi2 <- sort(unique(yi))
##     rate <- c(reRate(ti, rep(yi, m), wi, yi))
##     Lam <- exp(-rate)
##     keep <- !duplicated(yi)
##     Lam0 <- approxfun(yi[keep], Lam[keep],
##                       yleft = min(Lam), yright = max(Lam))
##     Haz <- c(temHaz(rep(0, p), rep(0, p), xi, yi, zi, di, wi, yi2))
##     Haz0 <- approxfun(yi2, Haz, yleft = min(Haz), yright = max(Haz))
##     zi <- (m + 0.01) / (Lam + 0.01)
##     out <- list(Lam0 = Lam0, Haz0 = Haz0, log.muZ = log(mean(zi)))
##     if (B > 0) {
##         n <- length(unique(DF$id))
##         E1 <- matrix(rexp(n * B), n)
##         E2 <- matrix(rexp(n * B), n)
##         rate <- apply(E1, 2, function(e) reRate(ti, rep(yi, m), rep(e, m), yi))
##         rate <- apply(rate, 1, quantile, c(.025, .975))
##         Lam <- exp(-rate)
##         Haz <- apply(E2, 2, function(e) temHaz(rep(0, p), rep(0, p), xi, yi, zi, di, e, yi2))
##         Haz <- apply(Haz, 1, quantile, c(.025, .975))
##         zi <- (m + 0.01) / (Lam + 0.01)
##         out$Lam0.lower <- approxfun(yi[keep], Lam[2, keep],
##                                     yleft = min(Lam[2,]), yright = max(Lam[2,]))
##         out$Lam0.upper <- approxfun(yi[keep], Lam[1, keep],
##                                     yleft = min(Lam[1,]), yright = max(Lam[1,]))
##         out$Haz0.lower <- approxfun(yi2, Haz[1,],
##                                     yleft = min(Haz[1,]), yright = max(Haz[1,]))
##         out$Haz0.upper <- approxfun(yi2, Haz[2,],
##                                     yleft = min(Haz[2,]), yright = max(Haz[2,]))
##     }
##     return(out)
## }
