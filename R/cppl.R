#' Functions adopted from M-Y
#'
#' @noRd
CPPL.EL <- function(dNit, Yit, Xit, w1t, w2t,
                    initial = NULL, bootstrap.weights = NULL, 
                    maxit1 = 100, maxit2 = 10, tol = 1e-7, trace = FALSE) {
  n <- nrow(dNit)
  p <- dim(Xit)[2]
  if (is.null(initial) | initial == 0) initial <- rep(0, p)
  if (is.null(bootstrap.weights)) bootstrap.weights <- rep(1 / n, n)
  sfun <- function(b, weights = bootstrap.weights) {   
    if (!is.null(w2t)) {
      return(PPLCscore2(dNit = dNit, Yit = Yit, Xit = Xit, beta = b,
                        w1_t = w1t, w2_t = w2t, rw = weights))
    } else if (!is.null(w1t)) {
      return(PPLCscore1(dNit = dNit, Yit = Yit, Xit = Xit, beta = b,
                        w_t = w1t, rw = weights))
    } else {
      return(PPLscore(dNit, Yit, Xit, beta = b, rw = weights))
    }
  }
  est <- EL_saddle(sfun, initial, maxit1, maxit2, tol, trace)
  bhat <- drop(est$beta)
  S <- sfun(bhat)
  Sigmahat <- sq_empirical(S$score.subject, bootstrap.weights = 1) * n
  Hess <- S$diff.score
  Vhat <- Ginv_eigen(t(Hess) %*% Ginv_eigen(Sigmahat) %*% Hess) / n
  rownames(Vhat) <- colnames(Vhat) <- dimnames(Xit)[[2]]
  names(bhat) <- dimnames(Xit)[[2]]
  if (sum(is.na(bhat))>0) {
    dLhat <- NA
  } else {
    S0t <- colSums(exp(apply(aperm(Xit, c(2, 1, 3)) * bhat, c(2, 3), sum)) * 
                   as.vector(Yit) * bootstrap.weights)
    dLhat <- colSums(dNit * bootstrap.weights) * (S0t!=0)/(S0t+(S0t==0))
  }
  list(coefficient = bhat,
       covariance = Vhat,
       diff.baseline = dLhat)
}

CPPL.GMM <- function(dNit, Yit, Xit, w1t, w2t,
                     initial = NULL, bootstrap.weights = NULL, 
                     maxit1 = 100, maxit2 = 10, tol = 1e-7, trace = FALSE) {
  n <- nrow(dNit)
  p <- dim(Xit)[2]
  if (is.null(initial) | initial == 0) initial <- rep(0, p)
  if (is.null(bootstrap.weights)) bootstrap.weights <- rep(1 / n, n)
  sfun0 <- function(b) PPLscore(dNit, Yit, Xit, beta = b, rw = bootstrap.weights)
  est0 <- GMM_NRline(sfun0, initial, diag(p), maxit1, tol, trace)
  GMM_initial <- drop(est0$beta)
  sfun <- function(b, weights = bootstrap.weights) {
    if (!is.null(w2t)) {
      return(PPLCscore2(dNit = dNit, Yit = Yit, Xit = Xit, beta = b,
                        w1_t = w1t, w2_t = w2t, rw = weights))
    } else if (!is.null(w1t)) {
      return(PPLCscore1(dNit = dNit, Yit = Yit, Xit = Xit, beta = b,
                        w_t = w1t, rw = weights))
    } else {
      return(PPLscore(dNit, Yit, Xit, beta = b, rw = weights))
    }
  }
  S <- sfun(GMM_initial)
  What <- Ginv_eigen(sq_empirical(S$score.subject, bootstrap.weights = 1) * n)
  est <- GMM_NRline(sfun, GMM_initial, What, maxit1, tol, trace)
  bhat <- drop(est$beta)
  names(bhat) <- dimnames(Xit)[[2]]
  S <- sfun(bhat)
  Hess <- S$diff.score
  Sigmahat <- (sq_empirical(S$score.subject, bootstrap.weights = 1) * n)
  Vhat <- Ginv_eigen(t(Hess) %*% Ginv_eigen(Sigmahat) %*% Hess) / n
  rownames(Vhat) <- colnames(Vhat) <- dimnames(Xit)[[2]]
  names(bhat) <- dimnames(Xit)[[2]]
  if (sum(is.na(bhat))>0) {
    dLhat <- NA
  } else {
    S0t <- colSums(exp(apply(aperm(Xit, c(2, 1, 3)) * bhat, c(2, 3), sum)) * 
                   as.vector(Yit) * bootstrap.weights)
    dLhat <- colSums(dNit * bootstrap.weights) * (S0t!=0)/(S0t+(S0t==0))
  }
  list(coefficient = bhat,
       covariance = Vhat,
       diff.baseline = dLhat)  
}

#' Other backgroud functions; I haven't touched those functions.
#' @noRd

PPLscore <- function(dNit, Yit, Xit, beta, rw = NULL) {
  number_n <- dim(dNit)[1]
  number_t <- dim(dNit)[2]
  number_p <- dim(Xit)[2]
  if (is.vector(rw) == FALSE) {
    rw <- rep(1/number_n, times = number_n)
  }
  Xsqit <- array(Xit[,rep(1:number_p, times = number_p),] * 
                 Xit[,rep(1:number_p, each = number_p),],
                 c(number_n, number_p^2, number_t))
  IXit <- apply(aperm(Xit, c(2, 1, 3)) * as.vector(beta), c(2, 3), sum) 
  IXit <- pmin(IXit, 700)
  YeIXit <- Yit * exp(IXit)
  S0 <- colSums(YeIXit * rw)
  S1 <- apply(aperm(Xit, c(1, 3, 2)) * as.vector(YeIXit) * rw, c(2, 3), sum)
  S2 <- apply(aperm(Xsqit, c(1, 3, 2)) * as.vector(YeIXit) * rw, c(2, 3), sum)
  S1std <- S1 * (S0!=0)/(S0+(S0==0))                #c(number_t,number_p)
  S2std <- S2 * (S0!=0)/(S0+(S0==0))                #c(number_t,number_p * number_p)
  dLhat <- colSums(dNit * rw) * (S0!=0)/(S0+(S0==0))  #c(number_t)
  dAit <- t(YeIXit) * dLhat                         #c(number_t,number_n)
  dMit <- dNit-t(dAit)                            #c(number_n,number_t)
  DXit <- aperm(Xit, c(3, 2, 1))-as.vector(S1std) #c(number_t,number_p,number_n)
  DXsqit <- array(DXit[ , rep(1:number_p, times = number_p),] *
                  DXit[ , rep(1:number_p, each = number_p), ],
                  c(number_t, number_p * number_p, number_n))
  S1stdsq <- array(S1std[ , rep(1:number_p, times = number_p)] *
                   S1std[ , rep(1:number_p, each = number_p)],
                   c(number_t, number_p * number_p))
  scorei <- apply(aperm(DXit, c(3, 1, 2)) * as.vector(dMit), c(1, 3), sum) * rw #c(number_n,number_p)
  informationi <- array(-apply(aperm(DXsqit, c(1, 3, 2)) * as.vector(dAit), c(2, 3), sum) +
                        rowSums(dMit[rep(1:number_n, times = number_p^2) , ] * 
                                t(S1stdsq-S2std)[rep(1:(number_p^2), each = number_n),]),
                        c(number_n, number_p, number_p)) * rw
  score <- colSums(scorei)
  information <- apply(informationi, c(2, 3), sum)
  results <- list(score.subject = scorei,
                  diff.score.subject = informationi,
                  score = score,
                  diff.score = information)

  return(results)
}

PPLCscore1 <- function(dNit, Yit, Xit, beta, w_t, rw = NULL) {
  number_n <- dim(dNit)[1]
  number_t <- dim(dNit)[2]
  number_p <- dim(Xit)[2]
  if (is.vector(rw)==FALSE) {
    rw <- rep(1/number_n, times = number_n)
  }
  Xsqit <- array(Xit[ , rep(1:number_p, times = number_p), ] *
                 Xit[ , rep(1:number_p, each = number_p), ],
                 c(number_n, number_p^2, number_t))
  IXit <- apply(aperm(Xit, c(2, 1, 3)) * as.vector(beta), c(2, 3), sum) #c(number_n,number_t)
  IXit <- pmin(IXit, 700)
  YeIXit <- Yit * exp(IXit)                                             #c(number_n,number_t)
  S0 <- colSums(YeIXit * rw)
  S1 <- apply(aperm(Xit, c(1, 3, 2)) * as.vector(YeIXit) * rw, c(2, 3), sum)
  S2 <- apply(aperm(Xsqit, c(1, 3, 2)) * as.vector(YeIXit) * rw, c(2, 3), sum)
  S1std <- S1 * (S0!=0)/(S0+(S0==0))                #c(number_t,number_p)
  S2std <- S2 * (S0!=0)/(S0+(S0==0))                #c(number_t,number_p * number_p)
  dLhat <- colSums(dNit * rw) * (S0!=0)/(S0+(S0==0))  #c(number_t)
  dAit <- t(YeIXit) * dLhat                         #c(number_t,number_n)
  tdAit <- dAit * w_t                               #c(number_t,number_n)
  dMit <- dNit-t(dAit)                            #c(number_n,number_t)
  tdMit <- t(t(dMit) * w_t)                         #c(number_n,number_t)
  DXit <- aperm(Xit, c(3, 2, 1))-as.vector(S1std) #c(number_t,number_p,number_n)
  DXsqit <- array(DXit[ , rep(1:number_p, times = number_p), ] *
                  DXit[ , rep(1:number_p, each = number_p), ],
                  c(number_t, number_p * number_p, number_n))
  S1stdsq <- array(S1std[ , rep(1:number_p, times = number_p)] *
                   S1std[ , rep(1:number_p, each = number_p)],
                   c(number_t, number_p * number_p))
  scorePi <- apply(aperm(DXit, c(3, 1, 2)) * as.vector(dMit), c(1, 3), sum) * rw  #c(number_n,number_p)
  score1i <- apply(aperm(DXit, c(3, 1, 2)) * as.vector(tdMit), c(1, 3), sum) * rw #c(number_n,number_p)
  scorei <- cbind(scorePi, score1i)
  informationPi <- array(-apply(aperm(DXsqit, c(1, 3, 2)) * as.vector(dAit), c(2, 3),sum)+
                         rowSums(dMit[rep(1:number_n, times = number_p^2),] * 
                                 t(S1stdsq-S2std)[rep(1:(number_p^2), each = number_n), ]),
                         c(number_n, number_p^2)) * rw
  information1i <- array(-apply(aperm(DXsqit, c(1, 3, 2)) * as.vector(tdAit), c(2, 3),sum)+
                         rowSums(tdMit[rep(1:number_n, times = number_p^2),] * 
                                 t(S1stdsq-S2std)[rep(1:(number_p^2), each = number_n), ]),
                         c(number_n, number_p^2)) * rw
  informationi <-
    cbind(informationPi,
          information1i)[ ,as.vector(aperm(array(1:(2 * number_p^2),
                                                 c(number_p, number_p, 2)), c(1, 3, 2)))]
  dim(informationi) <- c(number_n, 2 * number_p, number_p)
  score <- apply(scorei, 2, sum)
  information <- apply(informationi, c(2, 3), sum)
  results <- list(score.subject = scorei,
                  diff.score.subject = informationi,
                  score = score,
                  diff.score = information)

  return(results)
}

PPLCscore2 <- function(dNit, Yit, Xit, beta, w1_t, w2_t, rw = NULL) {
  number_n <- dim(dNit)[1]
  number_t <- dim(dNit)[2]
  number_p <- dim(Xit)[2]
  if (is.vector(rw)==FALSE) {
    rw <- rep(1/number_n, times = number_n)
  }
  Xsqit <- array(Xit[ , rep(1:number_p, times = number_p),] * 
                 Xit[ , rep(1:number_p, each = number_p),],
                 c(number_n, number_p^2, number_t))
  IXit <- apply(aperm(Xit, c(2, 1, 3)) * as.vector(beta), c(2, 3), sum) #c(number_n,number_t)
  IXit <- pmin(IXit, 700)
  YeIXit <- Yit * exp(IXit)                                             #c(number_n,number_t)
  S0 <- colSums(YeIXit * rw)
  S1 <- apply(aperm(Xit, c(1, 3, 2)) * as.vector(YeIXit) * rw, c(2, 3), sum)
  S2 <- apply(aperm(Xsqit, c(1, 3, 2)) * as.vector(YeIXit) * rw, c(2, 3), sum)
  S1std <- S1 * (S0!=0)/(S0+(S0==0))                #c(number_t,number_p)
  S2std <- S2 * (S0!=0)/(S0+(S0==0))                #c(number_t,number_p * number_p)
  dLhat <- colSums(dNit * rw) * (S0!=0)/(S0+(S0==0))  #c(number_t)
  dAit <- t(YeIXit) * dLhat                         #c(number_t,number_n)
  t1dAit <- dAit * w1_t                             #c(number_t,number_n)
  t2dAit <- dAit * w2_t                             #c(number_t,number_n)
  dMit <- dNit-t(dAit)                            #c(number_n,number_t)
  t1dMit <- t(t(dMit) * w1_t)                       #c(number_n,number_t)
  t2dMit <- t(t(dMit) * w2_t)                       #c(number_n,number_t)
  DXit <- aperm(Xit, c(3, 2, 1))-as.vector(S1std) #c(number_t,number_p,number_n)
  DXsqit <- array(DXit[ , rep(1:number_p, times = number_p),] * 
                  DXit[,rep(1:number_p, each = number_p),],
                  c(number_t, number_p * number_p, number_n))
  S1stdsq <- array(S1std[ , rep(1:number_p, times = number_p)] *
                   S1std[ , rep(1:number_p, each = number_p)],
                   c(number_t, number_p * number_p))
  scorePi <- apply(aperm(DXit, c(3, 1, 2)) * as.vector(dMit), c(1, 3), sum) * rw
  score1i <- apply(aperm(DXit, c(3, 1, 2)) * as.vector(t1dMit), c(1, 3), sum) * rw
  score2i <- apply(aperm(DXit, c(3, 1, 2)) * as.vector(t2dMit), c(1, 3), sum) * rw
  scorei <- cbind(scorePi, score1i, score2i)
  informationPi <- array(-apply(aperm(DXsqit, c(1, 3, 2)) * as.vector(dAit), c(2, 3), sum)+
                         rowSums(dMit[rep(1:number_n, times = number_p^2),] * 
                                 t(S1stdsq-S2std)[rep(1:(number_p^2), each = number_n), ]),
                         c(number_n, number_p^2)) * rw
  information1i <- array(-apply(aperm(DXsqit, c(1, 3, 2)) * as.vector(t1dAit), c(2, 3), sum)+
                         rowSums(t1dMit[rep(1:number_n, times = number_p^2),] * 
                                 t(S1stdsq-S2std)[rep(1:(number_p^2), each = number_n), ]),
                         c(number_n, number_p^2)) * rw
  information2i <- array(-apply(aperm(DXsqit, c(1, 3, 2)) * as.vector(t2dAit), c(2, 3), sum)+
                         rowSums(t2dMit[rep(1:number_n, times = number_p^2),] * 
                                 t(S1stdsq-S2std)[rep(1:(number_p^2), each = number_n), ]),
                         c(number_n, number_p^2)) * rw
  informationi <-
    cbind(informationPi,
          information1i,
          information2i)[ , as.vector(aperm(array(1:(3 * number_p^2),
                                                  c(number_p, number_p, 3)), c(1, 3, 2)))]
  dim(informationi) <- c(number_n, 3 * number_p, number_p)
  score <- apply(scorei, 2, sum)
  information <- apply(informationi, c(2, 3), sum)
  results <- list(score.subject = scorei,
                  diff.score.subject = informationi,
                  score = score,
                  diff.score = information)
  return(results)
}

## EL background functions

inner_loop <- function(score.subject,
                       iteration.max = 100, step.max = 10, tolerance = 1e-7) {
  number_n <- dim(score.subject)[1]
  number_m <- dim(score.subject)[2]
  eta_step <- rep(0, times = number_m)
  denominator_step <- as.vector(1-score.subject %*% eta_step)
  L_step <- -sum(log(denominator_step))
  score_eta_step <- colSums(score.subject/denominator_step)
  information_eta_step <- matrix(colSums(
    as.matrix(score.subject[ , rep(1:number_m, times = number_m)] * 
              score.subject[ , rep(1:number_m , each = number_m)])/(denominator_step^2)),
    nrow = number_m, ncol = number_m)
  for (iter_inner in 1:iteration.max) {
    step_size <- 1
    direction_step <- Ginv_eigen(information_eta_step) %*% score_eta_step
    eta_step_new <- eta_step-direction_step * step_size
    for (iter_step in 1:step.max) {
      denominator_step_new <- as.vector(1-score.subject %*% eta_step_new)
      if (sum(denominator_step_new<=0)>=1) {
        L_step_new <- L_step
      } else {
        L_step_new <- -sum(log(denominator_step_new))
      }
      if ((L_step_new>=L_step)+sum(denominator_step_new<=1/number_n)>=1) {
        step_size <- step_size/2
        eta_step_new <- eta_step-direction_step * step_size
      }else
        break
    }
    if (L_step_new<L_step-tolerance) {
      eta_step <- eta_step_new
      denominator_step <- denominator_step_new
      L_step <- L_step_new
      score_eta_step <- colSums(score.subject/denominator_step)
      information_eta_step <- matrix(colSums(
        as.matrix(score.subject[ , rep(1:number_m, times = number_m)] * 
                  score.subject[ , rep(1:number_m, each = number_m)])/(denominator_step^2)),
        nrow = number_m, ncol = number_m)
    } else
      break
  }
  results <- list(eta = eta_step,
                  denominator = denominator_step,
                  L = L_step,
                  score = score_eta_step,
                  information = information_eta_step)
  return(results)
}

EL_saddle <- function(score.subject.function, initial,
                      iteration.max = 100, step.max = 10,
                      tolerance = 1e-7, trace = FALSE) {
  number_p <- length(initial)
  beta_step <- initial
  scoref_step <- score.subject.function(beta_step)
  scorei_step <- scoref_step$score.subject
  informationi_step <- scoref_step$diff.score.subject
  inner <- inner_loop(score.subject = scorei_step,
                      iteration.max = iteration.max, step.max = step.max,
                      tolerance = tolerance)
  eta_step <- inner$eta
  denominator_step <- inner$denominator
  L_step <- inner$L
  information_eta_step <- inner$information
  T_step <- apply(informationi_step/denominator_step, c(2, 3), sum)
  if (trace) {
    cat(paste("Iteration", 0, " EL=", L_step, sep = ""))
    cat("\n")
  }
  for (iter_outer in 1:iteration.max) {
    step_size <- 1
    direction_step <- Ginv_eigen(-t(T_step) %*% Ginv_eigen(information_eta_step) %*% T_step) %*%
      t(T_step) %*% eta_step
    beta_step_new <- beta_step-direction_step * step_size
    for (iter_step in 1:step.max) {
      scoref_step_new <- score.subject.function(beta_step_new)
      scorei_step_new <- scoref_step_new$score.subject
      informationi_step_new <- scoref_step_new$diff.score.subject
      inner <- inner_loop(score.subject = scorei_step_new,
                          iteration.max = iteration.max, step.max = step.max,
                          tolerance = tolerance)
      eta_step_new <- inner$eta
      denominator_step_new <- inner$denominator
      L_step_new <- inner$L
      information_eta_step_new <- inner$information
      if (L_step_new<=L_step) {
        step_size <- step_size/2
        beta_step_new <- beta_step-direction_step * step_size
      } else
        break
    }
    if (L_step_new>L_step+tolerance) {
      beta_step <- beta_step_new
      scorei_step <- scorei_step_new
      informationi_step <- informationi_step_new
      eta_step <- eta_step_new
      denominator_step <- denominator_step_new
      L_step <- L_step_new
      information_eta_step <- information_eta_step_new
      T_step <- apply(informationi_step/denominator_step, c(2, 3), sum)
      if (trace) {
        cat(paste("Iteration", iter_outer, " EL=", L_step, sep = ""))
        cat("\n")
      }
    } else
      break
  }
  results <- list(beta = beta_step,
                  eta = eta_step,
                  L = L_step)
  return(results)
}

## GMM functions

GMM_NRline <- function(score.function, initial, weight.matrix,
                       iteration.max = 100, tolerance = 1e-7,
                       trace = FALSE) {
  f0 <- function(theta) {
    scoref <- score.function(theta)
    SS <- sum(colSums(weight.matrix * scoref$score) * scoref$score)
    return(SS)
  }
  f2 <- function(theta) {
    scoref <- score.function(theta)
    SS <- sum(colSums(weight.matrix * scoref$score) * scoref$score)
    gradient <- t(scoref$diff.score) %*% weight.matrix %*% scoref$score * 2
    Hessian <- t(scoref$diff.score) %*% weight.matrix %*% scoref$diff.score * 2
    results  <- list(value = SS,
                     gradient = gradient,
                     Hessian = Hessian)
    return(results)
  }
  estimation <- NRline(f.diff2 = f2, f.diff0 = f0, initial = initial,
                       iteration.max = iteration.max, tolerance = tolerance,
                       trace = trace)
  results <- list(beta = estimation$minimizer,
                  score = estimation$gradient,
                  information = estimation$Hessian,
                  value = estimation$value)
  return(results)
}

#' @importFrom stats nlminb
NRline <- function(f.diff2, f.diff0 = NULL, initial,
                   iteration.max = 100, tolerance = 1e-7,
                   trace = FALSE) {
  if (!is.function(f.diff0)) {
    f.diff0 <- function(theta) {
      value <- f.diff2(theta)$value
      return(value)
    }
  }
  theta_step <- initial
  evaluation_step <- f.diff2(theta_step)
  obj_step <- evaluation_step$value
  grad_step <- evaluation_step$gradient
  Hess_step <- evaluation_step$Hessian
  Newton_step <- Ginv_eigen(Hess_step) %*% grad_step
  if (trace) {
    cat(paste("s", 0, " value=", obj_step, sep = ""))
    cat("\n")
  }
  for (s in 1:iteration.max) {
    obj_line <- function(alpha) {
      theta_step_new <- theta_step-Newton_step * alpha
      objective <- f.diff0(theta_step_new)
      return(objective)
    }
    line_search <- nlminb(start = 0, objective = obj_line)
    if (line_search$objective<obj_step-tolerance) {
      theta_step <- theta_step-Newton_step * line_search$par
      evaluation_step <- f.diff2(theta_step)
      obj_step <- evaluation_step$value
      grad_step <- evaluation_step$gradient
      Hess_step <- evaluation_step$Hessian
      Newton_step <- Ginv_eigen(Hess_step) %*% grad_step
      if (trace) {
        cat(paste("s", s, " value=", obj_step, sep = ""))
        cat("\n")
      }
    } else
      break
  }
  results <- list(minimizer = theta_step,
                  value = obj_step,
                  gradient = grad_step,
                  Hessian = Hess_step)
  return(results)
}

Gauss_row <- function(A) {
  M <- A
  m <- nrow(A)
  n <- ncol(A)
  k_max <- min(m, n)
  npivot <- 0
  for (k in 1:n) {
    i_max <- which.max(abs(M[(npivot+1):m, k]))[1]+npivot
    if (M[i_max,k]==0)
      next
    npivot <- npivot+1
    M[c(npivot,i_max), ] <- M[c(i_max,npivot), ]
    for (i in 1:m) {
      if (i==npivot)
        next
      M[i, ] <- M[i, ]-M[npivot, ] * M[i, k]/M[npivot, k]
    }
    M[npivot, ] <- M[npivot, ]/M[npivot, k]
    if (npivot==k_max)
      break
  }
  return(M)
}

Ginv_eigen <- function(S, tolerance = 1e-7) {
  S <- as.matrix(S)
  decomp <- eigen(S)
  inv.values <- (abs(decomp$values)>tolerance)/(decomp$values+(abs(decomp$values)<=tolerance))
  D.inv <- diag(dim(S)[1])
  diag(D.inv) <- inv.values
  S.ginv <- decomp$vectors %*% D.inv %*% solve(decomp$vectors)
  return(S.ginv)
}

sq_empirical <- function(X, bootstrap.weights = NULL) {
  X <- as.matrix(X)
  number_n <- dim(X)[1]
  number_p <- dim(X)[2]
  if (is.vector(bootstrap.weights)==FALSE) {
    bootstrap.weights <- rep(1/number_n, times = number_n)
  }
  S <- matrix(colSums(matrix(X[,rep(1:number_p, times = number_p)] * 
                             X[,rep(1:number_p, each = number_p)],
                             nrow = number_n, ncol = number_p^2) * bootstrap.weights),
              nrow = number_p, ncol = number_p)
  return(S)
}
