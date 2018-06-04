print.reSurv <- function(x, ...) {
   print(x$reTb)
}

print.reReg <- function(x, ...) {
    if (!is.reReg(x))
        stop("Must be a reReg x")
    cat("Call:\n")
    print(x$call)
    if (all(!is.na(x$alpha))) {
        if (x$method != "cox.LWYY") {
            if (x$method == "am.XC") 
                cat("\nMethod: Joint Scale-Change Model \n")
            if (x$method == "am.GL")
                cat("\nMethod: Ghosh-Lin Model \n")
            if (x$method == "cox.HW")
                cat("\nMethod: Huang-Wang Model \n")
            cat("\nCoefficients:\n")
            coefMat <- rbind(x$alpha, x$beta)
            rownames(coefMat) <- c("alpha", "beta")
            colnames(coefMat) <- x$varNames
        } else {
            cat("\nMethod: Lin-Wei-Yang-Ying Model \n")
            cat("\nCoefficients:\n")
            coefMat <- rbind(x$alpha, NULL)
            rownames(coefMat) <- c("alpha")
            colnames(coefMat) <- x$varNames        
        }
        print(coefMat)
    } else {
        n <- length(unique(x$DF$id))
        nevent <- sum(x$DF$event)
        avg.event <- sum(x$DF$event) / n
        tab <- data.frame(n = n, event = nevent, avg.event = avg.event)
        colnames(tab) <- c("n", "number of recurrent events", "average recurrent events per subject")
        cat("\n")
        print(tab, row.names = FALSE)
        cat("\n")            
        ## cat(capture.output(tab), sep = "\n")
    }
}

summary.reReg <- function(object, ...) {
    if (!is.reReg(object)) stop("Must be a reReg x")
    if (all(!is.na(object$alpha))) {
        if (object$se == "NULL") object$alphaSE <- object$betaSE <- NA
        tabA <- cbind(Estimate = round(object$alpha, 3),
                      StdErr = round(object$alphaSE, 3),
                      z.value = round(object$alpha / object$alphaSE, 3),
                      p.value = round(2 * pnorm(-abs(object$alpha / object$alphaSE)), 3))
        rownames(tabA) <- object$varNames
        tabB <- cbind(Estimate = round(object$beta, 3),
                      StdErr = round(object$betaSE, 3),
                      z.value = round(object$beta / object$betaSE, 3),
                      p.value = round(2 * pnorm(-abs(object$beta / object$betaSE)), 3))
        rownames(tabB) <- object$varNames
        out <- list(call = object$call, method = object$method, tabA = tabA, tabB = tabB)
    } else {
        out <- list(call = object$call, method = object$method, 
                    tabB = data.frame(Time = object$t0, cum.rate = object$lam,
                                      cum.rate.L = object$lamL, cum.rate.U = object$lamU,
                                      cum.haz = object$haz, cum.haz.L = object$hazL,
                                      cum.haz.U = object$hazU), tabA = NA)
        ## assuming tabA is na if fit with ~1
    }
    if (object$method == "sc.XCYH" & object$se == "resampling") {
        p <- length(object$alpha)
        out$HA.chi <- object$alpha %*% solve(object$varMat[1:p, 1:p]) %*% object$alpha
        out$HB.chi <- object$beta %*%
            solve(object$varMat[1:p, 1:p] + object$varMat[(p+2):(2*p+1), (p+2):(2*p+1)] +
                  2 * object$varMat[1:p, (p+2):(2*p+1)]) %*% object$beta
        out$HG.chi <- object$gamma[-1] %*%
            solve(object$varMat[(p+2):(2*p+1), (p+2):(2*p+1)]) %*% object$gamma[-1]
        out$HA.pval <- 1 - pchisq(out$HA.chi, length(object$alpha))
        out$HB.pval <- 1 - pchisq(out$HB.chi, length(object$beta))
        out$HG.pval <- 1 - pchisq(out$HG.chi, length(object$beta))
    }
    class(out) <- "summary.reReg"
    out
}

print.summary.reReg <- function(x, ...) {
    cat("Call: ")
    print(x$call)
    if (!is.na(x$tabA)[1]) {
        if (x$method != "sc.XCYH" & x$method != "cox.LWYY") {
            if (x$method == "am.XC") 
                cat("\nMethod: Joint Scale-Change Model \n")
            if (x$method == "am.GL")
                cat("\nMethod: Ghosh-Lin Model \n")
            if (x$method == "am.XCHWY")
                cat("\nMethod: Xu et al. (2016) Model \n")
            if (x$method == "cox.HW")
                cat("\nMethod: Huang-Wang Model \n")
            cat("\nCoefficients (rate):\n")
            printCoefmat(as.data.frame(x$tabA), P.values = TRUE, has.Pvalue = TRUE)
            cat("\nCoefficients (hazard):\n")
            printCoefmat(as.data.frame(x$tabB), P.values = TRUE, has.Pvalue = TRUE)
        }
        if (x$method == "sc.XCYH") {
            p <- nrow(x$tabA)
            cat("\nMethod: Generalized Scale-Change Model \n")
            cat("\nScale-change effect:\n")
            printCoefmat(as.data.frame(x$tabA), P.values = TRUE, has.Pvalue = TRUE)
            cat("\nMultiplicative coefficients:\n")
            printCoefmat(as.data.frame(x$tabB), P.values = TRUE, has.Pvalue = TRUE)
            if (!is.null(x$HA.chi)) {
                cat("\nHypothesis tests:\n")
                cat("\nH0 Cox-type model:")
                cat(paste("\n     X-squared = ", round(x$HA.chi, 4), ", df = ", p,
                          ", p-value = ", round(x$HA.pval, 4), sep = ""))
                cat("\nH0 Accelerated rate model:")
                cat(paste("\n     X-squared = ", round(x$HB.chi, 4), ", df = ", p,
                          ", p-value = ", round(x$HB.pval, 4), sep = ""))
                cat("\nH0 Accelerated mean model:")
                cat(paste("\n     X-squared = ", round(x$HG.chi, 4), ", df = ", p,
                          ", p-value = ", round(x$HG.pval, 4), sep = ""))
            }
        }
        if (x$method == "cox.LWYY") {
            cat("\nMethod: Lin-Wei-Yang-Ying method \n")
            cat("\nCoefficients effect:\n")
            printCoefmat(as.data.frame(x$tabA), P.values = TRUE, has.Pvalue = TRUE)
        }
    }
    if (is.na(x$tabA)[1]) {
        colnames(x$tabB) <- c("time", "cum.rate", "95% LCL", "95% UCL",
                                   "cum.hazard", "95% LCL", "95% UCL")
        x$tabB[,sapply(1:ncol(x$tabB), function(y) all(is.na(x$tabB[,y])))] <- NULL
        cat("\n")
        print(round(unique(x$tabB), 4), row.names = FALSE)
    }
    cat("\n")
}

coef.reReg <- function(object, ...) {
    as.numeric(c(object$alpha, object$beta))
}

