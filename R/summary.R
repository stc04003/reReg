#' @export
print.reReg <- function(x, ...) {
    if (!is.reReg(x))
        stop("Must be a reReg object")
    cat("Call: \n")
    dput(x$call)
    if (all(!is.na(x$alpha)) & !is.null(x$alpha)) {
        if(x$recType == "cox.LWYY")
            cat("\nFitted with the Cox model of Lin et al. (2000):")
        if(x$recType == "cox.GL")
            cat("\nFitted with the Cox model of Ghosh and Lin (2002):")
        if(x$recType == "am.GL")
            cat("\nFitted with the accelerated mean model of Ghosh and Lin (2003):")
        cat("\nRecurrent event process:")
        if (x$recType == "sc") {
            p <- length(x$alpha) / 2
            mat <- rbind(c("Shape", rep("", p - 1), "Size", rep("", p - 1)),
                         rep(x$varNames, 2), format(x$alpha, digits = 5))
            mat <- cbind(mat[,1:p], "    ", mat[,1:p])
            prmatrix(mat, rowlab = rep("", nrow(mat)), collab = rep("", 1 + ncol(mat)), quote = FALSE)
        } else {
            mat <- rbind(x$varNames, format(x$alpha, digits = 5))
            prmatrix(mat, rowlab = rep("", nrow(mat)), collab = rep("", ncol(mat)), quote = FALSE)
        }
        if (length(x$beta) > 0) {
            cat("\nTerminal event:")
            if (x$recType == "sc") {
                p <- length(x$beta) / 2
                mat <- rbind(c("Shape", rep("", p - 1), "Size", rep("", p - 1)),
                             rep(x$varNames, 2), format(x$beta, digits = 5))
                mat <- cbind(mat[,1:p], "    ", mat[,1:p])
                prmatrix(mat, rowlab = rep("", nrow(mat)), collab = rep("", 1 + ncol(mat)),
                         quote = FALSE)
            } else {
                mat <- rbind(x$varNames, format(x$beta, digits = 5))
                prmatrix(mat, rowlab = rep("", nrow(mat)), collab = rep("", ncol(mat)), quote = FALSE)
            }
        }
        ## print.default(format(x$alpha, digits = digits), print.gap = 2L, quote = FALSE)
    } else {
        n <- length(unique(x$DF$id))
        nevent <- sum(x$DF$event)
        avg.event <- sum(x$DF$event) / n
        cat("\nNumber of subjects:", n)
        cat("\nNumber of recurrent events:", nevent)
        cat("\nAverage recurrent events per subject:", avg.event)
        cat("\n")
    }
}

pvalTab <- function(pe, se) {
    if (is.null(se)) se <- NA
    cbind(Estimate = round(pe, 3), StdErr = round(se, 3),
          z.value = round(pe / se, 3), p.value = round(2 * pnorm(-abs(pe / se)), 3))    
}
    
#' @export
summary.reReg <- function(object, test = FALSE, ...) {
    if (!is.reReg(object)) stop("Must be a reReg object")
    if (object$method == "nonparametric") {
        t0 <- sort(unique(c(object$DF$time1, object$DF$time2)))
        t0 <- t0[t0 > 0]
        out <- list(call = object$call, method = object$method, 
                    tabA = data.frame(time = t0, rate = object$Lam0(t0), hazard = object$Haz0(t0)))
        out
    }
    if (object$method != "nonparametric" & all(!is.na(object$alpha))) {
        tabA <- pvalTab(object$alpha, object$alphaSE)
        if (object$recType == "sc") {
            p <- length(object$alpha) / 2
            rownames(tabA) <- rep(object$varNames, 2)
            tabA <- list(tabA1 = tabA[1:p,, drop = FALSE], tabA2 = tabA[-(1:p),, drop = FALSE])
        } else rownames(tabA) <- object$varNames
        out <- list(call = object$call, method = object$method, tabA = tabA)
        if (!is.null(object$beta)) {
            tabB <- pvalTab(object$beta, object$betaSE)
            if (object$temType == "sc") {
                p <- length(object$beta) / 2
                rownames(tabB) <- rep(object$varNames, 2)
                tabB <- list(tabB1 = tabB[1:p,, drop = FALSE], tabB2 = tabB[-(1:p),, drop = FALSE])
            } else rownames(tabB) <- object$varNames
            out$tabB <- tabB
        }
        if (object$recType == "sc" & object$se == "resampling") {
            p <- length(object$alpha) / 2
            out$HA.chi <- object$alpha[1:p] %*%
                solve(object$varMat[1:p, 1:p, drop = FALSE]) %*% object$alpha[1:p]
            out$HB.chi <- object$alpha[-(1:p)] %*%
                solve(object$varMat[1:p, 1:p, drop = FALSE] +
                      object$varMat[(p+2):(2*p+1), (p+2):(2*p+1), drop = FALSE] +
                      2 * object$varMat[1:p, (p+2):(2*p+1), drop = FALSE]) %*%
                object$alpha[-(1:p)]
            g <- object$alpha[-(1:p)] - object$alpha[1:p]
            out$HG.chi <- g %*% solve(object$varMat[(p+2):(2*p+1), (p+2):(2*p+1), drop = FALSE]) %*% g
            out$HA.pval <- 1 - pchisq(out$HA.chi, p)
            out$HB.pval <- 1 - pchisq(out$HB.chi, p)
            out$HG.pval <- 1 - pchisq(out$HG.chi, p)
        }
        out$recType <- object$recType
        out$temType <- object$temType
        out$test <- test
    }
    class(out) <- "summary.reReg"
    out
}

printCoefmat2 <- function(tab) 
    printCoefmat(as.data.frame(tab), P.values = TRUE, has.Pvalue = TRUE, signif.legend = FALSE)

#' @export
print.summary.reReg <- function(x, ...) {
    cat("Call: \n")
    dput(x$call)
    if (x$method != "nonparametric" & !is.na(x$tabA)[1]) {
        if(x$recType == "cox.LWYY")
            cat("\nFitted with the Cox model of Lin et al. (2000):")
        if(x$recType == "cox.GL")
            cat("\nFitted with the Cox model of Ghosh and Lin (2002):")
        if(x$recType == "am.GL")
            cat("\nFitted with the accelerated mean model of Ghosh and Lin (2003):")
        if (x$recType == "sc") {
            p <- nrow(x$tabA$tabA1)
            cat("\nRecurrent event process (shape):\n")
            printCoefmat2(x$tabA[[1]])
            cat("\nRecurrent event process (size):\n")
            printCoefmat2(x$tabA[[2]])
            if (x$test) {
                cat("\nHypothesis tests:")
                cat("\nHo: shape = 0 (Cox-type model):")
                cat(paste("\n     X-squared = ", round(x$HA.chi, 4), ", df = ", p,
                          ", p-value = ", round(x$HA.pval, 4), sep = ""))
                cat("\nHo: shape = size (Accelerated rate model):")
                cat(paste("\n     X-squared = ", round(x$HB.chi, 4), ", df = ", p,
                          ", p-value = ", round(x$HB.pval, 4), sep = ""))
                cat("\nHo: size = 0 (Accelerated mean model):")
                cat(paste("\n     X-squared = ", round(x$HG.chi, 4), ", df = ", p,
                          ", p-value = ", round(x$HG.pval, 4), sep = ""))
            }
        } else {
            cat("\nRecurrent event process:\n")
            printCoefmat(x$tabA)
        }
        ## Lin-Wei-Yang-Ying method (fitted with coxph with robust variance)
        if (x$temType != ".") {
            if (x$temType == "sc") {
                p <- nrow(x$tabB$tabB1)
                cat("\nTerminal event (shape):\n")
                printCoefmat2(x$tabB[[1]])
                cat("\nTerminal event (size):\n")
                printCoefmat2(x$tabB[[2]])
            } else {
                  cat("\nTerminal event:\n")
                  printCoefmat(x$tabB)
            }
        }
    }    
    if (x$method == "nonparametric") {
        cat("\n")
        print(round(unique(x$tabA), 4), row.names = FALSE)
    }
    cat("\n")
}

#' @export
coef.reReg <- function(object, ...) {
    if (is.null(object$beta)) return(object$alpha)
    return(as.numeric(c(object$alpha, object$beta)))
}

#' @export
vcov.reReg <- function(object, ...) {
    if (is.null(object$betaVar)) return(object$alphaVar)
    return(list(alpha.vcov = object$alphaVar, beta.vcov = object$betaVar))
}
