#' @exportS3Method print reReg
print.reReg <- function(x, ...) {
    if (!is.reReg(x))
        stop("Must be a reReg object")
    cat("Call: \n")
    dput(x$call)
    if(x$typeRec == "cox.LWYY")
        cat("\nFitted with the Cox model of Lin et al. (2000):")
    if(x$typeRec == "cox.GL")
        cat("\nFitted with the Cox model of Ghosh and Lin (2002):")
    if(x$typeRec == "am.GL")
        cat("\nFitted with the accelerated mean model of Ghosh and Lin (2003):")
    if (x$typeTem != ".") cat("\nRecurrent event process:")
    if (x$typeRec == "gsc") {
        p <- length(x$par1)
        mat <- rbind(c("Shape", rep("", p - 1), "Size", rep("", p - 1)),
                     rep(x$varNames, 2), format(c(x$par1, x$par2), digits = 5))
        prmatrix(mat, rowlab = rep("", nrow(mat)), collab = rep("", 1 + ncol(mat)), quote = FALSE)
    } else {
        mat <- rbind(x$varNames, format(x$par1, digits = 5))
        prmatrix(mat, rowlab = rep("", nrow(mat)), collab = rep("", ncol(mat)), quote = FALSE)
    }    
    if (length(x$par3) > 0) {
        cat("\nTerminal event:")
        if (x$typeRec == "gsc") {
            p <- length(x$par3)
            mat <- rbind(c("Shape", rep("", p - 1), "Size", rep("", p - 1)),
                         rep(x$varNames, 2), format(c(x$par3, x$par4), digits = 5))
            ## mat <- cbind(mat[,1:p], "    ", mat[,1:p])
            prmatrix(mat, rowlab = rep("", nrow(mat)), collab = rep("", 1 + ncol(mat)),
                     quote = FALSE)
        } else {
            mat <- rbind(x$varNames, format(x$par3, digits = 5))
            prmatrix(mat, rowlab = rep("", nrow(mat)), collab = rep("", ncol(mat)), quote = FALSE)
        }
        ## print.default(format(x$alpha, digits = digits), print.gap = 2L, quote = FALSE)
    } else {
        n <- length(unique(x$DF$id))
        nevent <- sum(x$DF$event > 0)
        avg.event <- nevent / n
        cat("\nNumber of subjects:", n)
        cat("\nNumber of recurrent events:", nevent)
        cat("\nAverage recurrent events per subject:", avg.event)
        cat("\n")
    }
}

pvalTab <- function(pe, se, names = NULL) {
    if (is.null(se)) se <- NA
    tab <- cbind(Estimate = round(pe, 5), StdErr = round(se, 5),
                 z.value = round(pe / se, 5), p.value = round(2 * pnorm(-abs(pe / se)), 5))
    if (!is.null(names)) rownames(tab) <- names
    return(tab)
}
    
#' @exportS3Method summary reReg
summary.reReg <- function(object, test = FALSE, ...) {
    if (!is.reReg(object)) stop("Must be a reReg object")
    if (object$typeRec == "nonparametric") {
        t0 <- sort(unique(c(object$DF$time1, object$DF$time2)))
        t0 <- t0[t0 > 0]
        out <- list(call = object$call, typeRec = object$typeRec, 
                    tabA = data.frame(time = t0, rate = object$Lam0(t0), hazard = object$Haz0(t0)))
        out
    }
    if (object$typeRec != "nonparametric") {
        tabA <- pvalTab(object$par1, object$par1.se, object$varNames)
        if (object$typeRec == "gsc") 
            tabA <- list(tabA1 = tabA,
                         tabA2 = pvalTab(object$par2, object$par2.se, object$varNames))
        out <- list(call = object$call, typeRec = object$typeRec, tabA = tabA)
        if (!is.null(object$par3))
            out$tabB <- pvalTab(object$par3, object$par3.se, object$varNames)
        if (object$typeTem == "gsc")
            out$tabB <- list(tabB1 = out$tabB,
                             tabB2 = pvalTab(object$par4, object$par4.se, object$varNames))
        if (object$typeRec == "gsc" & !is.null(object$par1.vcov) & !is.null(object$par2.vcov)) {
            p <- length(object$par1)
            out$HA.chi <- object$par1 %*% solve(object$par1.vcov) %*% object$par1
            out$HB.chi <- object$par2 %*%
                solve(object$par1.vcov + object$par2.vcov[-1, -1] +
                      2 * object$vcovRec[1:p, (p+2):(2*p+1), drop = FALSE]) %*%
                object$par2
            g <- object$par2 - object$par1
            out$HG.chi <- g %*% solve(object$par2.vcov[-1,-1]) %*% g
            out$HA.pval <- 1 - pchisq(out$HA.chi, p)
            out$HB.pval <- 1 - pchisq(out$HB.chi, p)
            out$HG.pval <- 1 - pchisq(out$HG.chi, p)
        }
        out$typeRec <- object$typeRec
        out$typeTem <- object$typeTem
        out$test <- test
    }
    class(out) <- "summary.reReg"
    return(out)
}

printCoefmat2 <- function(tab) 
    printCoefmat(as.data.frame(tab), P.values = TRUE, has.Pvalue = TRUE, signif.legend = FALSE)

#' @exportS3Method print summary.reReg
print.summary.reReg <- function(x, ...) {
    cat("Call: \n")
    dput(x$call)
    if (x$typeRec != "nonparametric" & !is.na(x$tabA)[1]) {
        if(x$typeRec == "cox.LWYY")
            cat("\nFitted with the Cox model of Lin et al. (2000):")
        if(x$typeRec == "cox.GL")
            cat("\nFitted with the Cox model of Ghosh and Lin (2002):")
        if(x$typeRec == "am.GL")
            cat("\nFitted with the accelerated mean model of Ghosh and Lin (2003):")
        if (x$typeRec == "gsc") {
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
        if (x$typeTem != ".") {
            if (x$typeTem == "gsc") {
                p <- nrow(x$tabB$tabB1)
                cat("\n\nTerminal event (shape):\n")
                printCoefmat2(x$tabB[[1]])
                cat("\nTerminal event (size):\n")
                printCoefmat2(x$tabB[[2]])
            } else {
                  cat("\nTerminal event:\n")
                  printCoefmat(x$tabB)
            }
        }
    }    
    if (x$typeRec == "nonparametric") {
        cat("\n")
        print(round(unique(x$tabA), 4), row.names = FALSE)
    }
    cat("\n")
}

#' @exportS3Method coef reReg
coef.reReg <- function(object, ...) {
    as.numeric(c(object$par1, object$par2, object$par3, object$par4))
}

#' @exportS3Method vcov reReg
vcov.reReg <- function(object, ...) {
    list(vcovRec = object$vcovRec, vcovTem = object$vcovTem)
}
