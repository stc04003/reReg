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
                    coefficients.rec = data.frame(time = t0, rate = object$Lam0(t0)))
        if (!is.null(object$Haz0))
            out$coefficients.rec$hazard <- object$Haz0(t0)
        out$typeRec <- object$typeRec
        out
    }
    if (object$typeRec != "nonparametric") {
        coefficients.rec <- pvalTab(object$par1, object$par1.se, object$varNames)
        if (object$typeRec == "gsc") 
            coefficients.rec <-
                list(coefficients.shape = coefficients.rec,
                     coefficients.size = pvalTab(object$par2, object$par2.se, object$varNames))
        out <- list(call = object$call, typeRec = object$typeRec, coefficients.rec = coefficients.rec)
        if (!is.null(object$par3))
            out$coefficients.haz <- pvalTab(object$par3, object$par3.se, object$varNames)
        if (object$typeTem == "gsc")
            out$coefficients.haz <-
                list(coefficients.shape = out$coefficients.haz,
                     coefficients.size = pvalTab(object$par4, object$par4.se, object$varNames))
        if (object$typeRec == "gsc" & !is.null(object$par1.vcov) & !is.null(object$par2.vcov)) {
            p <- length(object$par1)
            out$HA.chi <- bAib(object$par1.vcov, object$par1)
            out$HB.chi <- bAib(object$par2.vcov, object$par2)
            g <- object$par2 - object$par1
            out$HG.chi <- bAib(object$par2.vcov - object$par1.vcov - 2 * object$vcovRec12, g)
            out$HA.pval <- 1 - pchisq(out$HA.chi, p)
            out$HB.pval <- 1 - pchisq(out$HB.chi, p)
            out$HG.pval <- 1 - pchisq(out$HG.chi, p)
        }
        if (object$typeTem == "gsc" & !is.null(object$par3.vcov) & !is.null(object$par4.vcov)) {
            p <- length(object$par3)
            out$tem.HA.chi <- bAib(object$par3.vcov, object$par3)
            out$tem.HB.chi <- bAib(object$par4.vcov, object$par4)
            g <- object$par4 - object$par3
            out$tem.HG.chi <- bAib(object$par4.vcov - object$par3.vcov - 2 * object$vcovTem12, g)
            out$tem.HA.pval <- 1 - pchisq(out$tem.HA.chi, p)
            out$tem.HB.pval <- 1 - pchisq(out$tem.HB.chi, p)
            out$tem.HG.pval <- 1 - pchisq(out$tem.HG.chi, p)
        }        
        out$typeRec <- object$typeRec
        out$typeTem <- object$typeTem
        out$test <- test
    }
    out$vcov <- vcov(object)
    class(out) <- "summary.reReg"
    return(out)
}

#' @exportS3Method print summary.reReg
print.summary.reReg <- function(x, ...) {
    if (x$typeRec != "nonparametric" & !is.na(x$coefficients.rec)[1]) {
        cat("Call: \n")
        dput(x$call)
        if(x$typeRec == "cox.LWYY")
            cat("\nFitted with the Cox model of Lin et al. (2000):")
        if(x$typeRec == "cox.GL")
            cat("\nFitted with the Cox model of Ghosh and Lin (2002):")
        if(x$typeRec == "am.GL")
            cat("\nFitted with the accelerated mean model of Ghosh and Lin (2003):")
        if (x$typeRec == "gsc") {
            p <- nrow(x$coefficients.rec$coefficients.shape)
            cat("\nRecurrent event process (shape):\n")
            printCoefmat2(x$coefficients.rec$coefficients.shape)
            cat("\nRecurrent event process (size):\n")
            printCoefmat2(x$coefficients.rec$coefficients.size)
            if (x$test) {
                cat("\nHypothesis tests:")
                cat("\nHo: shape = 0 (Cox-type model):")
                cat(paste("\n     X-squared = ", round(x$HA.chi, 4), ", df = ", p,
                          ", p-value = ", round(x$HA.pval, 4), sep = ""))
                cat("\nHo: size = 0 (Accelerated rate model):")
                cat(paste("\n     X-squared = ", round(x$HB.chi, 4), ", df = ", p,
                          ", p-value = ", round(x$HB.pval, 4), sep = ""))
                cat("\nHo: shape = size (Accelerated mean model):")
                cat(paste("\n     X-squared = ", round(x$HG.chi, 4), ", df = ", p,
                          ", p-value = ", round(x$HG.pval, 4), sep = ""))
            }
        } else {
            cat("\nRecurrent event process:\n")
            printCoefmat2(x$coefficients.rec)
        }
        ## Lin-Wei-Yang-Ying method (fitted with coxph with robust variance)
        if (x$typeTem != ".") {
            if (x$typeTem == "gsc") {
                p <- nrow(x$coefficients.haz$coefficients.shape)
                cat("\n\nTerminal event (shape):\n")
                printCoefmat2(x$coefficients.haz$coefficients.shape)
                cat("\nTerminal event (size):\n")
                printCoefmat2(x$coefficients.haz$coefficients.size)
                if (x$test) {
                    cat("\nHypothesis tests:")
                    cat("\nHo: shape = 0 (Cox-type model):")
                    cat(paste("\n     X-squared = ", round(x$tem.HA.chi, 4), ", df = ", p,
                              ", p-value = ", round(x$tem.HA.pval, 4), sep = ""))
                    cat("\nHo: size = 0 (Accelerated rate model):")
                    cat(paste("\n     X-squared = ", round(x$tem.HB.chi, 4), ", df = ", p,
                              ", p-value = ", round(x$tem.HB.pval, 4), sep = ""))
                    cat("\nHo: shape = size (Accelerated mean model):")
                    cat(paste("\n     X-squared = ", round(x$tem.HG.chi, 4), ", df = ", p,
                              ", p-value = ", round(x$tem.HG.pval, 4), sep = ""))
                }
            } else {
                cat("\nTerminal event:\n")
                printCoefmat2(x$coefficients.haz)
            }
        }
    }    
    if (x$typeRec == "nonparametric") {
        cat("\nNonparametric estimation:\n")
        print(head(x$coefficients.rec))
    }
    cat("\n")
}

#' @exportS3Method coef reReg
coef.reReg <- function(object, ...) {
    as.numeric(c(object$par1, object$par2, object$par3, object$par4))
}

#' @exportS3Method vcov reReg
vcov.reReg <- function(object, ...) {
    vcovRec <- vcovTem <- NULL
    if (is.null(object$par1.vcov))
        return(list(vcovRec = vcovRec, vcovTem = vcovTem))
    if (object$typeRec == "cox") {
        vcovRec <- object$par1.vcov[-1, -1, drop = FALSE]
        attr(vcovRec, "dimnames") <- list(object$varNames, object$varNames)
    }        
    if (object$typeRec %in% c("am", "ar")) {
        vcovRec <- object$par1.vcov
        attr(vcovRec, "dimnames") <- list(object$varNames, object$varNames)
    }
    if (object$typeRec == "gsc") {
        vcovRec.shape <- object$par1.vcov
        vcovRec.size <- object$par2.vcov
        attr(vcovRec.shape, "dimnames") <- list(object$varNames, object$varNames)
        attr(vcovRec.size, "dimnames") <- list(object$varNames, object$varNames)
        vcovRec <- list(vcovRec.shape = vcovRec.shape, vcovRec.size = vcovRec.size)
    }
    if (object$typeTem == ".") return(vcovRec)    
    if (object$typeTem != "gsc") {
        vcovTem <- object$par3.vcov
        attr(vcovTem, "dimnames") <- list(object$varNames, object$varNames)
    }
    if (object$typeTem == "gsc") {
        vcovTem.shape <- object$par3.vcov
        vcovTem.size <- object$par4.vcov
        attr(vcovTem.shape, "dimnames") <- list(object$varNames, object$varNames)
        attr(vcovTem.size, "dimnames") <- list(object$varNames, object$varNames)
        vcovTem <- list(vcovTem.shape = vcovTem.shape, vcovTem.size = vcovTem.size)        
    }
    list(vcovRec = vcovRec, vcovTem = vcovTem)
}

#' @exportS3Method coef summary.reReg
coef.summary.reReg <- function(object, ...) {
    if (is.null(object$coefficients.haz))
        return(object$coefficients.rec)
    else
        return(list(coefficients.rec = object$coefficients.rec,
                    coefficients.haz = object$coefficients.haz))
}


#' @exportS3Method vcov summary.reReg
vcov.summary.reReg <- function(object, ...) {
    object$vcov
}
