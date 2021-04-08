#' Calculate Residuals for a `reReg' Fit
#'
#' Calculates residuals for a joint frailty scale-change model fitted by `reReg`.
#' Under the recurrent event model, at each observation time, \eqn{t},
#' the residual is calculated as
#' \deqn{\mbox{observed number of recurrent events at } t- 
#' \mbox{expected number of recurrent events at} t.}
#' The expected number of recurrent events at \eqn{t} is calculated by the
#' cumulative rate function at \eqn{t}.
#' Under the failure time model, the residual is calculated as
#' \deqn{\Delta - H(t),}
#' where \eqn{\Delta} is the terminal event indicator and
#' \eqn{H(t)} is the cumulative hazard function at \eqn{t}.
#'
#' @param object an object of class \code{reReg} returned by the \code{reReg()} function.
#' @param model a character string specifying whether the residuals will be calculated
#' under the recurrent event model or the failure time model.
#' @param ... additional parameters for future development.
#' 
#' @exportS3Method residuals reReg
residuals.reReg <- function(object, model = c("recurrent", "failure"), ...) {
    if (!is.reReg(object)) stop("Must be a reReg object")
    model <- match.arg(model)
    if (object$typeRec %in% c("cox.LWYY", "cox.GL", "am.GL"))
        stop("Residuals calculation not available.")
    if (object$typeTem == "." & model == "failure")
        stop("Residuals on failure model not available. ")
    t0 <- object$DF$time2
    X <- as.matrix(subset(object$DF, select = object$varNames))
    mi <- unlist(lapply(split(object$DF$id, object$DF$id), length))
    O <-  unlist(sapply(mi, function(x) 1:x)) - 1
    names(O) <- NULL
    p <- length(object$varNames)
    if (model == "recurrent") {
        exa1 <- exa2 <- 1
        if (object$typeRec == "cox") exa2 <- exp(X %*% object$par1)
        if (object$typeRec == "ar") exa1 <- exp(X %*% object$par1)
        if (object$typeRec == "am") exa1 <- exa2 <- exp(X %*% object$par1)
        if (object$typeRec == "gsc") {
            exa1 <- exp(X %*% object$par1)
            exa2 <- exp(X %*% object$par2)
        }
        E <- drop(rep(object$zi, mi) * object$Lam0(t0 * exa1) * exa2 / exa1)
        names(E) <- NULL
        return(O - E)
    }
    if (model == "failure") {
        exb1 <- exb2 <- 1
        if (object$typeTem == "cox") exb2 <- exp(X %*% object$par3)
        if (object$typeTem == "ar") exb1 <- exp(X %*% object$par3)
        if (object$typeTem == "am") exb1 <- exb2 <- exp(X %*% object$par3)
        if (object$typeTem == "gsc") {
            exb1 <- exp(X %*% object$par3)
            exb2 <- exp(X %*% object$par4)
        }
        E <- drop(rep(object$zi, mi) * object$Haz0(t0 * exb1) * exb2 / exb1)
        names(E) <- NULL
        return((O - E) * (1 - object$DF$event))
    }
}
