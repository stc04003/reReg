globalVariables(c("m", "x2", "Z")) ## global variables for simSC

###################################################################################
## Codes to generate simulated data
###################################################################################

Lam <- function(t, z, exa, exb) z * exb * log(1 + t * exa) / exa / .5
Lam0 <- function(t) log(1 + t) / .5
invLam <- function(t, z, exa, exb) (exp(.5 * t * exa / exb / z) - 1) / exa
invHaz <- function(t, z, exa, exb) (exp(8 * t * exa / exb / z) - 1) / exa

## k is the scaling parameter
## Lam <- function(t, z, exa, exb) z * exb * log(1 + t * exa) / exa / k
## Lam0 <- function(t) log(1 + t) / k
## invLam <- function(t, z, exa, exb) (exp(k * t * exa / exb / z) - 1) / exa

#' Function to generate simulated recurrent event data
#'
#' The function \code{simSC} generates simulated recurrent event data from either
#' a Cox-type model, an accelerated mean model, or a scale-change model.
#' For all scenarios, two covariates are considered; \eqn{X = (X_1, X_2)}, where
#' \eqn{X_1} follows a Bernoulli distribution with probability 0.5 and
#' \eqn{X_2} follows a standard normal distribution.
#' The censoring time could be either independent (given covariates) or informative.
#' The simulated data is used for illustration.
#'
#' The function \code{simSC()} generates simulated recurrent event data over
#' the interval \eqn{(0, \tau)} based on the specification of the recurrent process and
#' the terminal events.
#' Specifically, the rate function, \eqn{\lambda(t)}, of the recurrent process
#' can be specified as one of the following model:
#' \describe{
#'   \item{Cox-type model}{
#' \deqn{\lambda(t) = Z \lambda_0(t) e^{X^\top a}.}}
#'   \item{Accelerated mean model}{
#' \deqn{\lambda(t) = Z \lambda_0(te^{X^\top a}) e^{X^\top a}.}}
#'   \item{Scale-change model}{
#' \deqn{\lambda(t) = Z \lambda_0(te^{X^\top a}) e^{X^\top b}.}}
#' }
#' Similarly, the hazard function, \eqn{h(t)}, of the terminal event 
#' can be specified as one of the following model:
#' \describe{
#'   \item{Cox-type model}{
#' \deqn{h(t) = Zh_0(t)e^{X^\top b}.}}
#'   \item{Accelerated mean model}{
#' \deqn{h(t) = Zh_0(te^{X^\top b})e^{X^\top b}.}}
#'   \item{Scale-change model}{
#' \deqn{h(t) = Zh_0(te^{X^\top a})e^{X^\top b}.}}
#' }
#' An informative censoring time, \eqn{C}, is generated separately from an
#' exponential distribution with a rate parameter of 1 / 60 if \eqn{X_1} is 1,
#' or \eqn{Z^2 / 30} if \eqn{X_1} is 0.
#' The observed recurrent events is then observed up to the minimum of \eqn{C},
#' terminal event, and \eqn{\tau}.
#' Lastly, we assume the baseline functions
#' \deqn{\lambda_0(t) = \frac{2}{1 + t}, h_0(t) = \frac{1}{8(1 + t)}.}
#' 
#' @param n number of observation.
#' @param a a numeric vector of parameter of length 2.
#' @param b a numeric vector of parameter of length 2.
#' @param type a character string specifying the underlying model.
#' See \bold{Details}
#' @param zVar a numeric variable specifying the variance of the fraility variable,\eqn{Z}, when \code{zVar} > 0.
#' When \code{zVar} = 0, \eqn{Z} is set to a fixed constant 1. The default value is 0.25.
#' @param tau a numeric value specifying the maximum observation time.
#' @param summary a logical value indicating whether a brief data summary will be printed.
#'
#' @seealso \code{\link{reReg}}
#' @export
#'
#'
#' @example inst/examples/ex_simu.R
simSC <- function(n, a, b, type = "cox", zVar = .25, tau = 60, summary = FALSE) {
    if (length(a) != 2L) stop("Require length(a) = 2.")
    if (length(b) != 2L) stop("Require length(b) = 2.")
    allcomb <- apply(expand.grid(c("cox", "am", "sc"), c("cox", "am", "sc")), 1, paste, collapse = "|")
    type <- match.arg(type, c("cox", "am", "sc", allcomb))
    if (grepl("|", type, fixed = TRUE)) {
        recType <- substring(type, 1, regexpr("[|]", type) - 1)
        temType <- substring(type, regexpr("[|]", type) + 1)
    } else {
        recType <- temType <- type
    }
    if (zVar <= 0) Z <- rep(1, n)
    else Z <- rgamma(n, 1/zVar, 1/zVar)
    X <- cbind(sample(0:1, n, TRUE), rnorm(n))
    Cen <- sapply(1:n, function(x) rexp(1, X[x, 1] / 60 + (1 - X[x, 1]) * Z[x]^2 / 30))
    rr <- rexp(n)
    simOne <- function(id, z, x, cen, rr) {
        exa <- c(exp(x %*% a))
        exb <- c(exp(x %*% b))
        if (temType == "cox") D <- invHaz(rr, z, 1, exb)
        if (temType == "am") D <- invHaz(rr, z, exb, exb)
        if (temType == "sc") D <- invHaz(rr, z, exa, exb)
        y <- min(cen, tau, D) 
        status <- 1 * (y == D)
        m <- -1
        tij <- NULL
        if (recType == "cox") up <- Lam(y, z, 1, exa)
        if (recType == "am") up <- Lam(y, z, exa, exa)
        if (recType == "sc") up <- Lam(y, z, exa, exb)
        while(sum(tij) < up) {
            tij <- c(tij, rexp(1))
            m <- m + 1
        }
        if (m > 0) {
            if (recType == "cox") tij <- invLam(cumsum(tij[1:m]), z, 1, exa)
            if (recType == "am") tij <- invLam(cumsum(tij[1:m]), z, exa, exa)
            if (recType == "sc") tij <- invLam(cumsum(tij[1:m]), z, exa, exb)
            return(data.frame(id = id, Time = c(sort(tij), y),
                              event = c(rep(1, m), 0), status = c(rep(0, m), status),
                              Z = z, m = m, x1 = x[1], x2 = x[2]))
        } else {
            return(data.frame(id = id, Time = y, event = 0, status = status,
                              Z = z, m = m, x1 = x[1], x2 = x[2]))
        }
    }
    dat <- data.frame(do.call(rbind, lapply(1:n, function(y) simOne(y, Z[y], X[y,], Cen[y], rr[y]))))
    if (summary) {
        cat("\n")
        cat("Summary results for number of recurrent event per subject:\n")
        dat <- dat[order(dat$id),]
        base <- dat[cumsum(table(dat$id)),] 
        d <- base$status
        x1 <- base$x1
        print(summary(base$m))
        cat("\n")
        cat(paste("Number of failures: ", sum(d), " (", round(100 * sum(d) / length(d), 2), "%); ",
                  "Number of censored events: ", sum(d < 1), " (", round(100 * sum(d < 1) / length(d), 2), "%)\n\n",
                  sep = ""))
        cat(paste("Number of x1 == 1: ", sum(x1), " (", round(100 * sum(x1) / length(x1), 2), "%); ",
                  "Number of x1 == 0: ", sum(x1 < 1), " (", round(100 * sum(x1 < 1) / length(x1), 2), "%)\n", sep = ""))
        cat("Summary results for x2:\n")
        print(summary(base$x2))
        cat("\n\n")
    }
    return(subset(dat, select = c(-m, -Z)))
}
