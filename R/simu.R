###################################################################################
## Codes to generate simulated data
###################################################################################

Lam <- function(t, z, exa, exb) z * exb * log(1 + t * exa) / exa / 2
Lam0 <- function(t) log(1 + t) / 2
invLam <- function(t, z, exa, exb) exp(2 * t * exb / exa / z) / exa
invHaz <- function(t, z, exa, exb) exp(4 * t * exb / exa / z) / exa

#' Function to generate simulated data
#'
#' The function \code{simDat} generates simulated recurrent event data from either
#' a Cox-type model, an accelerated mean model, or a scale-change model.
#' The censoring time could be either independent or informative.
#' The simulated data is used for illustraction.
#'
#' The function \code{simDat} generates simulated recurrent event data under different
#' scenarios based on the following assumptions.
#' Suppose the recurrent event process and the failure events are observed in the time
#' interval \eqn{t\in[0,\tau]}, for some constant \eqn{\tau}.
#' Let ...
#' 
#' When \code{type == "cox"}, we assume the rate function of the counting process takes the form
#' \deqn{\lambda(t) = Z \lambda_0(t) X^\top\alpha},
#' and the hazard function of the failure time satisfies
#' \deqn{h(t) = Z h_0(t) X^\top\alpha}.
#'
#' When \code{type == "am"}, we assume the rate function and the hazard function are
#' \deqn{\lambda(t) = Z \lambda_0(te^{X^\top\alpha}) X^\top\beta},
#' and
#' \deqn{h(t) = Z h_0(te^{X^\top\alpha}) X^\top\beta},
#'
#' When \code{type = "sc"}, we assume the rate function has the form
#' \deqn{\lambda(t) = Z \lambda_0(tX^\top\alpha) X^\top\beta},
#' and ...
#'
#' For all three types, we assume
#' \deqn{\lambda_0(t) = \frac{1}{2 (1 + t)}}.
#' For illustration, we generate \eqn{X} from a uncorrelated bivariate standard normal distribution.
#' 
#' @param n number of observation.
#' @param a a numeric vector of parameter of length 2.
#' @param b a numeric vector of parameter of length 2.
#' @param indCen a logical value indicating whether independent censoring assumption is imposed.
#' When \code{indCen = TRUE}, we set \eqn{Z = 1}.
#' Otherwise, \eqn{Z} is generated from a gamma distribution with mean 1 and variance 0.25
#' (e.g., \code{rgamma(1, 4, 4)}).
#' @param type a character string specifying the underlying model. See \code{Details}.
#' @param tau a numeric value specifying the maximum observation time.
#'
#' @export
#' 
simDat <- function(n, a, b, indCen = TRUE, type = c("cox", "am", "sc"), tau = 60) {
    type <- match.arg(type)
    if (length(a) != 2L) stop("Require length(a) = 2.")
    if (length(b) != 2L) stop("Require length(b) = 2.")
    simOne <- function(id) {
        z <- ifelse(indCen, 1, rgamma(1, 4, 4))
        x <- rnorm(2)
        exa <- exp(x %*% a)
        exb <- exp(x %*% b)
        cen <- rexp(1, z * exa / 60)
        if (type == "cox") D <- invHaz(runif(1), z, 1, exb)
        if (type == "am") D <- invHaz(runif(1), z, exb, exb)
        if (type == "sc") D <- tau
        y <- min(cen, tau, D) 
        status <- 1 * (y == D)
        m <- 0
        tij <- NULL
        if (type == "cox") up <- Lam(y, z, 1, exa)
        if (type == "am") up <- Lam(y, z, exa, exa)
        if (type == "sc") up <- Lam(y, z, exa, exb)
        while(sum(tij) < up) {
            tij <- c(tij, rexp(1))
            m <- m + 1
        }
        if (m > 0) {
            if (type == "cox") tij <- invLam(cumsum(tij[1:m]), z, 1, exa)
            if (type == "am") tij <- invLam(cumsum(tij[1:m]), z, exa, exa)
            if (type == "sc") tij <- invLam(cumsum(tij[1:m]), z, exa, exb)
            return(data.frame(id = id, Time = c(sort(tij), y),
                              event = c(rep(1, m), 0), status = c(rep(0, m), status),
                              Z = z, m = m, X1 = x[1], X2 = x[2]))
        } else {
            return(data.frame(id = id, Time = y, event = 0, status = status,
                              Z = z, m = m, X1 = x[1], X2 = x[2]))
        }
    }
    tibble(do.call(rbind, lapply(1:n, simOne)))
}
