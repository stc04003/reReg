globalVariables(c("m", "x2", "Z")) ## global variables for simDat

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

#' Function to generate simulated data
#'
#' The function \code{simDat} generates simulated recurrent event data from either
#' a Cox-type model, an accelerated mean model, or a scale-change model.
#' The censoring time could be either independent (given covariates) or informative.
#' The simulated data is used for illustration.
#'
#' The function \code{simDat} generates simulated recurrent event data under different
#' scenarios based on the following assumptions.
#' See \bold{Details} in \code{\link{reReg}} for a more complete model assumptions.
#' \describe{
#'   \item{\code{type = "cox"}}{generates recurrent event data from a Cox-type model with
#' \deqn{\lambda(t) = Z \lambda_0(t) e^{X^\top a}, h(t) = Zh_0(t)e^{X^\top b}.}}
#'   \item{\code{type = "am"}}{generates recurrent event data from an accelerated mean model with
#' \deqn{\lambda(t) = Z \lambda_0(te^{X^\top a}) e^{X^\top a}, h(t) = Zh_0(te^{X^\top b})e^{X^\top b}.}}
#'   \item{\code{type = "sc"}}{generates recurrent event data from a generalized scale-change model with
#' \deqn{\lambda(t) = Z \lambda_0(te^{X^\top a}) e^{X^\top b}, h(t) = Zh_0(te^{X^\top a})e^{X^\top b}.}}
#' }
#' Let \eqn{D} be the informative failure time with the above hazard function.
#' An non-informative failure time, \eqn{C}, is generated separately from an
#' exponential distribution with mean 80.
#' The observed follow-up time is then taken to be \eqn{min(D, C, \tau)}.
#' We further assume
#' \deqn{\lambda_0(t) = \frac{2}{1 + t}, h_0(t) = \frac{1}{8(1 + t)}.}
#' Two covariates are considered; \code{x1} follows a Bernoulli distribution
#' with probability 0.5 and \code{x2} follows a standard normal distribution.
#' 
#' @param n number of observation.
#' @param a a numeric vector of parameter of length 2.
#' @param b a numeric vector of parameter of length 2.
#' @param indCen a logical value indicating whether the censoring assumption is imposed.
#' When \code{indCen = TRUE}, we set \eqn{Z = 1}.
#' Otherwise, \eqn{Z} is generated from a gamma distribution with mean 1 and variance 0.25
#' (e.g., \code{rgamma(1, 4, 4)}). See \bold{Details}.
#' @param type a character string specifying the underlying model. See \bold{Details}
#' @param tau a numeric value specifying the maximum observation time.
#' @param summary a logical value indicating whether a brief data summary will be printed.
#'
#' @seealso \code{\link{reReg}}
#' @export
#'
#' @importFrom tibble as_tibble
#'
#' @examples
#' set.seed(123)
#' simDat(200, c(-1, 1), c(-1, 1), summary = TRUE)
simDat <- function(n, a, b, indCen = TRUE, type = c("cox", "am", "sc"), tau = 60, summary = FALSE) {
    type <- match.arg(type)
    if (length(a) != 2L) stop("Require length(a) = 2.")
    if (length(b) != 2L) stop("Require length(b) = 2.")
    simOne <- function(id) {
        z <- ifelse(indCen, 1, rgamma(1, 4, 4))
        x <- c(sample(0:1, 1), rnorm(1))
        exa <- c(exp(x %*% a))
        exb <- c(exp(x %*% b))
        ## cen <- rexp(1, 1 / 80)
        cen <- rexp(1, x[1] / 60 + (1 - x[1]) * z^2 / 30)
        if (type == "cox") D <- invHaz(rexp(1), z, 1, exb)
        if (type == "am") D <- invHaz(rexp(1), z, exb, exb)
        if (type == "sc") D <- invHaz(rexp(1), z, exa, exb)
        y <- min(cen, tau, D) 
        status <- 1 * (y == D)
        m <- -1
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
                              Z = z, m = m, x1 = x[1], x2 = x[2]))
        } else {
            return(data.frame(id = id, Time = y, event = 0, status = status,
                              Z = z, m = m, x1 = x[1], x2 = x[2]))
        }
    }
    dat <- as_tibble(do.call(rbind, lapply(1:n, simOne)))
    if (summary) {
        cat("\n")
        cat("Summary results for number of recurrent event per subject:\n")
        base <- dat %>% group_by(id) %>% summarize(m = unique(m), d = max(status),
                                                   x1 = unique(x1), x2 = unique(x2))
        d <- base$d
        x1 <- base$x1
        print(summary(base$m))
        cat("\n")
        cat(paste("Number of failures: ", sum(d), " (", 100 * sum(d) / length(d), "%); ",
                  "Number of censored events: ", sum(d < 1), " (", 100 * sum(d < 1) / length(d), "%)\n\n",
                  sep = ""))
        cat(paste("Number of x1 == 1: ", sum(x1), " (", 100 * sum(x1) / length(x1), "%); ",
                  "Number of x1 == 0: ", sum(x1 < 1), " (", 100 * sum(x1 < 1) / length(x1), "%)\n", sep = ""))
        cat("Summary results for x2:\n")
        print(summary(base$x2))
        cat("\n\n")
    }
    return(dat %>% select(-m, -Z))
}
