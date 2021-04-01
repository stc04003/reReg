###################################################################################
## Codes to generate simulated data
###################################################################################

#' Function to find inverse of a given Lam
#' @noRd
inv <- function (t, z, exa, exb, fn) {
    mapply(t, FUN = function(u) {
        uf <- function(x) u - fn(x, z, exa, exb) ## / Lam.f(10, r, b, model)
        exp(optim(par = 1, fn = function(y) uf(exp(y))^2, 
                  control = list(warn.1d.NelderMead = FALSE))$par)
    })
}

#' Function to generate simulated recurrent event data
#'
#' The function \code{simSC()} generates simulated recurrent event data from either
#' a Cox-type model, an accelerated mean model, an accelerated rate model, or a scale-change model.
#'
#' 
#' The function \code{simSC()} generates simulated recurrent event data over
#' the interval \eqn{(0, \tau)} based on the specification of the recurrent process and
#' the terminal events.
#' Specifically, the rate function, \eqn{\lambda(t)}, of the recurrent process
#' can be specified as one of the following model:
#' \deqn{\lambda(t) = Z \lambda_0(te^{X^\top\alpha}) e^{X^\top\beta}, h(t) = Z h_0(te^{X^\top\eta})e^{X^\top\theta}, }
#' where \eqn{\lambda_0(t)} is the baseline rate function, \eqn{h_0(t)} is the baseline hazard function,
#' \eqn{X} is a \eqn{n} by \eqn{p} covariate matrix and \eqn{\alpha},
#' \eqn{Z} is an unobserved shared frailty variable, and
#' \eqn{(\alpha, \eta)} and \eqn{(\beta, \theta)} correspond to the shape and size parameters of the
#' rate function and the hazard function, respectively.
#'
#' For all scenarios, two covariates are considered; \eqn{X = (X_1, X_2)}, where
#' \eqn{X_1} follows a Bernoulli distribution with probability 0.5 and
#' \eqn{X_2} follows a standard normal distribution.
#' The censoring time could be either independent (given covariates) or informative.
#' The simulated data is used for illustration.
#' An informative censoring time, \eqn{C}, is generated separately from an
#' exponential distribution with a rate parameter of 1 / 60 if \eqn{X_1} is 1,
#' or \eqn{Z^2 / 30} if \eqn{X_1} is 0.
#' The observed recurrent events is then observed up to the minimum of \eqn{C},
#' terminal event, and \eqn{\tau}.
#' Lastly, when \code{lam0} is not specified, we assume the baseline rate function
#' to be \deqn{\lambda_0(t) = \frac{2}{1 + t}}.
#' On the other hand, when \code{haz0} is not specified, we assume
#' the baseline hazard function to be \deqn{h_0(t) = \frac{1}{8(1 + t)}}.
#' 
#' @param n number of observation.
#' @param alpha,beta,eta,theta are numerical vectors correspond to the \eqn{\alpha}, \eqn{\beta}, \eqn{\eta}, and \eqn{\theta} in the joint model, respectively. See \bold{Details}
#' @param censoring a numeric variable specifying the censoring times for each of the \eqn{n} observation.
#' @param xmat an optional matrix specifying the design matrix.
#' @param frailty a numeric variable specifying the frailty variable.
#' @param tau a numeric value specifying the maximum observation time.
#' @param origin a numeric value specifying the time origin.
#' @param Lam0 is an optional function that specifies the baseline cumulative rate function.
#' When left-unspecified, the recurrent events are generated using the
#' baseline rate function of \deqn{\lambda_0(t) = \frac{2}{1 + t},} or equivalently,
#' the cumulative rate function of \deqn{\Lambda_0(t) = 2\log(1 + t).}
#' @param Haz0 is an optional function that specifies the baseline hazard function.
#' When left-unspecified, the recurrent events are generated using the baseline hazard function
#' \deqn{h_0(t) = \frac{1}{8(1 + t)},} or equivalently,
#' the cumulative hazard function of \deqn{H_0(t) = \log(1 + t) / 8.}
#' @param summary a logical value indicating whether a brief data summary will be printed.
#'
#' @seealso \code{\link{reReg}}
#' @export
#'
#' @example inst/examples/ex_simu.R
simSC <- function(n, summary = FALSE,
                  alpha, beta, eta, theta,
                  xmat, censoring, frailty, tau, origin,
                  Lam0, Haz0) {
    call <- match.call()
    if (missing(tau)) tau <- 60
    if (missing(origin)) origin <- 0
    if (missing(frailty)) Z <- rgamma(n, 4, 4)
    else Z <- frailty
    if (missing(xmat)) {
        X <- cbind(sample(0:1, n, TRUE), rnorm(n, sd = .5))
        Cen <- runif(n, 0, X[,1] * tau * 2 + (1 - X[,1]) * 2 * Z^2 * tau)
    } else X <- xmat
    if (!missing(censoring)) Cen <- censoring
    p <- ncol(X)
    if (missing(alpha)) alpha <- rep(0, p)
    if (missing(eta)) eta <- rep(0, p)    
    if (missing(beta)) beta <- rep(-1, p)
    if (missing(theta)) theta <- rep(1, p)
    msg.mismatch <- function(x)
        paste("Parameter", substitute(x), "does not match with the number of covariates.")
    if (length(alpha) != p) stop(msg.mismatch(alpha))
    if (length(eta) != p) stop(msg.mismatch(eta))
    if (length(beta) != p) stop(msg.mismatch(beta))
    if (length(theta) != p) stop(msg.mismatch(theta))
    ## lapply(list(alpha, eta, beta, theta), mismatch, p = p)
    if (missing(Lam0)) {
        Lam <- function(t, z, exa, exb) z * exb * log(1 + t * exa) / exa / .5
        invLam <- function(t, z, exa, exb) (exp(.5 * t * exa / exb / z) - 1) / exa
    } else {
        Lam <- function(t, z, exa, exb) z * Lam0(t * exa) * exb / exa
        invLam <- function(t, z, exa, exb) inv(t, z, exa, exb, Lam)
    }
    if (missing(Haz0)) {
        invHaz <- function(t, z, exa, exb) (exp(5 * t * exa / exb / z) - 1) / exa
    } else {
        Haz <- function(t, z, exa, exb) z * Haz0(t * exa) * exb / exa
        invHaz <- function(t, z, exa, exb) inv(t, z, exa, exb, Haz)
    }
    if (n != length(origin) & length(origin) > 1)
        stop("Invalid length for 'origin'. See '?simSC' for details.")
    simOne <- function(id, z, x, cen) {
        D <- invHaz(rexp(1), z, c(exp(x %*% eta)), c(exp(x %*% theta)))
        y <- min(cen, tau, D)
        status <- 1 * (y == D)
        m <- -1
        tij <- NULL
        up <- Lam(y, z, c(exp(x %*% alpha)), c(exp(x %*% beta)))
        while(sum(tij) < up) {
            tij <- c(tij, rexp(1))
            m <- m + 1
        }
        if (m > 0) {
            tij <- invLam(cumsum(tij[1:m]), z, c(exp(x %*% alpha)), c(exp(x %*% beta)))
            return(data.frame(id = id, Time = c(sort(tij), y),
                              event = c(rep(1, m), 0), status = c(rep(0, m), status),
                              Z = z, m = m, x = t(x)))
        } else {
            return(data.frame(id = id, Time = y, event = 0, status = status,
                              Z = z, m = m, x = t(x)))
        }
    }
    dat <- data.frame(do.call(rbind, lapply(1:n, function(i) simOne(i, Z[i], X[i,], Cen[i]))))
    names(dat)[grep("x.", names(dat))] <- paste0("x", 1:p)
    if (length(origin) > 1) origin <- rep(origin, unlist(lapply(split(dat$id, dat$id), length)))
    dat$t.start <- do.call(c, lapply(split(dat$Time, dat$id), function(x)
        c(0, x[-length(x)]))) + origin
    dat$t.stop <- dat$Time + origin
    dat$Time <- NULL
    if (summary) {
        dg <- min(3, getOption("digits"))
        cat("Call: \n")
        print(call)
        cat("\n")
        cat("Summary:\n")
        cat("Sample size:                                   ", n, "\n")
        cat("Number of recurrent event observed:            ", sum(dat$event), "\n")
        cat("Average number of recurrent event per subject: ", round(sum(dat$event) / n, dg), "\n")
        cat("Proportion of subjects with a terminal event:  ", round(sum(dat$status) / n, dg), "\n")
        base <- dat[cumsum(table(dat$id)), ]
        y <- base$t.stop
        d <- base$status
        oy <- order(y)
        d <- d[oy]
        y <- y[oy]
        r <- n - rank(y, ties.method = "min") + 1
        is_first_y <- ! duplicated(y)
        s <- cumprod(1 - (d / r)[is_first_y])
        medTem <- ifelse(s[length(s)] > .5, NA, as.numeric(y[is_first_y][which.max(s - .5 < 0)]))
        if (!is.na(medTem))
            cat("Median time-to-terminal event:                 ",
                round(medTem, dg), "\n")            
        ## cat("Proportion of subjects with a x1 = 1:  ", round(mean(base$x1), dg), "\n")
        cat("\n\n")
    }
    dat$m <- dat$Z <- NULL
    ## dat <- dat[,c(1, 6:7, 2:5)]
    ## Reorder columns
    ord <- c("id", "t.start", "t.stop", "event", "status")
    dat <- dat[,c(ord, setdiff(names(dat), ord))]
    attr(dat, "Call") <- call
    return(dat)
}
               
