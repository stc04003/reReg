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
#' The function \code{simGSC()} generates simulated recurrent event data from either
#' a Cox-type model, an accelerated mean model, an accelerated rate model, or a generalized scale-change model.
#'
#' 
#' The function \code{simGSC()} generates simulated recurrent event data over
#' the interval \eqn{(0, \tau)} based on the specification of the recurrent process and
#' the terminal events.
#' Specifically, the rate function, \eqn{\lambda(t)}, of the recurrent process
#' can be specified as one of the following model:
#' \deqn{\lambda(t) = Z \lambda_0(te^{X^\top\alpha}) e^{X^\top\beta}, h(t) = Z h_0(te^{X^\top\eta})e^{X^\top\theta}, }
#' where \eqn{\lambda_0(t)} is the baseline rate function,
#' \eqn{h_0(t)} is the baseline hazard function,
#' \eqn{X} is a \eqn{n} by \eqn{p} covariate matrix and \eqn{\alpha},
#' \eqn{Z} is an unobserved shared frailty variable, and
#' \eqn{(\alpha, \eta)} and \eqn{(\beta, \theta)} correspond to the shape and size parameters of the
#' rate function and the hazard function, respectively.
#'
#' Under the default settings, the \code{simGSC()} function assumes \eqn{p = 2}
#' and the regression parameters to be \eqn{\alpha = \eta = (0, 0)^\top},
#' and \eqn{\beta = \theta = (1, 1)^\top}.
#' When the \code{xmat} argument is not specified, the \code{simGSC()} function
#' assumes \eqn{X_i} is a two-dimensional vector \eqn{X_i = (X_{i1}, X_{i2}), i = 1, \ldots, n},
#' where \eqn{X_{i1}} is a Bernoulli variable with rate 0.5 and
#' \eqn{X_{i2}} is a standard normal variable.
#' With the default \code{xmat}, the censoring time $C$ is generated from
#' an independent uniform distribution in \eqn{[0, 2\tau X_{i1} + 2Z^2\tau(1 - X_{i1})]}.
#' Thus, the censoring distribution is covariate dependent and
#' is informative when \eqn{Z} is not a constant.
#' When the \code{frailty} argument is not specified, the frailty variable \eqn{Z} is generated
#' from a gamma distribution with a unit mean and a variance of 0.25.
#' The default values for \code{tau} and \code{origin} are 60 and 0, respectively.
#' When arguments \code{Lam0} and \code{Haz0} are left unspecified,
#' the \code{simGSC()} function uses \eqn{\Lambda_0(t) = 2\log(1 + t)}
#' and \eqn{H_0(t) = \log(1 + t) / 5}, respectively.
#' This is equivalent to setting
#' \code{Lam0 = function(x) 2 * log(1 + x)} and \code{Haz0 = function(x) log(1 + x) / 5}.
#' Overall, the default specifications generate the recurrent events and the terminal events
#' from the model:
#' \deqn{\lambda(t) = \displaystyle \frac{2Z}{1 + te^{-X_{i1} - X_{i2}}},
#' h(t) = \displaystyle \frac{Z}{5(1 + te^{X_{i1} + X_{i2}})},  t\in[0, 60].}
#'
#' 
#' @param n number of observation.
#' @param para a list of numerical vectors for the regression coefficients
#' in the joint scale-change model. 
#' The names of the list elements are \code{alpha}, \code{beta}, \code{eta}, and
#' \code{theta}, correspond to \eqn{\alpha}, \eqn{\beta}, \eqn{\eta}, and \eqn{\theta}
#' in the joint scale-change model, respectively.
#' See \bold{Details} for \code{\link{reReg}}.
#' @param censoring a numeric variable specifying the censoring times for each of the
#' \eqn{n} observation.
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
#' \deqn{h_0(t) = \frac{1}{5(1 + t)},} or equivalently,
#' the cumulative hazard function of \deqn{H_0(t) = \log(1 + t) / 5.}
#' @param summary a logical value indicating whether a brief data summary will be printed.
#'
#' @seealso \code{\link{reReg}}
#' @export
#'
#' @example inst/examples/ex_simu.R
simGSC <- function(n, summary = FALSE, para,
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
    } else {
        if (!missing(censoring)) Cen <- censoring
        if (missing(censoring)) Cen <- runif(n, 0, 2 * tau)
        X <- xmat
    }
    
    p <- ncol(X)
    para0 <- list(alpha = rep(0, p), beta = rep(-1, p), eta = rep(0, p), theta = rep(1, p))
    if (missing(para)) para <- para0
    namel <- names(para)
    para0[namel] <- para
    alpha <- para0$alpha
    beta <- para0$beta
    eta <- para0$eta
    theta <- para0$theta
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
        stop("Invalid length for 'origin'. See '?simGSC' for details.")
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
               
