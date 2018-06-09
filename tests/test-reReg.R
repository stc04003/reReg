library(reReg)

## -----------------------------------------------------------------------------------
## With simulated data
## -----------------------------------------------------------------------------------

simDat <- function(n, a, b, latent = FALSE) {
    ## setting rate function
    Lam.f <- function(t, z, x, w) .5 * z * exp(-x + w) * log(1 + t * exp(x))
    Lam.f0 <- function(t) .5 * log(1 + t)
    invLam.f  <- function(t, z, x, w) (exp((2 / z) * exp(x - w) * t )- 1) / exp(x)
    ## setting hazard funciton
    ## Haz.f0 <- function(t) .5 * log(1 + t) # assume constant hazard for now
    dat <- NULL
    for (id in 1:n) {
        z <- ifelse(latent, rgamma(1, 4, 4), 1)
        x1 <- rnorm(1)
        x2 <- rnorm(1)
        x <- c(x1, x2)
        cen <- rexp(1, z * exp(x %*% b) / 60) ## this gives constant hazard of 1/60
        y <- min(cen, 60)
        D <- 1 * (cen == y)
        tmpt <- NULL
        while(sum(tmpt) < Lam.f(y, z, c(x %*% a), c(x %*% b))) {
            tmpt <- c(tmpt, rexp(1))
        }
        m <- length(tmpt) - 1
        if (m > 0) {
            tt <- invLam.f(cumsum(tmpt[1:m]), z, c(x %*% a), c(x %*% b))
            dat <- rbind(dat, cbind(ID = id, Time = c(tt[order(tt)], y),
                                    event = c(rep(1, m), 0), status = c(rep(0, m), D),
                                    Z = z, M = m, X1 = x1, X2 = x2))
        } else {
            dat <- rbind(dat, cbind(ID = id, Time = y, event = 0, status = D,
                                    Z = z, M = m, X1 = x1, X2 = x2))
        }
    }
    return(data.frame(dat))
}

dat <- simDat(200, 


## -----------------------------------------------------------------------------------
## With readmission data
## -----------------------------------------------------------------------------------
data(readmission, package = "frailtypack")
fm <- reSurv(t.stop, id, event, death) ~ sex + chemo
