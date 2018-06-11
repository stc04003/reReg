library(reReg)

## -----------------------------------------------------------------------------------
## With simulated data
## -----------------------------------------------------------------------------------

args(simDat)
fm <- reSurv(Time, id, event, status) ~ x1 + x2

simDat(200, c(1, 1), c(1, 1))
simDat(200, c(1, 1), c(1, 1), type = "am")
simDat(200, c(1, 1), c(1, 1), type = "sc")

simDat(200, c(1, 1), c(1, 1), indCen = FALSE)
simDat(200, c(1, 1), c(1, 1), type = "am", indCen = FALSE)
simDat(200, c(1, 1), c(1, 1), type = "sc", indCen = FALSE)


dat <- simDat(200, c(1, 1), c(1, 1))
dat
with(dat, reSurv(Time, id, event, status)$reDF)
with(dat, reSurv(Time, id, event, status)$reTb)

coef(reReg(fm, data = dat))


do <- function(n = 200, a = c(1, 1), b = c(1, 1), type = "cox", indCen = TRUE) {
    dat <- simDat(n, a, b, indCen, type)
    c1 <- coef(reReg(fm, data = dat))
    c2 <- coef(reReg(fm, data = dat, method = "cox.HW"))
    c3 <- coef(reReg(fm, data = dat, method = "am.GL"))
    c4 <- coef(reReg(fm, data = dat, method = "am.XCHWY"))
    c5 <- coef(reReg(fm, data = dat, method = "sc.XCYH"))
    c(c1, c2, c3, c4, c5)
}

set.seed(1)
do(a = c(-1, -1), b = c(1, 1), indCen = TRUE, type = "cox")

e



set.seed(1527) 
do(a = c(-1, -1), b = c(1, 1), indCen = TRUE, type = "sc")
set.seed(1526)  ## and 1527
dat <- simDat(200, a = c(-1, -1), b = c(1, 1), indCen = TRUE, type = "sc")
coef(reReg(fm, data = dat)) ## ok
coef(reReg(fm, data = dat, method = "cox.HW")) ## ok
## coef(reReg(fm, data = dat, method = "am.GL"))
coef(reReg(fm, data = dat, method = "am.XCHWY")) ## errors
coef(reReg(fm, data = dat, method = "sc.XCYH")) ## ok



rp <- 100
f1 <- f2 <- f3 <- f4 <- f5 <- f6 <- matrix(NA, rp, 20)
for (i in 1:rp) {
    set.seed(i)
    f1[i,] <- do(200, c(1, 1), c(1, 1), indCen = TRUE)
    set.seed(i)
    f2[i,] <- do(200, c(1, 1), c(1, 1), indCen = FALSE)
    set.seed(i)
    f3[i,] <- do(200, c(1, 1), c(1, 1), indCen = TRUE, type = "am")
    set.seed(i)
    f4[i,] <- do(200, c(1, 1), c(1, 1), indCen = FALSE, type = "am")
    set.seed(i)
    f5[i,] <- do(200, c(1, 1), c(1, 1), indCen = TRUE, type = "sc")
    set.seed(i)
    f6[i,] <- do(200, c(1, 1), c(1, 1), indCen = FALSE, type = "sc")
    ## do(200, -c(1, 1), c(1, 1), indCen = TRUE)
    ## do(200, -c(1, 1), c(1, 1), indCen = FALSE)
    ## do(200, c(-1, 1), c(-1, 1), indCen = TRUE, type = "am")
    ## do(200, c(-1, 1), c(-1, 1), indCen = FALSE, type = "am")
    ## do(200, c(-1, 1), c(1, -1), indCen = TRUE, type = "sc")
    ## do(200, c(-1, 1), c(1, -1), indCen = FALSE, type = "sc")
    print(i)
}

matrix(colMeans(f1), 4)
matrix(colMeans(f2), 4)
matrix(colMeans(f3), 4)
matrix(colMeans(f4), 4)
matrix(colMeans(f5), 4)
matrix(colMeans(f6), 4)

set.seed(15)
do(200, c(1, 1), c(1, 1), indCen = TRUE)

e

## Testing results
## do(200, c(1, 1), c(1, 1), indCen = TRUE)
## do(200, c(1, 1), c(1, 1), indCen = FALSE)
## do(200, c(1, 1), c(1, 1), indCen = TRUE, type = "am")
## do(200, c(1, 1), c(1, 1), indCen = FALSE, type = "am")
## do(200, c(1, 1), c(1, 1), indCen = TRUE, type = "sc")
## do(200, c(1, 1), c(1, 1), indCen = FALSE, type = "sc")
## do(200, -c(1, 1), c(1, 1), indCen = TRUE)
## do(200, -c(1, 1), c(1, 1), indCen = FALSE)
## do(200, c(-1, 1), c(-1, 1), indCen = TRUE, type = "am")
## do(200, c(-1, 1), c(-1, 1), indCen = FALSE, type = "am")
## do(200, c(-1, 1), c(1, -1), indCen = TRUE, type = "sc")
## do(200, c(-1, 1), c(1, -1), indCen = FALSE, type = "sc")
e

#################################################################################

library(tidyverse)
Lam <- function(t, z, exa, exb) z * exb * log(1 + t * exa) / exa / 1
Lam0 <- function(t) log(1 + t) / 1
invLam <- function(t, z, exa, exb) (exp(1 * t * exa / exb / z) - 1) / exa
invHaz <- function(t, z, exa, exb) (exp(4 * t * exa / exb / z) - 1) / exa

simDat <- function(n, a, b, indCen = TRUE, type = c("cox", "am", "sc"), tau = 60, summary = FALSE) {
    type <- match.arg(type)
    if (length(a) != 2L) stop("Require length(a) = 2.")
    if (length(b) != 2L) stop("Require length(b) = 2.")
    simOne <- function(id) {
        z <- ifelse(indCen, 1, rgamma(1, 4, 4))
        x <- c(sample(0:1, 1), rnorm(1))
        exa <- c(exp(x %*% a))
        exb <- c(exp(x %*% b))
        cen <- rexp(1, z * exa / 60)
        if (type == "cox") D <- invHaz(rexp(1), z, 1, exb)
        if (type == "am") D <- invHaz(rexp(1), z, exb, exb)
        if (type == "sc") D <- cen
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
        base <- dat %>% group_by(id) %>% summarize(m = unique(m), d = max(status), x1 = unique(x1), x2 = unique(x2))
        d <- base$d
        x1 <- base$x1
        print(summary(base$m))
        cat("\n")
        cat(paste("Number of failures: ", sum(d), " (", 100 * sum(d) / length(d), "%); ",
                  "Number of censored events: ", sum(d < 1), " (", 100 * sum(d < 1) / length(d), "%)\n",
                  sep = ""))
        cat(paste("Number of x1 == 1: ", sum(x1), " (", 100 * sum(x1) / length(x1), "%); ",
                  "Number of x1 == 0: ", sum(x1 < 1), " (", 100 * sum(x1 < 1) / length(x1), "%)\n", sep = ""))
        cat("\nSummary results for x2:\n")
        print(summary(base$x2))
        cat("\n\n")
    }
    return(dat %>% select(-m, -Z))
}

simDat(1e4, c(1, 1), c(1, 1), indCen = TRUE, summary = TRUE)
simDat(1e4, c(1, 1), c(1, 1), indCen = FALSE, summary = TRUE)
simDat(1e4, -c(1, 1), c(1, 1), indCen = TRUE, summary = TRUE)
simDat(1e4, -c(1, 1), c(1, 1), indCen = FALSE, summary = TRUE)

simDat(1e4, c(1, 1), c(1, 1), indCen = TRUE, summary = TRUE, type = "am")
simDat(1e4, c(1, 1), c(1, 1), indCen = FALSE, summary = TRUE, type = "am")
simDat(1e4, c(-1, 1), c(-1, 1), indCen = TRUE, summary = TRUE, type = "am")
simDat(1e4, c(-1, 1), c(-1, 1), indCen = FALSE, summary = TRUE, type = "am")

simDat(1e4, c(1, 1), c(1, 1), indCen = TRUE, summary = TRUE, type = "sc")
simDat(1e4, c(1, 1), c(1, 1), indCen = FALSE, summary = TRUE, type = "sc")
simDat(1e4, c(-1, 1), c(1, -1), indCen = TRUE, summary = TRUE, type = "sc")
simDat(1e4, c(-1, 1), c(1, -1), indCen = FALSE, summary = TRUE, type = "sc")

dat <- simDat(1e4, c(1, 1), c(1, 1))
table(dat$status)


e

#################################################################################


## -----------------------------------------------------------------------------------
## With readmission data
## -----------------------------------------------------------------------------------
data(readmission, package = "frailtypack")
fm <- reSurv(t.stop, id, event, death) ~ sex + chemo
