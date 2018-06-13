library(reReg)

fm <- reSurv(Time, id, event, status) ~ x1 + x2

## Use simulated data
R0 <- function(x) log(1 + x)
H0 <- function(x) log(1 + x) / 4

## Cox type
dat <- simDat(200, c(1, 1), c(1, 1), type = "cox", indCen = TRUE)
fit <- reReg(fm, data = dat)
summary(fit)

plot(fit)
plotRate(fit)
plotHaz(fit)

fit$rate0(50:70)
fit$haz0(50:70)

fit2 <- reReg(fm, data = dat, method = "cox.HW")
plot(fit2)


t0 <- seq(0, 100, .01)[-1]

do <- function(indCen = TRUE) {
    dat <- simDat(200, c(1, 1), c(1, 1), type = "cox", indCen = indCen)
    fit1 <- reReg(fm, data = dat)
    fit2 <- reReg(fm, data = dat, method = "cox.HW")
    c(fit1$rate0(t0), fit2$rate0(t0), fit1$haz0(t0), fit2$haz0(t0))
}

## base1 <- base2 <- matrix(NA, 500, 4 * length(t0))
## for (i in 1:500) {
##     base1[i,] <- do(TRUE)
##     base2[i,] <- do(FALSE)
##     if (i %% 5 == 0) print(i)
## }

## plot(t0, R0(t0), 'l')
## lines(t0, colMeans(base1)[1:length(t0)], col = 2) ## okay
## lines(t0, colMeans(base1)[1:length(t0) + 1 * length(t0)], col = 3) ## okay

## plot(t0, H0(t0), 'l')
## lines(t0, colMeans(base1)[1:length(t0) + 2 * length(t0)], col = 2)  ## okay
## lines(t0, colMeans(base1)[1:length(t0) + 3 * length(t0) ], col = 3)  ## okay

## plot(t0, R0(t0), 'l')
## lines(t0, colMeans(base2)[1:length(t0)], col = 2) ## okay
## lines(t0, colMeans(base2)[1:length(t0) + 1 * length(t0)], col = 3) ## okay

## plot(t0, H0(t0), 'l')
## lines(t0, colMeans(base2)[1:length(t0) + 2 * length(t0)], col = 2)  ## okay
## lines(t0, colMeans(base2)[1:length(t0) + 3 * length(t0) ], col = 3)  ## okay

do2 <- function(indCen = TRUE) {
    dat <- simDat(200, c(-1, 1), c(-1, 1), type = "am", indCen = indCen)
    fit1 <- reReg(fm, data = dat, method = "am.GL")
    fit2 <- reReg(fm, data = dat, method = "am.XCHWY")
    c(fit1$rate0(t0), fit2$rate0(t0), fit1$haz0(t0), fit2$haz0(t0))
}

## base3 <- base4 <- matrix(NA, 100, 4 * length(t0))
## for (i in 1:100) {
##     base3[i,] <- do2(TRUE)
##     base4[i,] <- do2(FALSE)
##     if (i %% 5 == 0) print(i)
## }

## plot(t0, R0(t0), 'l')
## lines(t0, colMeans(base3)[1:length(t0) + 1 * length(t0)], col = 3) ## okay
## lines(t0, colMeans(base3)[1:length(t0)], col = 2) 
## lines(t0, apply(base3, 2, median)[1:length(t0)], col = 2)

## plot(t0, H0(t0), 'l')
## lines(t0, colMeans(base3)[1:length(t0) + 2 * length(t0)], col = 2)
## lines(t0, colMeans(base3)[1:length(t0) + 3 * length(t0) ], col = 3)  ## okay

## plot(t0, R0(t0), 'l')
## lines(t0, colMeans(base4)[1:length(t0) + 1 * length(t0)], col = 3) ## okay
## lines(t0, colMeans(base4)[1:length(t0)], col = 2) 
## lines(t0, apply(base3, 2, median)[1:length(t0)], col = 2)

## plot(t0, H0(t0), 'l')
## lines(t0, colMeans(base4)[1:length(t0) + 2 * length(t0)], col = 2)
## lines(t0, colMeans(base4)[1:length(t0) + 3 * length(t0) ], col = 3)  ## okay

dat <- simDat(200, c(-1, 1), c(-1, 1), type = "am", indCen = TRUE)
fit <- reReg(fm, data = dat)
plot(fit)
plot(fit, baseline = "rate")
plot(fit, baseline = "haz")

fit <- reReg(fm, data = dat, method = "sc.XCYH")
summary(fit)
plot(fit)
plot(fit, baseline = "rate")
plot(fit, baseline = "haz")

fit <- reReg(fm, data = dat, method = "sc.XCYH", plot.ci = TRUE)
summary(fit)
plot(fit)
plot(fit, baseline = "rate")
plot(fit, baseline = "haz")

plotRate(fit)
plotHaz(fit)

names(fit)
t0 <- seq(0, 80, .01)

plot(t0, fit$rate0(t0) * fit$log.muZ, 's')
## lines(t0, fit$rate0(t0) * exp(fit$log.muZ), 's')
lines(t0, fit$lam0(t0), 's', col = 2)
lines(t0, H0(t0), 's', col = 3)
fit$log.muZ
