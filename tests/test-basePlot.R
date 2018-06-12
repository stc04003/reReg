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
## lines(t0, colMeans(base3)[1:length(t0)], col = 2) 
## lines(t0, apply(base3, 2, median)[1:length(t0)], col = 2)
## lines(t0, colMeans(base3)[1:length(t0) + 1 * length(t0)], col = 3) ## okay

## plot(t0, H0(t0), 'l')
## lines(t0, colMeans(base3)[1:length(t0) + 2 * length(t0)], col = 2)
## lines(t0, colMeans(base3)[1:length(t0) + 3 * length(t0) ], col = 3)  ## okay

## plot(t0, R0(t0), 'l')
## lines(t0, colMeans(base4)[1:length(t0)], col = 2) 
## lines(t0, colMeans(base4)[1:length(t0) + 1 * length(t0)], col = 3) ## okay

## plot(t0, H0(t0), 'l')
## lines(t0, colMeans(base4)[1:length(t0) + 2 * length(t0)], col = 2)
## lines(t0, colMeans(base4)[1:length(t0) + 3 * length(t0) ], col = 3)  ## okay

## foo <- replicate(100, coef(reReg(fm, data = simDat(200, c(-1, 1), c(-1, 1), type = "am", indCen = TRUE), method = "am.GL")))
## rowMeans(foo)

set.seed(59)
dat <- simDat(200, c(-1, 1), c(-1, 1), type = "am", indCen = TRUE)
fit1 <- reReg(fm, data = dat, method = "am.GL")
coef(fit1) ## [1] -0.9743668  0.9328660 -1.0998299  0.8394684

fit1$rate0(1:60)
##  [1] 2.375244 3.645588 4.394820 4.676356 4.956797 5.140377 5.523711 5.796210
##  [9] 6.116759 6.314676 6.332858 6.372134 6.455467 6.596160 6.767956 6.826780
## [17] 7.157439 7.218045 7.318045 7.485862 7.524323 7.524323 7.574323 7.574323
## [25] 7.928490 8.061823 8.138746 8.292592 8.292592 8.292592 8.565320 8.838047
## [33] 8.838047 8.838047 8.838047 8.838047 9.038047 9.038047 9.038047 9.038047
## [41] 9.038047 9.038047 9.038047 9.038047 9.180904 9.180904 9.180904 9.180904
## [49] 9.180904 9.380904 9.380904 9.380904 9.380904 9.380904 9.380904 9.380904
## [57] 9.380904 9.380904 9.380904 9.380904

t0 <- seq(0, 60, .1)
plot(t0, R0(t0), 's')
lines(t0, fit1$rate0(t0), 's', col = 2)

plot(t0, H0(t0), 's')
lines(t0, fit1$haz0(t0), 's', col = 2)

plot(t0, R0(t0) / max(R0(t0)), 's')
lines(t0, fit1$rate0(t0) / max(fit1$rate0(t0)), 's', col = 2)

plot(t0, H0(t0) / max(H0(t0)), 's')
lines(t0, fit1$haz0(t0) / max(fit1$haz0(t0)), 's', col = 2)


##

baseId <- as.numeric(cumsum(unlist(lapply(split(dat$id, dat$id), length))))
X0 <- matrix(c(dat$x1[baseId], dat$x2[baseId]), length(unique(dat$id)))
Y <- dat$Time[baseId]

summary(Y * exp(X0 %*% coef(fit1)[3:4]))
summary(fit1$t0.haz)

summary(with(dat, Time * exp(cbind(x1, x2) %*% coef(fit1)[1:2])))

summary(fit1$t0.rate)


summary(X0 %*% (coef(fit1)[1:2] - coef(fit1)[3:4]))
summary(X0 %*% -(coef(fit1)[1:2] - coef(fit1)[3:4]))
