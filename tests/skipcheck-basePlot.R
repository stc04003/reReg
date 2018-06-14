library(reReg)

fm <- reSurv(Time, id, event, status) ~ x1 + x2

## Use simulated data
R0 <- function(x) log(1 + x)
H0 <- function(x) log(1 + x) / 4

## non-parametric
set.seed(1)
dat <- simDat(200, c(1, 1), c(1, 1), type = "cox", indCen = TRUE)
system.time(fit <- reReg(reSurv(Time, id, event, status) ~ x1 + x2,
                         method = "cox.HW", se = "bootstrap", data = dat))
## 87.268   0.000  87.282 
summary(fit)
## Call: reReg(formula = reSurv(Time, id, event, status) ~ x1 + x2, data = dat, 
##     method = "cox.HW", se = "bootstrap")
## Method: Huang-Wang Model 
## Coefficients (rate):
##    Estimate StdErr z.value   p.value    
## x1    1.064  0.204   5.205 < 2.2e-16 ***
## x2    1.114  0.086  12.902 < 2.2e-16 ***
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
## Coefficients (hazard):
##    Estimate StdErr z.value   p.value    
## x1    1.143  0.233   4.906 < 2.2e-16 ***
## x2    1.217  0.128   9.542 < 2.2e-16 ***
## ---
##     Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

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


t0 <- seq(0, 80, .05)[-1]
rateMat <- matrix(NA, 100, length(t0) * 5)
hazMat <- matrix(NA, 100, length(t0) * 4)
for (i in 1:100) {
    set.seed(i)
    dat <- simDat(200, c(-1, 1), c(-1, 1), type = "am", indCen = TRUE)
    fit1 <- reReg(fm, data = dat)
    fit2 <- reReg(fm, data = dat, method = "cox.HW")
    fit3 <- reReg(fm, data = dat, method = "am.GL")
    fit4 <- reReg(fm, data = dat, method = "am.XCHWY")
    fit5 <- reReg(fm, data = dat, method = "sc.XCYH")
    rateMat[i,] <- c(fit1$rate0(t0), fit2$rate0(t0), fit3$rate0(t0), fit4$rate0(t0), fit5$rate0(t0))
    hazMat[i,] <- c(fit1$rate0(t0), fit2$rate0(t0), fit1$haz0(t0), fit2$haz0(t0))
}

plot(t0, R0(t0), 's')
sapply(2:6, function(x)
    invisible(lines(t0, colMeans(rateMat)[1:length(t0) + (x - 2) * length(t0)], 's', col = x)))
sapply(2:6, function(x)
    invisible(lines(t0, apply(rateMat, 2, median)[1:length(t0) + (x - 2) * length(t0)], 's', col = x, lty = 2)))

lines(t0, colMeans(rateMat)[1:length(t0) + 2 * length(t0)], 's', col = 2)
e
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

fit2 <- reReg(fm, data = dat, method = "sc.XCYH", se = "resampling")
summary(fit2)
plot(fit2)
plot(fit2, baseline = "rate")
plot(fit2, baseline = "haz")

plotRate(fit)
plotHaz(fit)

names(fit)
t0 <- seq(0, 80, .01)

plot(t0, R0(t0), 's')
lines(t0, fit$rate0(t0), 's', col = 2)

fit <- reReg(fm, data = dat, method = "sc.XCYH")
debug(doNonpara.sc.XCYH)
doNonpara.sc.XCYH(DF = DF, alpha = fit$alpha, beta = fit$beta, engine = engine, stdErr = stdErr)


do <- function() {
    dat <- simDat(200, c(-1, 1), c(-1, 1), type = "am", indCen = TRUE)
    fit <- reReg(fm, data = dat, method = "sc.XCYH")
    fit$rate0(t0)
}

t0 <- seq(0, 80, .01)
rateMat <- matrix(NA, 100, length(t0))
for (i in 1:100) rateMat[i,] <- do()

plot(t0, R0(t0), 's')
lines(t0, colMeans(rateMat), 's', col = 2)
lines(t0, apply(rateMat, 2, median), 's', col = 3)
lines(t0, fit$rate0(t0), col = 4, 's')
lines(t0, fit$rate0(t0), col = 4, 's')




data(readmission, package = "frailtypack")
set.seed(123)
fit <- reReg(reSurv(t.stop, id, event, death) ~ sex + chemo,
             data = subset(readmission, id < 50),
             method = "am.XCHWY", se = "resampling", B = 20)
fit
summary(fit)

## Plot both the baseline cumulative rate and hazard function
plot(fit)
## Plot baseline cumulative hazard function
plotHaz(fit)
## Plot with user-specified labels
plotHaz(fit, control = list(xlab = "User xlab", ylab = "User ylab", title = "User title"))  
