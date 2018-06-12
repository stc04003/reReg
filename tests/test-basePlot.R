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

fit2 <- reReg(fm, data = dat, method = "cox.HW")
plot(fit2)
