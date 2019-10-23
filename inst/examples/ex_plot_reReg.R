library(reReg)

set.seed(1)
dat <- simSC(50, c(-1, 1), c(-1, 1))

fit <- reReg(reSurv(Time, id, event, status) ~ x1 + x2, data = dat, method = "cox.HW")
plot(fit, baseline = "rate")
plot(fit, baseline = "rate", xlab = "Time (days)")
