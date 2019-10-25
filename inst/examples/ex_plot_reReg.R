set.seed(1)
dat <- simSC(50, c(-1, 1), c(-1, 1))
fm <- Recur(Time, id, event, status) ~ x1 + x2

fit <- reReg(fm, data = dat, method = "cox.HW")
plot(fit)
plot(fit, baseline = "rate")
plot(fit, baseline = "rate", xlab = "Time (days)")
plot(fit, baseline = "rate", smooth = TRUE)


fit <- reReg(fm, data = dat, method = "cox.HW", se = "resampling", B = 20)
plot(fit)
