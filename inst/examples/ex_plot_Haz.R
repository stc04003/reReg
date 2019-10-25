set.seed(1)
dat <- simSC(50, c(-1, 1), c(-1, 1))
fm <- Recur(Time, id, event, status) ~ x1 + x2

fit <- reReg(fm, data = dat, method = "cox.HW")
## Plot both the baseline cumulative rate and hazard function
plot(fit)
## Plot baseline cumulative hazard function
plotHaz(fit)
plotHaz(fit, smooth = TRUE)
## Plot with user-specified labels
plotHaz(fit, xlab = "User xlab", ylab = "User ylab", main = "User title")
plotHaz(fit, control = list(xlab = "User xlab", ylab = "User ylab", main = "User title"))

## With 95% confidence interval when `se` is enabled
fit <- reReg(fm, data = dat, method = "cox.HW", se = "resampling", B = 20)
plotHaz(fit)
