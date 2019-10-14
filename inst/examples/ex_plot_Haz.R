set.seed(1)
dat <- simSC(50, c(-1, 1), c(-1, 1))
fit <- reReg(reSurv(Time, id, event, status) ~ x1 + x2, data = dat,
             method = "cox.HW", se = "resampling", B = 20)
## Plot both the baseline cumulative rate and hazard function
plot(fit)
## Plot baseline cumulative hazard function
plotHaz(fit)
## Plot with user-specified labels
plotHaz(fit, control = list(xlab = "User xlab", ylab = "User ylab", main = "User title"))  
