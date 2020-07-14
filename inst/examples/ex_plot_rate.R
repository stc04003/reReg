data(simDat)
fm <- Recur(Time, id, event, status) ~ x1 + x2

fit <- reReg(fm, data = simDat, method = "cox")
## Plot both the baseline cumulative rate and hazard function
plot(fit)
## Plot baseline cumulative rate function
plotRate(fit)
plotRate(fit, smooth = TRUE)
## Plot with user-specified labels
plotRate(fit, xlab = "User xlab", ylab = "User ylab", main = "User title")
plotRate(fit, control = list(xlab = "User xlab", ylab = "User ylab", main = "User title"))
