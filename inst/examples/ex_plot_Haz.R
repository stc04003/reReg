data(simDat)
fm <- Recur(Time, id, event, status) ~ x1 + x2

fit <- reReg(fm, data = simDat, method = "cox.HW")
## Plot both the baseline cumulative rate and hazard function
plot(fit)
## Plot baseline cumulative hazard function
plotHaz(fit)
plotHaz(fit, smooth = TRUE)
## Plot with user-specified labels
plotHaz(fit, xlab = "User xlab", ylab = "User ylab", main = "User title")
plotHaz(fit, control = list(xlab = "User xlab", ylab = "User ylab", main = "User title"))
