data(simDat)
fm <- Recur(t.stop, id, event, status) ~ x1 + x2

fit <- reReg(fm, data = simDat, method = "cox", B = 0)
plot(fit)
plot(fit, baseline = "rate")
plot(fit, baseline = "rate", xlab = "Time (days)", smooth = TRUE)

