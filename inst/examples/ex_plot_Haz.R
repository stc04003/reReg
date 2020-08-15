data(simDat)
fm <- Recur(t.stop, id, event, status) ~ x1 + x2

fit <- reReg(fm, data = simDat, method = "cox", B = 0)
## Plot both the baseline cumulative rate and hazard function
plot(fit)
## Plot baseline cumulative hazard function
plotHaz(fit)
plotHaz(fit, smooth = TRUE)
