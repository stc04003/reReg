data(simDat)
fm <- Recur(t.start %to% t.stop, id, event, status) ~ x1 + x2
fit <- reReg(fm, data = simDat, model = "cox|cox", B = 0)
## Plot both the baseline cumulative rate and hazard function
plot(fit)
## Plot baseline cumulative rate function
plotRate(fit)
plotRate(fit, smooth = TRUE)
