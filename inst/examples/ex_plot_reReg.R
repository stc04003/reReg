data(simDat)
fm <- Recur(t.start %to% t.stop, id, event, status) ~ x1 + x2

fit <- reReg(fm, data = simDat, B = 0)
plot(fit)
plot(fit, xlab = "Time (days)", smooth = TRUE)

## Predicted cumulative rate and hazard given covariates
newdata <- expand.grid(x1 = 0:1, x2 = mean(simDat$x2))
plot(fit, newdata = newdata, showName = TRUE)
