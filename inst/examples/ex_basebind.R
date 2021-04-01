data(simDat)
fm <- Recur(t.stop, id, event, status) ~ x1 + x2
fit1 <- reReg(fm, subset = x1 == 0, data = simDat, B = 200)
fit2 <- reReg(fm, subset = x1 == 1, data = simDat, B = 200)
basebind(plot(fit1), plot(fit2))
