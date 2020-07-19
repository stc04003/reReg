data(simDat)
plotMCF(Recur(Time, id, event, status) ~ 1, data = simDat)
plotMCF(Recur(Time, id, event, status) ~ x1, data = simDat, onePanel = TRUE)
