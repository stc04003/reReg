data(simDat)
plotMCF(Recur(t.stop, id, event, status) ~ 1, data = simDat)
plotMCF(Recur(t.stop, id, event, status) ~ x1, data = simDat, onePanel = TRUE)
