data(simDat)
plotCSM(Recur(Time, id, event, status) ~ 1, data = simDat)
plotCSM(Recur(Time, id, event, status) ~ x1, data = simDat)
plotCSM(Recur(Time, id, event, status) ~ x1, data = simDat, onePanel = TRUE)
