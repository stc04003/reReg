data(simDat)
m <- mcf(Recur(t.start %to% t.stop, id, event, status) ~ x1, data = simDat)
plot(m, conf.int = TRUE)
