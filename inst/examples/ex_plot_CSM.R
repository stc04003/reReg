set.seed(1)
dat <- simSC(30, c(-1, 1), c(-1, 1))
plotCSM(Recur(Time, id, event, status) ~ 1, data = dat)
plotCSM(Recur(Time, id, event, status) ~ x1, data = dat)
plotCSM(Recur(Time, id, event, status) ~ x1, data = dat, onePanel = TRUE)


