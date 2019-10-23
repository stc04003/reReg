set.seed(1)
dat <- simSC(200, c(-1, 1), c(-1, 1))
with(dat, Recur(Time, id, event, status))
