set.seed(1)
dat <- simSC(200, c(-1, 1), c(-1, 1))

## being deprecated in Verson 1.1.7
with(dat, reSurv(Time, id, event, status))
## Use Recur() instead
with(dat, Recur(Time, id, event, status))
