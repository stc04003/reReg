set.seed(1)
dat <- simSC(30, c(-1, 1), c(-1, 1))
plotEvents(reSurv(Time, id, event, status) ~ 1, data = dat)

## Separate plots by x1
plotEvents(reSurv(Time, id, event, status) ~ x1, data = dat)

## Separate plots by x1 and x3
dat$x3 <- ifelse(dat$x2 < 0, "x2 < 0", "x2 > 0")
plotEvents(reSurv(Time, id, event, status) ~ x1 + x3, data = dat)
## With multiple hypothetical event types
plotEvents(reSurv(Time, id, event * sample(1:3, nrow(dat), TRUE), status) ~ x1, data = dat)
