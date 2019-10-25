set.seed(1)
dat <- simSC(30, c(-1, 1), c(-1, 1))
plotEvents(Recur(Time, id, event, status) ~ 1, data = dat)
plotEvents(Recur(Time, id, event, status) ~ 1, data = dat,
           xlab = "Time in days", ylab = "Subjects arranged by terminal time")

## Separate plots by x1
plotEvents(Recur(Time, id, event, status) ~ x1, data = dat)

## For multiple recurrent events
dat$x3 <- ifelse(dat$x2 < 0, "x2 < 0", "x2 > 0")
plotEvents(Recur(Time, id, event, status) ~ x1 + x3, data = dat)
plotEvents(Recur(Time, id, event * sample(1:3, nrow(dat), TRUE), status) ~ x1, data = dat)
plotEvents(Recur(Time, id, event * sample(1:3, nrow(dat), TRUE), status) ~ x1 + x3, data = dat)
