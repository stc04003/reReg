set.seed(1)
dat <- simSC(30, c(-1, 1), c(-1, 1))
reObj <- with(dat, reSurv(Time, id, event, status))

## Event plots:
## Default labels
plot(reObj)
plot(reObj, order = FALSE)
## User specified labels
plot(reObj, control = list(xlab = "User xlab", ylab = "User ylab", main = "User title"))

## With hypothetical multiple event types
set.seed(1)
reObj2 <- with(dat, reSurv(Time, id, event * sample(1:3, nrow(dat), TRUE), status))
plot(reObj2)

## CSM plots
plot(reObj, CSM = TRUE)
