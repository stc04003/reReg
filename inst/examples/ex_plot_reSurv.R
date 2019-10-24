set.seed(1)
dat <- simSC(30, c(-1, 1), c(-1, 1))
reObj <- with(dat, Recur(Time, id, event, status))

## Event plots:
plot(reObj)
plot(reObj, event.result = "decreasing")
plot(reObj, event.result = "asis")
plot(reObj, control = list(xlab = "User xlab", ylab = "User ylab", main = "User title"))

## CSM plots
plot(reObj, CSM = TRUE)
plot(reObj, CSM = TRUE, csm.adjrisk = FALSE)
plot(reObj, CSM = TRUE, csm.smooth = TRUE)
plot(reObj, CSM = TRUE, control = list(xlab = "User xlab", ylab = "User ylab", main = "User title"))

## With (hypothetical) multiple event types
set.seed(1)
reObj2 <- with(dat, Recur(Time, id, event * sample(1:3, nrow(dat), TRUE), status))
plot(reObj2)
plot(reObj2, event.result = "decreasing")
plot(reObj2, event.result = "asis")

plot(reObj2, CSM = TRUE)
plot(reObj2, CSM = TRUE, csm.adjrisk = FALSE)
plot(reObj2, CSM = TRUE, csm.smooth = TRUE)

