data(simDat)
reObj <- with(simDat, Recur(Time, id, event, status))

## Event plots:
plot(reObj)
plot(reObj, event.result = "decreasing")
plot(reObj, control = list(xlab = "User xlab", ylab = "User ylab", main = "User title"))

## MCF plots
plot(reObj, mcf = TRUE)
plot(reObj, mcf = TRUE, mcf.adjrisk = FALSE)
plot(reObj, mcf = TRUE, mcf.smooth = TRUE)
plot(reObj, mcf = TRUE, control = list(xlab = "User xlab", ylab = "User ylab", main = "User title"))

## With (hypothetical) multiple event types
set.seed(1)
reObj2 <- with(simDat, Recur(Time, id, event * sample(1:3, nrow(simDat), TRUE), status))
plot(reObj2)
plot(reObj2, mcf = TRUE)

