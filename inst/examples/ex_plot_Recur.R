data(simDat)
reObj <- with(simDat, Recur(t.start %to% t.stop, id, event, status))

## Event plots:
plot(reObj)
plot(reObj, event.result = "decreasing")

## With (hypothetical) multiple event types
simDat$event2 <- with(simDat, ifelse(t.stop > 10 & event > 0, 2, event))
reObj2 <- with(simDat, Recur(t.start %to% t.stop, id, event2, status))
plot(reObj2)

## With (hypothetical) calendar times
simDat2 <- simDat
simDat2$t.start <- as.Date(simDat2$t.start + simDat2$x2 * 5, origin = "20-01-01")
simDat2$t.stop <- as.Date(simDat2$t.stop + simDat2$x2 * 5, origin = "20-01-01")
reObj3 <- with(simDat2, Recur(t.start %to% t.stop, id, event, status))
plot(reObj3, event.calendarTime = TRUE)

## MCF plots
plot(reObj, mcf = TRUE)
plot(reObj, mcf = TRUE, mcf.adjustRiskset = FALSE)

