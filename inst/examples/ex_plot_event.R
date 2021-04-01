data(simDat)
plotEvents(Recur(t.start %to% t.stop, id, event, status) ~ 1, data = simDat,
           xlab = "Time in days", ylab = "Subjects arranged by terminal time")

## Separate plots by x1
plotEvents(Recur(t.start %to% t.stop, id, event, status) ~ x1, data = simDat)

## For multiple recurrent events
simDat$x3 <- ifelse(simDat$x2 < 0, "x2 < 0", "x2 > 0")
simDat$event <- simDat$event * sample(1:3, nrow(simDat), TRUE)
plotEvents(Recur(t.start %to% t.stop, id, event, status) ~ x1 + x3, data = simDat)
