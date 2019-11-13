library(reReg)
library(ggplot2)
library(gridExtra)

packageVersion("reReg")

data(readmission, package = "frailtypack")
head(readmission)
readmission <- subset(readmission, !(id %in% c(60, 109, 280)))
attach(readmission)

Recur(t.start %to% t.stop, id = id, event = event, status = death)
Recur(time = t.stop, id = id, event = event, status = death)

Recur(t.start %to% t.stop, id, event, death)
Recur(t.stop, id, event, death)
Recur(t.stop, id, event, death)

reObj <- Recur(t.stop, id, event, death)
plot(reObj)

plot(reObj, cex = 1.5, xlab = "Time in days", ylab = "Patients", 
     main = "Event plot for readmission data", 
     terminal.name = "Death", 
     recurrent.name = "Hospital readmission")

detach(readmission)

plotEvents(reObj)
plotEvents(reObj, data = readmission)
plotEvents(reObj ~ 1, data = readmission)
plotEvents(Recur(t.stop, id, event, death) ~ 1, data = readmission)

plotEvents(Recur(t.stop, id, event, death) ~ 1, data = readmission,
           cex = 1.5, xlab = "Time in days", ylab = "Patients", 
           main = "Event plot for readmission data", 
           terminal.name = "Death", recurrent.name = "Hospital readmission")

plotEvents(Recur(t.stop, id, event, death) ~ sex, data = readmission)
plotEvents(Recur(t.stop, id, event, death) ~ sex + chemo, data = readmission)
plot(reObj, CSM = TRUE)

p1 <- plotCSM(Recur(t.stop, id, event, death) ~ 1, data = readmission, main = "")
p2 <- plotCSM(Recur(t.stop, id, event, death) ~ 1, data = readmission, 
              adjrisk = FALSE, main = "")
grid.arrange(p1, p2, ncol=2)

p1 <- plotCSM(Recur(t.stop, id, event, death) ~ sex + chemo, data = readmission, main = "")
p2 <- plotCSM(Recur(t.stop, id, event, death) ~ sex + chemo, data = readmission, 
              adjrisk = FALSE, main = "")
grid.arrange(p1, p2, ncol=2)

set.seed(1)
readmission$event2 <- readmission$event * sample(1:3, 852, TRUE)
plotEvents(Recur(t.stop, id, event2, death) ~ sex, data = readmission)

plotCSM(Recur(t.stop, id, event2, death) ~ sex, adjrisk = FALSE, data = readmission,
        recurrent.name = "Event types", recurrent.type = c("Type 1", "Type 2", "Type 3"))


p1 <- plotCSM(Recur(t.stop, id, event2, death) ~ sex, data = readmission, main = "") +
  theme(legend.position="none")
p2 <- plotCSM(Recur(t.stop, id, event2, death) ~ sex, data = readmission, 
              adjrisk = FALSE, main = "") +
  theme(legend.position="none")
grid.arrange(p1, p2, ncol = 2)
