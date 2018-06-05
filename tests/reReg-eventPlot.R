library(reReg)
data(readmission)
head(readmission)

## checking reSurv
attach(readmission)
reSurv(t.stop)
reSurv(t.stop, id)
reSurv(t.stop, id, event)
reSurv(t.stop, id, event, death)

reSurv(t.start, t.stop)
reSurv(t.start, t.stop, id)
reSurv(t.start, t.stop, id, event)
reSurv(t.start, t.stop, id, event, death)


reSurv(t.stop)$reDF
reSurv(t.start, t.stop)$reDF
identical(reSurv(t.stop)$reTb, reSurv(t.start, t.stop)$reTb) # FALSE

identical(reSurv(t.stop, id)$reTb, reSurv(t.start, t.stop, id)$reTb) # TRUE
identical(reSurv(t.stop, id, event)$reTb, reSurv(t.start, t.stop, id, event)$reTb) # TRUE
identical(reSurv(t.stop, id, event, death)$reTb, reSurv(t.start, t.stop, id, event, death)$reTb) # TRUE
detach(readmission)


## plotEvents examples

reObj <- with(readmission, reSurv(t.stop, id, event, death))
plot(reObj)

plotEvents(reSurv(t.stop, id, event, death) ~ 1, data = readmission)
plotEvents(reSurv(t.stop, id, event, death) ~ sex, data = readmission)
plotEvents(reSurv(t.stop, id, event, death) ~ sex + chemo, data = readmission)
plotEvents(reSurv(t.stop, id, event * sample(1:3, 861, TRUE), death) ~
               sex + chemo, data = readmission)

## plotCMF examples

plot(reObj, CMF = TRUE)

plotCMF(reObj)
plotCMF(reObj ~ 1, data = readmission)
plotCMF(reObj ~ sex, data = readmission)
plotCMF(reObj ~ sex + chemo, data = readmission)


plotCMF(reSurv(t.stop, id, event, death) ~ 1, data = readmission)
plotCMF(reSurv(t.stop, id, event, death) ~ sex, data = readmission)
plotCMF(reSurv(t.stop, id, event, death) ~ sex + chemo, data = readmission)
plotCMF(reSurv(t.stop, id, event * sample(1:3, 861, TRUE), death) ~
            sex + chemo, data = readmission)
plotCMF(reSurv(t.stop, id, event * sample(1:3, 861, TRUE), death) ~
            sex + chemo, data = readmission,
        control = list(recurrent.type = letters[1:3]))

plotCMF(reSurv(t.stop, id, event, death) ~ sex + chemo, data = readmission,
        control = list(title = "Some title"))

## Same plot with collapse = TRUE
plotCMF(reSurv(t.stop, id, event, death) ~ sex, data = readmission, onePanel = TRUE)
plotCMF(reSurv(t.stop, id, event * sample(1:3, 861, TRUE), death) ~ sex,
        data = readmission, onePanel = TRUE)
plotCMF(reSurv(t.stop, id, event, death) ~ sex + chemo, data = readmission, onePanel = TRUE)

## with data from reda
library(reda)
library(gridExtra)
data(simuDat)

simuMcf <- mcf(Survr(ID, time, event) ~ group + gender, data = simuDat)


plotCMF(reSurv(time, ID, event) ~ group + gender, data = simuDat, onePanel=T)
plot(simuMcf, lty = 1 : 4, legendName = "Treatment & Gender")


plot1 <- plotCMF(reSurv(time, ID, event) ~ group + gender, data = simuDat, onePanel=T)
plot2 <- plot(simuMcf, legendName = "Treatment & Gender")
grid.arrange(plot1, plot2, ncol = 1)

simuMcf <- mcf(Survr(ID, time, event) ~ 1, data = simuDat)

plot1 <- plotCMF(reSurv(time, ID, event) ~ 1, data = simuDat)
plot2 <- plot(simuMcf, legendName = "Treatment & Gender")
grid.arrange(plot1, plot2, ncol = 1)


simuDat0 <- subset(simuDat, ID <= 4)
simuDat0 <- subset(simuDat, ID <= 5)

simuMcf <- mcf(Survr(ID, time, event) ~ 1, data = simuDat0)
p1 <- plotCMF(reSurv(time, ID, event) ~ 1, data = simuDat0)
p2 <- plot(simuMcf)
grid.arrange(p1, p2, ncol = 1)

subset(simuDat0, event == 1)
dim(subset(simuDat0, event == 1))

library(survival)
data(rats)
