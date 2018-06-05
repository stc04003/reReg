library(reReg)

## plotEvents examples
data(readmission)

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
data(simuDat)

simuMcf <- mcf(Survr(ID, time, event) ~ group + gender, data = simuDat)


plotCMF(reSurv(time, ID, event) ~ group + gender, data = simuDat, onePanel=T)
plot(simuMcf, lty = 1 : 4, legendName = "Treatment & Gender")

require(gridExtra)
plot1 <- plotCMF(reSurv(time, ID, event) ~ group + gender, data = simuDat, onePanel=T)
plot2 <- plot(simuMcf, lty = 1 : 4, legendName = "Treatment & Gender")
grid.arrange(plot1, plot2, ncol = 1)




simuMcf <- mcf(Survr(ID, time, event) ~ 1, data = simuDat)

require(gridExtra)
plot1 <- plotCMF(reSurv(time, ID, event) ~ 1, data = simuDat)
plot2 <- plot(simuMcf, legendName = "Treatment & Gender")
grid.arrange(plot1, plot2, ncol = 1)


simuDat0 <- subset(simuDat, ID <= 4)
plotCMF(reSurv(time, ID, event) ~ 1, data = simuDat0, onePanel=T)
plotCMF(reSurv(time, ID, event) ~ group + gender, data = simuDat0, onePanel=T)
