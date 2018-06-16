#############################################################################################
## load packages and data readmission
#############################################################################################

library(reReg)
library(reda)
library(gridExtra)
data(readmission, package = "frailtypack")

## ------------------------------------------------------------------------------------------
## checking reSurv
## ------------------------------------------------------------------------------------------
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

## ------------------------------------------------------------------------------------------
## checking plotEvents
## ------------------------------------------------------------------------------------------
reObj <- with(readmission, reSurv(t.stop, id, event, death))

plot(reObj)
plot(reObj, xlab = "User X", ylab = "User Y", title = "User title")
plot(reObj, control = list(xlab = "User X", ylab = "User Y", title = "User title"))

plot(reObj, order = FALSE)
plot(reObj, order = FALSE, xlab = "User X", ylab = "User Y", title = "User title")
plot(reObj, order = FALSE, control = list(xlab = "User X", ylab = "User Y", title = "User title"))

plotEvents(reSurv(t.stop, id, event, death) ~ 1, data = readmission)
plotEvents(reSurv(t.stop, id, event, death) ~ 1, data = readmission,
           xlab = "User X", ylab = "User Y", title = "User title")
           
plotEvents(reSurv(t.stop, id, event, death) ~ sex, data = readmission)
plotEvents(reSurv(t.stop, id, event, death) ~ sex, data = readmission,
           xlab = "User X", ylab = "User Y", title = "User title")           
plotEvents(reSurv(t.stop, id, event, death) ~ sex + chemo, data = readmission)
plotEvents(reSurv(t.stop, id, event, death) ~ sex + chemo, data = readmission,
           xlab = "User X", ylab = "User Y", title = "User title",
           control = list(terminal.name = "User terminal",
                          recurrent.name = "User event"))

plot(reObj, CSM = TRUE)
plot(reObj, CSM = TRUE, xlab = "User X", ylab = "User Y", title = "User title")
plot(reObj, CSM = TRUE, control = list(xlab = "User X", ylab = "User Y", title = "User title"))


## multiple event types
reObj <- with(readmission, reSurv(t.stop, id, event * sample(1:3, 861, TRUE), death))
plot(reObj)
plot(reObj, xlab = "User X", ylab = "User Y", title = "User title")
plot(reObj, control = list(xlab = "User X", ylab = "User Y", title = "User title"))
plot(reObj, xlab = "User X", ylab = "User Y", title = "User title",
     control = list(recurrent.name = "User event"))
plot(reObj, xlab = "User X", ylab = "User Y", title = "User title",
     control = list(recurrent.type = letters[1:3]))

plotEvents(reSurv(t.stop, id, event, death) ~ 1, data = readmission)
plotEvents(reSurv(t.stop, id, event, death) ~ sex, data = readmission,
           control = list(xlab = "User X", ylab = "User Y", title = "User title"))
plotEvents(reSurv(t.stop, id, event, death) ~ sex + chemo, data = readmission,
           xlab = "User X", ylab = "User Y", title = "User title")

set.seed(123)
fm <- reSurv(t.stop, id, event * sample(1:3, 861, TRUE), death) ~ sex + chemo
plotEvents(fm, data = readmission)
plotEvents(fm, data = readmission, xlab = "User X", ylab = "User Y", title = "User title",
           control = list(recurrent.type = letters[1:3]))

plotCSM(fm, data = readmission)
plotCSM(fm, data = readmission, xlab = "User X", ylab = "User Y", title = "User title",
           control = list(recurrent.type = letters[1:3]))
plotCSM(fm, data = readmission, xlab = "User X", ylab = "User Y", title = "User title",
           control = list(recurrent.name = "Types", recurrent.type = letters[1:3]))
plotCSM(fm, data = readmission, xlab = "User X", ylab = "User Y", title = "User title",
           control = list(recurrent.type = letters[1:3])) + ggplot2::theme(legend.position = "bottom")

## ------------------------------------------------------------------------------------------
## checking simulated data generator
## ------------------------------------------------------------------------------------------

simDat(5e3, c(1, -1), c(1, -1), indCen = TRUE, summary = TRUE)
simDat(5e3, c(1, -1), c(1, -1), indCen = FALSE, summary = TRUE)
simDat(5e3, c(1, -1), c(1, -1), indCen = TRUE, summary = TRUE, type = "am")
simDat(5e3, c(1, -1), c(1, -1), indCen = FALSE, summary = TRUE, type = "am")

simDat(5e3, c(1, -1), c(1, -1), indCen = TRUE, summary = TRUE, type = "sc")
simDat(5e3, c(1, -1), c(1, -1), indCen = FALSE, summary = TRUE, type = "sc")
simDat(5e3, c(0, 0), c(1, -1), indCen = TRUE, summary = TRUE, type = "sc")
simDat(5e3, c(0, 0), c(1, -1), indCen = FALSE, summary = TRUE, type = "sc")
simDat(5e3, c(1, -1), c(0, 0), indCen = TRUE, summary = TRUE, type = "sc")
simDat(5e3, c(1, -1), c(0, 0), indCen = FALSE, summary = TRUE, type = "sc")
simDat(5e3, c(-1, -1), c(1, 1), indCen = TRUE, summary = TRUE, type = "sc")
simDat(5e3, c(-1, -1), c(1, 1), indCen = FALSE, summary = TRUE, type = "sc")
e

## ------------------------------------------------------------------------------------------
## checking point esitmation
## ------------------------------------------------------------------------------------------

fm <- reSurv(Time, id, event, status) ~ x1 + x2
B <- 100

## Under Cox model; independent censoring
f1 <- f2 <- f3 <- f4 <- f5 <- matrix(NA, B, 4)
for (i in 1:B) {
    set.seed(i)
    dat <- simDat(200, a = c(1, -1), b = c(1, -1), type = "cox", indCen = TRUE)
    f1[i,] <- coef(reReg(fm, data = dat))
    f2[i,] <- coef(reReg(fm, data = dat, method = "cox.HW"))
    f3[i,] <- coef(reReg(fm, data = dat, method = "am.GL"))
    f4[i,] <- coef(reReg(fm, data = dat, method = "am.XCHWY"))
    f5[i,] <- coef(reReg(fm, data = dat, method = "sc.XCYH"))
    if (i %% 20 == 0) print(i)
}

sapply(1:5, function(x) eval(parse(text = paste("matrix(colMeans(f", x, "), 4)", sep = ""))))
sapply(1:5, function(x) eval(parse(text = paste("matrix(apply(f", x, ", 2, median), 4)", sep = ""))))

## > sapply(1:5, function(x) eval(parse(text = paste("matrix(colMeans(f", x, "), 4)", sep = ""))))
##            [,1]       [,2]       [,3]      [,4]        [,5]
## [1,]  1.0045515  0.9866267 -0.6408149  1.693991 -0.06327675
## [2,] -0.9989104 -0.9846831 33.9799462 -1.702738  0.02135911
## [3,]  0.0000000  0.9831873  1.9082346  1.697227  0.98241508
## [4,]  0.0000000 -0.9866988 -1.8796832 -1.762501 -0.99980678
## > sapply(1:5, function(x) eval(parse(text = paste("matrix(apply(f", x, ", 2, median), 4)", sep = ""))))
##            [,1]       [,2]       [,3]      [,4]         [,5]
## [1,]  1.0089382  0.9811865 -0.7239646  1.719878 -0.064861179
## [2,] -0.9999909 -0.9734948  1.9427521 -1.700137 -0.003709682
## [3,]  0.0000000  0.9921877  1.9117672  1.674758  0.958277053
## [4,]  0.0000000 -0.9835988 -1.8527776 -1.658847 -1.013114033


## Under Cox model: informative censoring
f1 <- f2 <- f3 <- f4 <- f5 <- matrix(NA, B, 4)
for (i in 1:B) {
    set.seed(i)
    dat <- simDat(200, a = c(1, -1), b = c(1, -1), type = "cox", indCen = FALSE)
    f1[i,] <- coef(reReg(fm, data = dat))
    f2[i,] <- coef(reReg(fm, data = dat, method = "cox.HW"))
    f3[i,] <- coef(reReg(fm, data = dat, method = "am.GL"))
    f4[i,] <- coef(reReg(fm, data = dat, method = "am.XCHWY"))
    f5[i,] <- coef(reReg(fm, data = dat, method = "sc.XCYH"))
    if (i %% 20 == 0) print(i)
}

sapply(1:5, function(x) eval(parse(text = paste("matrix(colMeans(f", x, "), 4)", sep = ""))))
sapply(1:5, function(x) eval(parse(text = paste("matrix(apply(f", x, ", 2, median), 4)", sep = ""))))
