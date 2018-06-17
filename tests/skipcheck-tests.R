#############################################################################################
## load packages and data readmission
#############################################################################################

library(reReg)
library(reda)
library(gridExtra)
data(readmission, package = "frailtypack")

R0 <- function(x) log(1 + x)
H0 <- function(x) log(1 + x) / 4

fm <- reSurv(Time, id, event, status) ~ x1 + x2    

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
simDat(5e3, c(1, 1), c(-1, -1), indCen = TRUE, summary = TRUE, type = "sc")
simDat(5e3, c(1, 1), c(-1, -1), indCen = FALSE, summary = TRUE, type = "sc")

e

## ------------------------------------------------------------------------------------------
## checking point esitmation
## ------------------------------------------------------------------------------------------

fm <- reSurv(Time, id, event, status) ~ x1 + x2
B <- 200

## -----------------------------------------
## Under Cox model; independent censoring
## -----------------------------------------
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

## -----------------------------------------
## Under Cox model: informative censoring
## -----------------------------------------
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

## > sapply(1:5, function(x) eval(parse(text = paste("matrix(colMeans(f", x, "), 4)", sep = ""))))
##            [,1]       [,2]      [,3]      [,4]        [,5]
## [1,]  0.8689491  0.9512708 -2.124821  1.681021 -0.04609916
## [2,] -0.8869352 -0.9785694 -4.683311 -1.689699  0.01385236
## [3,]  0.0000000  0.9672966  1.894967 13.441088  0.94418619
## [4,]  0.0000000 -1.0041152 -1.923131 -1.828242 -0.99715960
## > sapply(1:5, function(x) eval(parse(text = paste("matrix(apply(f", x, ", 2, median), 4)", sep = ""))))
##            [,1]       [,2]       [,3]      [,4]        [,5]
## [1,]  0.8777029  0.9576853  0.2604682  1.680410  0.02014634
## [2,] -0.8885664 -0.9664461 -1.2535982 -1.659235  0.01751235
## [3,]  0.0000000  0.9570289  1.9145469  1.686200  0.93996583
## [4,]  0.0000000 -0.9913744 -1.9149213 -1.718060 -1.00247617

## ----------------------------------------------------
## Under accelerated mean model: independent censoring
## ----------------------------------------------------
f1 <- f2 <- f3 <- f4 <- f5 <- matrix(NA, B, 4)
for (i in 1:B) {
    set.seed(i)
    dat <- simDat(200, a = c(1, -1), b = c(1, -1), type = "am", indCen = TRUE)
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
##            [,1]       [,2]       [,3]       [,4]       [,5]
## [1,]  0.3293329  0.4744924  1.0534743  1.0170697  0.9939016
## [2,] -0.3181031 -0.4684926  2.1080002 -0.9919633 -0.9970334
## [3,]  0.0000000  0.4766892  1.0221018  1.0179075  1.0041706
## [4,]  0.0000000 -0.4704039 -0.9810203 -0.9850167 -0.9899014
## > sapply(1:5, function(x) eval(parse(text = paste("matrix(apply(f", x, ", 2, median), 4)", sep = ""))))
##            [,1]       [,2]       [,3]       [,4]       [,5]
## [1,]  0.3243937  0.4722361  0.9887925  1.0117680  0.9860829
## [2,] -0.3174810 -0.4458255 -0.9756880 -0.9862884 -0.9939761
## [3,]  0.0000000  0.4868080  1.0327419  1.0195803  1.0114809
## [4,]  0.0000000 -0.4447878 -0.9911142 -0.9908719 -0.9898207

## ----------------------------------------------------
## Under accelerated mean model: informative censoring
## ----------------------------------------------------

f1 <- f2 <- f3 <- f4 <- f5 <- matrix(NA, B, 4)
for (i in 1:B) {
    set.seed(i)
    dat <- simDat(200, a = c(1, -1), b = c(1, -1), type = "am", indCen = FALSE)
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
##            [,1]       [,2]        [,3]       [,4]       [,5]
## [1,]  0.2914423  0.4734534   7.4112595  0.9963321  0.9594337
## [2,] -0.2941699 -0.4972407 -46.0428575 -1.0120288 -0.9812084
## [3,]  0.0000000  0.4645129   0.9418970  0.9548308  0.9784432
## [4,]  0.0000000 -0.4924252  -0.9800455 -0.9943737 -0.9857935
## > sapply(1:5, function(x) eval(parse(text = paste("matrix(apply(f", x, ", 2, median), 4)", sep = ""))))
##            [,1]       [,2]       [,3]       [,4]       [,5]
## [1,]  0.2967115  0.4855938  0.9751950  0.9834565  0.9374652
## [2,] -0.2940373 -0.4852154 -0.9806555 -1.0167114 -0.9538417
## [3,]  0.0000000  0.4453158  0.9591747  0.9907931  0.9679560
## [4,]  0.0000000 -0.4776840 -0.9818601 -1.0148674 -0.9870059

## -------------------------------------------------------
## Under scale-change model case 1 : independent censoring
## -------------------------------------------------------

f1 <- f2 <- f3 <- f4 <- f5 <- matrix(NA, B, 4)
for (i in 1:B) {
    set.seed(i)
    dat <- simDat(200, a = c(1, -1), b = c(1, -1), type = "sc", indCen = TRUE)
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
##            [,1]       [,2]       [,3]       [,4]       [,5]
## [1,]  0.3376109  0.4776484  2.1180249  1.0417487  0.9493306
## [2,] -0.3171822 -0.4654617  3.4339798 -0.9979099 -1.0091965
## [3,]  0.0000000  0.4602072  0.9843028  0.9817925  0.9933374
## [4,]  0.0000000 -0.4593015 -0.9580159 -0.9593768 -0.9986040
## > sapply(1:5, function(x) eval(parse(text = paste("matrix(apply(f", x, ", 2, median), 4)", sep = ""))))
##            [,1]       [,2]       [,3]       [,4]       [,5]
## [1,]  0.3320632  0.4847381  1.0028684  1.0546360  0.9194267
## [2,] -0.3176486 -0.4428428 -0.9940145 -1.0031604 -1.0032721
## [3,]  0.0000000  0.4525045  0.9764476  0.9822547  0.9972095
## [4,]  0.0000000 -0.4390052 -0.9422122 -0.9656356 -1.0070125

## -------------------------------------------------------
## Under scale-change model, case 1: informative censoring
## -------------------------------------------------------

f1 <- f2 <- f3 <- f4 <- f5 <- matrix(NA, B, 4)
for (i in 1:B) {
    set.seed(i)
    dat <- simDat(200, a = c(1, -1), b = c(1, -1), type = "sc", indCen = FALSE)
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
##            [,1]       [,2]        [,3]       [,4]       [,5]
## [1,]  0.2914423  0.4734534   7.4112595  0.9963321  0.9594337
## [2,] -0.2941699 -0.4972407 -46.0428575 -1.0120288 -0.9812084
## [3,]  0.0000000  0.4645129   0.9418970  0.9548308  0.9784432
## [4,]  0.0000000 -0.4924252  -0.9800455 -0.9943737 -0.9857935
## > sapply(1:5, function(x) eval(parse(text = paste("matrix(apply(f", x, ", 2, median), 4)", sep = ""))))
##            [,1]       [,2]       [,3]       [,4]       [,5]
## [1,]  0.2967115  0.4855938  0.9751950  0.9834565  0.9374652
## [2,] -0.2940373 -0.4852154 -0.9806555 -1.0167114 -0.9538417
## [3,]  0.0000000  0.4453158  0.9591747  0.9907931  0.9679560
## [4,]  0.0000000 -0.4776840 -0.9818601 -1.0148674 -0.9870059

## -------------------------------------------------------
## Under scale-change model case 2 : independent censoring
## -------------------------------------------------------

f1 <- f2 <- f3 <- f4 <- f5 <- matrix(NA, B, 4)
for (i in 1:B) {
    set.seed(i)
    dat <- simDat(200, a = c(0, 0), b = c(1, -1), type = "sc", indCen = TRUE)
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
## [1,]  1.0095395  1.0030558 -12.258359  1.701714 -0.07104285
## [2,] -0.9987734 -0.9828919  70.284788 -1.694750  0.00794948
## [3,]  0.0000000  0.9925946   1.901420  1.693173  0.98412444
## [4,]  0.0000000 -0.9903248  -1.889879 -1.793086 -0.99590091
## > sapply(1:5, function(x) eval(parse(text = paste("matrix(apply(f", x, ", 2, median), 4)", sep = ""))))
##            [,1]       [,2]      [,3]      [,4]        [,5]
## [1,]  1.0137409  0.9860589 -0.569459  1.713592 -0.07785388
## [2,] -0.9986787 -0.9726049  1.349517 -1.694036 -0.02879766
## [3,]  0.0000000  0.9926598  1.898220  1.674650  0.96859297
## [4,]  0.0000000 -0.9838892 -1.870775 -1.676171 -1.01655726

## -------------------------------------------------------
## Under scale-change model case 2 : informative censoring
## -------------------------------------------------------

f1 <- f2 <- f3 <- f4 <- f5 <- matrix(NA, B, 4)
for (i in 1:B) {
    set.seed(i)
    dat <- simDat(200, a = c(0, 0), b = c(1, -1), type = "sc", indCen = FALSE)
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
## [1,]  0.8831710  0.9882635  -6.678685  1.662889 -0.03127153
## [2,] -0.8862846 -0.9741403 -19.448102 -1.670255  0.02254210
## [3,]  0.0000000  0.9996500   1.928973  7.548137  0.98554706
## [4,]  0.0000000 -0.9970690  -1.919424 -1.779842 -0.99167877
## > sapply(1:5, function(x) eval(parse(text = paste("matrix(apply(f", x, ", 2, median), 4)", sep = ""))))
##            [,1]       [,2]       [,3]      [,4]        [,5]
## [1,]  0.8817936  0.9755894  0.1947878  1.682990  0.03008002
## [2,] -0.8844487 -0.9655540 -0.9124529 -1.658278  0.04240854
## [3,]  0.0000000  0.9750789  1.9610650  1.675099  0.99290055
## [4,]  0.0000000 -0.9856637 -1.9264643 -1.716020 -0.99817651

## -------------------------------------------------------
## Under scale-change model, case 3: independent censoring
## -------------------------------------------------------

f1 <- f2 <- f3 <- f4 <- f5 <- matrix(NA, B, 4)
for (i in 1:B) {
    set.seed(i)
    dat <- simDat(200, a = c(1, -1), b = c(0, 0), type = "sc", indCen = TRUE)
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
##            [,1]       [,2]      [,3]      [,4]        [,5]
## [1,] -0.6273514 -0.5257213 -1.511254 -1.247291  1.03154219
## [2,]  0.5963205  0.4809668 17.324138  1.269810 -1.03648721
## [3,]  0.0000000 -0.5219780 -1.411929 -1.208746 -0.01009765
## [4,]  0.0000000  0.5151750  1.433242  1.360865 -0.01261398
## > sapply(1:5, function(x) eval(parse(text = paste("matrix(apply(f", x, ", 2, median), 4)", sep = ""))))
##            [,1]       [,2]       [,3]      [,4]         [,5]
## [1,] -0.6166167 -0.5152762 -0.2131766 -1.271713  0.964870296
## [2,]  0.5953391  0.4940052  0.2089009  1.326962 -1.010374144
## [3,]  0.0000000 -0.5353492 -1.4207072 -1.281072 -0.025677149
## [4,]  0.0000000  0.5091196  1.4120230  1.343471  0.004867741

## -------------------------------------------------------
## Under scale-change model case 3 : informative censoring
## -------------------------------------------------------

f1 <- f2 <- f3 <- f4 <- f5 <- matrix(NA, B, 4)
for (i in 1:B) {
    set.seed(i)
    dat <- simDat(200, a = c(1, -1), b = c(0, 0), type = "sc", indCen = FALSE)
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
##            [,1]       [,2]      [,3]      [,4]         [,5]
## [1,] -0.5922499 -0.4974843  2.055044 -1.103509  0.989257876
## [2,]  0.5641647  0.4710069  2.593204  1.158041 -1.019649499
## [3,]  0.0000000 -0.5193405 -1.494127 -1.133016 -0.018755820
## [4,]  0.0000000  0.4897615  1.419987  1.186428 -0.006947674
## > sapply(1:5, function(x) eval(parse(text = paste("matrix(apply(f", x, ", 2, median), 4)", sep = ""))))
##            [,1]       [,2]       [,3]      [,4]         [,5]
## [1,] -0.5915571 -0.4922137 -0.1857666 -1.178810  0.925451397
## [2,]  0.5606655  0.4785538  0.2642478  1.237901 -1.000159504
## [3,]  0.0000000 -0.5159088 -1.4670673 -1.223989 -0.024519623
## [4,]  0.0000000  0.4969365  1.4249257  1.236251 -0.002014639

## -------------------------------------------------------
## Under scale-change model, case 4: independent censoring
## -------------------------------------------------------

f1 <- f2 <- f3 <- f4 <- f5 <- matrix(NA, B, 4)
for (i in 1:B) {
    set.seed(i)
    dat <- simDat(200, a = c(1, 1), b = c(-1, -1), type = "sc", indCen = TRUE)
    f1[i,] <- coef(reReg(fm, data = dat))
    f2[i,] <- coef(reReg(fm, data = dat, method = "cox.HW"))
    f3[i,] <- coef(reReg(fm, data = dat, method = "am.GL"))
    f4[i,] <- tryCatch(coef(reReg(fm, data = dat, method = "am.XCHWY")),
                       error = function(e) rep(NA, 4))
    f5[i,] <- coef(reReg(fm, data = dat, method = "sc.XCYH"))
    if (i %% 20 == 0) print(i)
}

sapply(1:5, function(x) eval(parse(text = paste("matrix(colMeans(f", x, "), 4)", sep = ""))))
sapply(1:5, function(x) eval(parse(text = paste("matrix(apply(f", x, ", 2, median), 4)", sep = ""))))

## > sapply(1:5, function(x) eval(parse(text = paste("matrix(colMeans(f", x, "), 4)", sep = ""))))
##           [,1]      [,2]       [,3]        [,4]       [,5]
## [1,] -1.445607 -1.282392  54.095636   -2.056268  1.0115391
## [2,] -1.429322 -1.197263 629.029911   -2.472715  1.0446197
## [3,]  0.000000 -1.314208  -2.704081 1022.070326 -0.8843080
## [4,]  0.000000 -1.212871  -2.648869 -658.625437 -0.8646925
## > sapply(1:5, function(x) eval(parse(text = paste("matrix(apply(f", x, ", 2, median), 4)", sep = ""))))
##           [,1]      [,2]       [,3]      [,4]       [,5]
## [1,] -1.443116 -1.263079 -0.6560799 -2.366647  0.9857615
## [2,] -1.433454 -1.178175 -0.9972435 -2.603102  1.0137420
## [3,]  0.000000 -1.302996 -2.7414862 -2.444653 -1.0306269
## [4,]  0.000000 -1.192750 -2.6541133 -2.718752 -0.9753301

## -------------------------------------------------------
## Under scale-change model, case 4: informative censoring
## -------------------------------------------------------

f1 <- f2 <- f3 <- f4 <- f5 <- matrix(NA, B, 4)
for (i in 1:B) {
    set.seed(i)
    dat <- simDat(200, a = c(1, 1), b = c(-1, -1), type = "sc", indCen = FALSE)
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
##           [,1]      [,2]      [,3]      [,4]       [,5]
## [1,] -1.374198 -1.233134 12.949477 -1.729938  1.0530883
## [2,] -1.327007 -1.182547 52.996464 -2.189577  1.0982379
## [3,]  0.000000 -1.275683 -2.806613 -6.174569 -0.8793427
## [4,]  0.000000 -1.221069 -2.702266 -3.484050 -0.8328464
## > sapply(1:5, function(x) eval(parse(text = paste("matrix(apply(f", x, ", 2, median), 4)", sep = ""))))
##           [,1]      [,2]       [,3]      [,4]       [,5]
## [1,] -1.373814 -1.234607 -0.4802198 -1.988484  0.9796008
## [2,] -1.336483 -1.175824 -1.0299474 -2.360519  0.9997686
## [3,]  0.000000 -1.274581 -2.7381195 -2.158546 -0.9409923
## [4,]  0.000000 -1.218858 -2.6801064 -2.443495 -0.9009280

## ------------------------------------------------------------------------------------------
## checking variance estimation
## ------------------------------------------------------------------------------------------

set.seed(1)
dat <- simDat(200, a = c(1, -1), b = c(1, -1), type = "cox", indCen = TRUE)

summary(reReg(fm, data = dat, method = "cox.HW", se = "boot"))
## Call: reReg(formula = fm, data = dat, method = "cox.HW", se = "boot")
## Method: Huang-Wang Model 
## Coefficients (rate):
##    Estimate StdErr z.value   p.value    
## x1    0.943  0.160   5.899 < 2.2e-16 ***
## x2   -1.089  0.089 -12.228 < 2.2e-16 ***
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
## Coefficients (hazard):
##    Estimate StdErr z.value   p.value    
## x1    0.782  0.211   3.711 < 2.2e-16 ***
## x2   -1.136  0.134  -8.458 < 2.2e-16 ***
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(reReg(fm, data = dat, method = "cox.HW", se = "boot", control = list(parallel = TRUE)))
## Method: Huang-Wang Model 
## Coefficients (rate):
##    Estimate StdErr z.value   p.value    
## x1    0.943  0.165   5.713 < 2.2e-16 ***
## x2   -1.089  0.091 -11.940 < 2.2e-16 ***
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
## Coefficients (hazard):
##    Estimate StdErr z.value   p.value    
## x1    0.782  0.221   3.547 < 2.2e-16 ***
## x2   -1.136  0.137  -8.294 < 2.2e-16 ***
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(reReg(fm, data = dat, method = "cox.HW", se = "resam"))
## Call: reReg(formula = fm, data = dat, method = "cox.HW", se = "resam")
## Method: Huang-Wang Model 
## Coefficients (rate):
##    Estimate StdErr z.value   p.value    
## x1    0.943  0.139   6.805 < 2.2e-16 ***
## x2   -1.089  0.087 -12.525 < 2.2e-16 ***
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
## Coefficients (hazard):
##    Estimate StdErr z.value   p.value    
## x1    0.782  0.186   4.196 < 2.2e-16 ***
## x2   -1.136  0.104 -10.903 < 2.2e-16 ***
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


set.seed(1)
dat <- simDat(200, a = c(1, 1), b = c(-1, -1), type = "cox", indCen = FALSE)
fit1 <- reReg(fm, data = dat, se = "bootstrap")
fit2 <- reReg(fm, data = dat, method = "cox.HW", se = "res")
fit3 <- reReg(fm, data = dat, method = "am.GL", se = "boot")
fit4 <- reReg(fm, data = dat, method = "am.XCHWY", se = "res")
fit5 <- reReg(fm, data = dat, method = "sc.XCYH", se = "resa")

summary(fit1)
summary(fit2)
summary(fit3)
summary(fit4)
summary(fit5)

do <- function() {
    dat <- simDat(200, a = c(1, -1), b = c(1, -1), type = "cox", indCen = TRUE)
    fit <- reReg(fm, data = dat, method = "cox.HW", se = "resam")
    as.numeric(c(coef(fit), fit$alphaSE, fit$betaSE))
}

system.time(print(do()))

foo <- replicate(200, do())
