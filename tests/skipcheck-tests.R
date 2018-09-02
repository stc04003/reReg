#############################################################################################
## load packages and data readmission
#############################################################################################

library(parallel)
library(reReg)
library(reda)
library(gridExtra)
library(xtable)
data(readmission, package = "frailtypack")

R0 <- function(x) log(1 + x) / .5
H0 <- function(x) log(1 + x) / 8

fm <- reSurv(Time, id, event, status) ~ x1 + x2    

meanOut <- function(dat, na.rm = TRUE) {
    dat[which(dat %in% boxplot(dat, plot = FALSE)$out)] <- NA
    dat <- dat[complete.cases(dat)]
    mean(dat, na.rm = na.rm)
}

sdOut <- function(dat, na.rm = TRUE) {
    dat[which(dat %in% boxplot(dat, plot = FALSE)$out)] <- NA
    dat <- dat[complete.cases(dat)]
    sd(dat, na.rm = na.rm)
}

varOut <- function(dat, na.rm = TRUE) {
    dat[which(dat %in% boxplot(dat, plot = FALSE)$out)] <- NA
    dat <- dat[complete.cases(dat)]
    var(dat, na.rm = na.rm)
}

varMatOut <- function(dat, na.rm = TRUE) {
    dat[which(dat %in% boxplot(dat, plot = FALSE)$out)] <- NA
    dat <- dat[complete.cases(dat),]
    var(dat, na.rm = na.rm)
}


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

reSurv(time1 = t.start, time2 = t.stop, id = id, event = event, status = death)
reSurv(time1 = t.stop, id = id, event = event, status = death)

detach(readmission)

## ------------------------------------------------------------------------------------------
## checking plotEvents
## ------------------------------------------------------------------------------------------
reObj <- with(readmission, reSurv(t.stop, id, event, death))

plot(reObj)
plot(reObj, xlab = "User X", ylab = "User Y", main = "User title")
plot(reObj, control = list(xlab = "User X", ylab = "User Y", main = "User title"))

plot(reObj, order = FALSE)
plot(reObj, order = FALSE, xlab = "User X", ylab = "User Y", main = "User title")
plot(reObj, order = FALSE, control = list(xlab = "User X", ylab = "User Y", main = "User title"))

plotEvents(reSurv(t.stop, id, event, death) ~ 1, data = readmission)
plotEvents(reSurv(t.stop, id, event, death) ~ 1, data = readmission,
           xlab = "User X", ylab = "User Y", main = "User title")
           
plotEvents(reSurv(t.stop, id, event, death) ~ sex, data = readmission)
plotEvents(reSurv(t.stop, id, event, death) ~ sex, data = readmission,
           xlab = "User X", ylab = "User Y", main = "User title")           
plotEvents(reSurv(t.stop, id, event, death) ~ sex + chemo, data = readmission)
plotEvents(reSurv(t.stop, id, event, death) ~ sex + chemo, data = readmission,
           xlab = "User X", ylab = "User Y", main = "User title",
           control = list(terminal.name = "User terminal",
                          recurrent.name = "User event"))

## multiple event types
reObj <- with(readmission, reSurv(t.stop, id, event * sample(1:3, 861, TRUE), death))
plot(reObj)
plot(reObj, xlab = "User X", ylab = "User Y", main = "User title")
plot(reObj, control = list(xlab = "User X", ylab = "User Y", main = "User title"))
plot(reObj, xlab = "User X", ylab = "User Y", main = "User title",
     control = list(recurrent.name = "User event"))
plot(reObj, xlab = "User X", ylab = "User Y", main = "User title",
     control = list(recurrent.type = letters[1:3]))

plotEvents(reSurv(t.stop, id, event, death) ~ 1, data = readmission)
plotEvents(reSurv(t.stop, id, event, death) ~ sex, data = readmission,
           control = list(xlab = "User X", ylab = "User Y", main = "User title"))
plotEvents(reSurv(t.stop, id, event, death) ~ sex + chemo, data = readmission,
           xlab = "User X", ylab = "User Y", main = "User title")

set.seed(123)
fm <- reSurv(t.stop, id, event * sample(1:3, 861, TRUE), death) ~ sex + chemo
plotEvents(fm, data = readmission)
plotEvents(fm, data = readmission, xlab = "User X", ylab = "User Y", main = "User title",
           control = list(recurrent.type = letters[1:3]))

## ------------------------------------------------------------------------------------------
## checking plotCSM
## ------------------------------------------------------------------------------------------
reObj <- with(readmission, reSurv(t.stop, id, event, death))

plot(reObj, CSM = TRUE)
plot(reObj, CSM = TRUE, xlab = "User X", ylab = "User Y", main = "User title")
plot(reObj, CSM = TRUE, control = list(xlab = "User X", ylab = "User Y", main = "User title"))

plotCSM(reSurv(t.stop, id, event, death) ~ 1, data = readmission)
plotCSM(reSurv(t.stop, id, event, death) ~ sex, data = readmission)
plotCSM(reSurv(t.stop, id, event, death) ~ sex, data = readmission, onePanel = TRUE)

set.seed(123)
fm <- reSurv(t.stop, id, event * sample(1:3, 861, TRUE), death) ~ sex + chemo
plotCSM(fm, data = readmission)
plotCSM(fm, data = readmission, xlab = "User X", ylab = "User Y", main = "User title",
           control = list(recurrent.type = letters[1:3]))
plotCSM(fm, data = readmission, xlab = "User X", ylab = "User Y", main = "User title",
           control = list(recurrent.name = "Types", recurrent.type = letters[1:3]))

plotCSM(fm, data = readmission, recurrent.name = "Types", recurrent.type = letters[1:3])

plotCSM(fm, data = readmission, xlab = "User X", ylab = "User Y", main = "User title",
           control = list(recurrent.type = letters[1:3])) + ggplot2::theme(legend.position = "bottom")

## ------------------------------------------------------------------------------------------
## checking simulated data generator
## ------------------------------------------------------------------------------------------

simDat(1e4, c(1, -1), c(1, -1), indCen = TRUE, summary = TRUE)
simDat(1e4, c(1, -1), c(1, -1), indCen = FALSE, summary = TRUE)
simDat(1e4, c(1, -1), c(1, -1), indCen = TRUE, summary = TRUE, type = "am")
simDat(1e4, c(1, -1), c(1, -1), indCen = FALSE, summary = TRUE, type = "am")

simDat(1e4, c(1, -1), c(1, -1), indCen = TRUE, summary = TRUE, type = "sc")
simDat(1e4, c(1, -1), c(1, -1), indCen = FALSE, summary = TRUE, type = "sc")
simDat(1e4, c(0, 0), c(1, -1), indCen = TRUE, summary = TRUE, type = "sc")
simDat(1e4, c(0, 0), c(1, -1), indCen = FALSE, summary = TRUE, type = "sc")
simDat(1e4, c(1, -1), c(0, 0), indCen = TRUE, summary = TRUE, type = "sc")
simDat(1e4, c(1, -1), c(0, 0), indCen = FALSE, summary = TRUE, type = "sc")
simDat(1e4, c(1, 1), c(-1, -1), indCen = TRUE, summary = TRUE, type = "sc")
simDat(1e4, c(1, 1), c(-1, -1), indCen = FALSE, summary = TRUE, type = "sc")


## ------------------------------------------------------------------------------------------
## checking point esitmation
## ------------------------------------------------------------------------------------------

fm <- reSurv(Time, id, event, status) ~ x1 + x2
B <- 200

do <- function(n = 100, a = c(1, -1), b = c(1, -1), type = "cox", indCen = TRUE) {
    dat <- simDat(n = n, a = a, b = b, type = type, indCen = indCen)
    f1 <- reReg(fm, data = dat, se = "boot", B = 300)
    f2 <- reReg(fm, data = dat, method = "cox.HW", se = "resam", B = 300)
    invisible(capture.output(f3 <- reReg(fm, data = dat, method = "am.GL", se = "boot", B = 300)))
    f4 <- reReg(fm, data = dat, method = "am.XCHWY", se = "res", B = 300)
    f5 <- reReg(fm, data = dat, method = "sc.XCYH", se = "resam", B = 300)
    c(coef(f1), f1$alphaSE, f1$betaSE,
      coef(f2), f2$alphaSE, f2$betaSE,
      coef(f3), f3$alphaSE, f3$betaSE,
      coef(f4), f4$alphaSE, f4$betaSE,
      coef(f5), f5$alphaSE, f5$betaSE)
}

cl <- makePSOCKcluster(8)
setDefaultCluster(cl)
clusterExport(NULL, c("do", "fm"))
clusterEvalQ(NULL, library(reReg))
f1 <- parSapply(NULL, 1:500, function(z) tryCatch(do(), error = function(e) rep(NA, 40)))
f2 <- parSapply(NULL, 1:500, function(z) tryCatch(do(indCen = FALSE), error = function(e) rep(NA, 40)))
f3 <- parSapply(NULL, 1:500, function(z) tryCatch(do(type = "am"), error = function(e) rep(NA, 40)))
f4 <- parSapply(NULL, 1:500, function(z)
    tryCatch(do(type = "am", indCen = FALSE), error = function(e) rep(NA, 40)))
f5 <- parSapply(NULL, 1:500, function(z)
    tryCatch(do(a = c(1, 1), b = c(-1, -1), type = "sc", indCen = TRUE),
             error = function(e) rep(NA, 40)))
f6 <- parSapply(NULL, 1:500, function(z)
    tryCatch(do(a = c(1, 1), b = c(-1, -1), type = "sc", indCen = FALSE),
             error = function(e) rep(NA, 40)))
stopCluster(cl)

## sim <- list(f1 = f1, f2 = f2, f3 = f3, f4 = f4, f5 = f5, f6 = f6)
## save(sim, file = "sim.RData")

datPre <- function(f0) {
    PE <- apply(f0, 1, mean, na.rm = T)[c(1:4, 9:12, 17:20, 25:28, 33:36)]
    ESE <- apply(f0, 1, sd, na.rm = T)[c(1:4, 9:12, 17:20, 25:28, 33:36)]
    ASE <- apply(f0, 1, mean, na.rm = T)[c(1:4, 9:12, 17:20, 25:28, 33:36) + 4]
    matrix(rbind(matrix(PE, 4), matrix(ESE, 4), matrix(ASE, 4)), 4)
}

datPre2 <- function(f0) {
    PE <- apply(f0, 1, meanOut)[c(1:4, 9:12, 17:20, 25:28, 33:36)]
    ESE <- apply(f0, 1, sdOut)[c(1:4, 9:12, 17:20, 25:28, 33:36)]
    ASE <- apply(f0, 1, meanOut)[c(1:4, 9:12, 17:20, 25:28, 33:36) + 4]
    matrix(rbind(matrix(PE, 4), matrix(ESE, 4), matrix(ASE, 4)), 4)
}

set.seed(1)
do(a = c(1, 1), b = c(-1, -1), type = "sc", indCen = TRUE)

set.seed(123)
dat <- simDat(200, a = c(1, 1), b = c(-1, -1), type = "sc", indCen = TRUE)
system.time(fit2 <- reReg(fm, data = dat, method = "am.XC", se = "res"))
summary(fit2)
system.time(fit22 <- reReg(fm, data = dat, method = "am.XC", se = "boo"))
summary(fit22)

system.time(fit3 <- reReg(fm, data = dat, method = "sc.X", se = "res"))
summary(fit3)





## ------------------------------------------------------------------------------------------
## checking am.GL
## ------------------------------------------------------------------------------------------

do <- function(n = 100, a = c(1, -1), b = c(1, -1), type = "cox", indCen = TRUE) {
    fm <- reSurv(Time, id, event, status) ~ x1 + x2
    dat <- simDat(n = n, a = a, b = b, type = type, indCen = indCen)
    invisible(capture.output(f3 <- reReg(fm, data = dat, method = "am.GL", se = "boot", B = 200)))
    f4 <- reReg(fm, data = dat, method = "am.XCHWY", se = "re", B = 200)
    ## c(coef(f3), f3$alphaSE, f3$betaSE,
    c(coef(f3),
      sqrt(diag$varMatOut(f3$SEmat))[1:2],
      sqrt(diag$varMatOut(f3$SEmat))[3:4],
      coef(f4), f4$alphaSE, f4$betaSE)
}

do <- function(n = 100, a = c(1, -1), b = c(1, -1), type = "cox", indCen = TRUE) {
    fm <- reSurv(Time, id, event, status) ~ x1 + x2
    dat <- simDat(n = n, a = a, b = b, type = type, indCen = indCen)
    invisible(capture.output(f3 <- reReg(fm, data = dat, method = "am.GL", se = "boot", B = 200)))
    c(coef(f3), f3$alphaSE, f3$betaSE)
}

do <- function(n = 100, a = c(1, -1), b = c(1, -1), type = "cox", indCen = TRUE) {
    fm <- reSurv(Time, id, event, status) ~ x1 + x2
    dat <- simDat(n = n, a = a, b = b, type = type, indCen = indCen)
    f3 <- reReg(fm, data = dat, method = "am.GL", se = NULL)
    f4 <- reReg(fm, data = dat, method = "am.XCHWY", se = NULL)
    c(coef(f3), coef(f4))
}

cl <- makePSOCKcluster(8)
## cl <- makePSOCKcluster(16)
setDefaultCluster(cl)
invisible(clusterExport(NULL, "do"))
invisible(clusterEvalQ(NULL, library(reReg)))
fo <- parSapply(NULL, 1:10, function(z) {
    set.seed(z); do(n = 200, type = "am", indCen = FALSE)})
    ## set.seed(z); do(type = "am", indCen = TRUE)})
stopCluster(cl)

cbind(apply(fo, 1, mean, na.rm = TRUE),
      apply(fo, 1, meanOut),
      apply(fo, 1, sdOut))


set.seed(6)
dat <- simDat(n = 100, a = c(1, -1), b = c(1, -1), type = "am", indCen = FALSE)
f3 <- reReg(fm, data = dat, method = "am.GL", se = "boot", B = 200)
summary(f3)

