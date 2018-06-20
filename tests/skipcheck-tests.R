#############################################################################################
## load packages and data readmission
#############################################################################################

library(parallel)
library(reReg)
library(reda)
library(gridExtra)
library(xtable)
data(readmission, package = "frailtypack")

R0 <- function(x) log(1 + x)
H0 <- function(x) log(1 + x) / 4

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

## t2 <- sapply(1:5, function(x)
##     eval(parse(text = paste("matrix(apply(f", x, ", 2, sd), 4)", sep = ""))))
## t1 <- sapply(1:5, function(x)
##     eval(parse(text = paste("matrix(apply(f", x, ", 2, meanOut), 4)", sep = ""))))
## t3 <- sapply(1:5, function(x)
##     eval(parse(text = paste("matrix(apply(f", x, ", 2, sdOut), 4)", sep = ""))))
## tab <- rbind(matrix(sapply(1:2, function(x) c(t1[,x], t2[,x], t3[,x])), 4),
##              matrix(sapply(3:4, function(x) c(t1[,x], t2[,x], t3[,x])), 4),
##              cbind(matrix(sapply(5, function(x) c(t1[,x], t2[,x], t3[,x])), 4), NA, NA, NA))
## xtable(tab, digits = 3)

do <- function(n = 200, a = c(1, -1), b = c(1, -1), type = "cox", indCen = TRUE) {
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

debug(doREFit.sc.XCYH)
debug(doREFit.sc.XCYH.resampling)
doREFit.sc.XCYH(DF = DF, engine = engine, stdErr = stdErr)
doREFit.sc.XCYH.resampling(DF = DF, engine = engine, stdErr = stdErr)

system.time(foo <- do())
foo

cl <- makePSOCKcluster(8)
setDefaultCluster(cl)
clusterExport(NULL, c("do", "fm"))
clusterEvalQ(NULL, library(reReg))
f1 <- parSapply(NULL, 1:500, function(z) tryCatch(do(), error = function(e) rep(NA, 40)))
f2 <- parSapply(NULL, 1:500, function(z) tryCatch(do(indCen = FALSE), error = function(e) rep(NA, 40)))
f3 <- parSapply(NULL, 1:500, function(z) tryCatch(do(type = "am"), error = function(e) rep(NA, 40)))
f4 <- parSapply(NULL, 1:500, function(z)
    tryCatch(do(type = "am", indCen = FALSE), error = function(e) rep(NA, 40)))
f5 <- parSapply(NULL, 1:100, function(z)
    tryCatch(do(a = c(1, 1), b = c(-1, -1), type = "sc", indCen = TRUE),
             error = function(e) rep(NA, 40)))
f6 <- parSapply(NULL, 1:100, function(z)
    tryCatch(do(a = c(1, 1), b = c(-1, -1), type = "sc", indCen = FALSE),
             error = function(e) rep(NA, 40)))
stopCluster(cl)

sim <- list(f1 = f1, f2 = f2, f3 = f3, f4 = f4, f5 = f5, f6 = f6)
save(sim, file = "sim.RData")

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
debug(do)
do(a = c(1, 1), b = c(-1, -1), type = "sc", indCen = TRUE)

set.seed(1)
dat <- simDat(200, a = c(1, -1), b = c(1, -1), type = "cox", indCen = TRUE)

system.time(fit1 <- reReg(fm, data = dat, method = "cox.HW", se = "boo"))
system.time(fit2 <- reReg(fm, data = dat, method = "cox.HW", se = "res"))

summary(fit1)

## Call: reReg(formula = fm, data = dat, method = "cox.HW", se = "boo")
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

summary(fit2)

## Call: reReg(formula = fm, data = dat, method = "cox.HW", se = "res")
## Method: Huang-Wang Model 
## Coefficients (rate):
##    Estimate StdErr z.value   p.value    
## x1    0.943  0.142   6.625 < 2.2e-16 ***
## x2   -1.089  0.095 -11.422 < 2.2e-16 ***
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
## Coefficients (hazard):
##    Estimate StdErr z.value   p.value    
## x1    0.782  0.187   4.177 < 2.2e-16 ***
## x2   -1.136  0.091 -12.518 < 2.2e-16 ***
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

tmp <- matrix(NA, 100, 8)
for (i in 1:100) {
    set.seed(i)
    dat <- simDat(200, a = c(1, -1), b = c(1, -1), type = "cox", indCen = TRUE)
    fit2 <- reReg(fm, data = dat, method = "cox.HW", se = "res")
    tmp[i,] <- c(coef(fit2), fit2$alphaSE, fit2$betaSE)
}

tmp <- matrix(NA, 100, 8)
for (i in 1:100) {
    set.seed(i)
    dat <- simDat(200, a = c(1, -1), b = c(1, -1), type = "am", indCen = TRUE)
    fit2 <- reReg(fm, data = dat, method = "am.XCHWY", se = "res")
    tmp[i,] <- c(coef(fit2), fit2$alphaSE, fit2$betaSE)
    if (i %% 10 == 0) print(i)
}

set.seed(1)
dat <- simDat(200, a = c(1, -1), b = c(1, -1), type = "am", indCen = TRUE)
fit1 <- reReg(fm, data = dat, method = "am.XCHWY", se = "boo")
fit2 <- reReg(fm, data = dat, method = "am.XCHWY", se = "res")

summary(fit1)
summary(fit2)


debug(doREFit.am.XCHWY.resampling)
doREFit.am.XCHWY.resampling(DF = DF, engine = engine, stdErr = stdErr)


coef(reReg(fm, data = dat))
coef(reReg(fm, data = dat, method = "cox.HW"))
coef(reReg(fm, data = dat, method = "am.GL"))
coef(reReg(fm, data = dat, method = "am.XCHWY"))
coef(reReg(fm, data = dat, method = "sc.XCYH"))

coef(reReg(fm, data = dat, method = "cox.HW", se = "res"))
coef(reReg(fm, data = dat, method = "sc.XCYH", se = "res"))

