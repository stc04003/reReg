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
e
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


## simDat(1e4, c(1, -1), c(1, -1), indCen = TRUE, summary = TRUE, type = "am")

do <- function() {
    ## dat <- simDat(200, a = c(1, -1), b = c(1, -1), type = "am", indCen = TRUE)
    ## dat <- simDat(400, a = c(1, -1), b = c(1, -1), type = "am", indCen = FALSE)
    dat <- simDat(200, a = c(1, -1), b = c(1, -1), type = "cox", indCen = FALSE)
    fit1 <- reReg(fm, data = dat, method = "cox.HW", se = "res", B = 300)
    c(coef(fit1), fit1$alphaSE, fit1$betaSE)
}

do <- function() {
    dat <- simDat(200, a = c(1, -1), b = c(1, -1), type = "am", indCen = TRUE)
    ## dat <- simDat(400, a = c(1, -1), b = c(1, -1), type = "am", indCen = FALSE)
    ## dat <- simDat(200, a = c(1, -1), b = c(1, -1), type = "cox", indCen = FALSE)
    fit1 <- reReg(fm, data = dat, method = "am.X", se = "res", B = 200)
    c(coef(fit1), fit1$alphaSE, fit1$betaSE)
}

system.time(print(do()))

cl <- makePSOCKcluster(7)
setDefaultCluster(cl)
clusterExport(NULL, c("do", "fm"))
clusterEvalQ(NULL, library(reReg))
f1 <- parSapply(NULL, 1:100, function(z) tryCatch(do(), error = function(e) rep(NA, 8)))
stopCluster(cl)

apply(f1, 1, mean)[1:4]
apply(f1, 1, meanOut)[1:4]

rbind(apply(f1, 1, mean)[5:8],
      apply(f1, 1, meanOut)[5:8],
      apply(f1, 1, sd)[1:4],
      apply(f1, 1, sdOut)[1:4])

##  dat <- simDat(500, a = c(1, -1), b = c(1, -1), type = "am", indCen = TRUE)
##             x1         x2        x1        x2
## [1,] 0.2985686 0.07557614 0.5164932 0.2615754
## [2,] 0.2904914 0.07186844 0.5161991 0.2615754
## [3,] 0.1411801 0.06733162 0.5332168 0.2603365
## [4,] 0.1294701 0.06733162 0.5332168 0.2471383




## ------------------------------------------------------------------------------------------
## fitting cox.LWYY
## ------------------------------------------------------------------------------------------

set.seed(1)
fm <- reSurv(Time, id, event, status) ~ x1 + x2
dat <- simDat(n = 100, c(1, -1), c(1, -1), type = "cox", indCen = TRUE)

system.time(f1 <- reReg(fm, data = dat, se = NULL, method = "cox.GL"))
f1
summary(f1)
coef(f1)  # 1.142801 -1.081793  1.062973 -1.306367
plot(f1)

system.time(f2 <- reReg(fm, data = dat, se = "boot", method = "cox.GL", B = 20))
f2
summary(f2)
plot(f2)




e
## ------------------------------------------------------------------------------------------

do <- function(n = 100, a = c(1, -1), b = c(1, -1), type = "cox", indCen = TRUE) {
    dat <- simDat(n = n, a = a, b = b, type = type, indCen = indCen)
    f1 <- reReg(fm, data = dat, method = "cox.GL", se = "boot")
    c(coef(f1), sqrt(diag(vcov(f1)$alpha.vcov)), sqrt(diag(vcov(f1)$beta.vcov)))
}

set.seed(160)
do()

foo <- replicate(100, do(indCen = FALSE))
rowMeans(foo)


cl <- makePSOCKcluster(16)
setDefaultCluster(cl)
invisible(clusterExport(NULL, "do"))
invisible(clusterExport(NULL, "fm"))
invisible(clusterEvalQ(NULL, library(reReg)))
invisible(clusterEvalQ(NULL, library(survival)))

sim1 <- parSapply(NULL, 1:500, function(z) do())
sim2 <- parSapply(NULL, 1:500, function(z) do(indCen = FALSE))

stopCluster(cl)
