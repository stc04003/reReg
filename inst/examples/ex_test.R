library(reReg)

args(Recur)
?Recur
Recur(3:5)
Recur(6:1, id = rep(1:2, 3))
is.Recur(Recur(3:5))
    
time1 <- c(1, 5, 7)
time2 <- c(3, 7, 9)
ex3 <- Recur(time1 %to% time2, id = c("A1", "A1", "A2"))
ex3
head(ex3)

## #########################################################
## Testing readmission data
## Some subjects in this data have multiple terminal events 
## #########################################################
data(readmission, package = "frailtypack")
head(readmission)
attach(readmission)

Recur(t.start %2% t.stop, id, event, death)
Recur(t.stop, id, event, death)

reSurv(t.start, t.stop, id, event, death)
reSurv(t.stop, id, event, death)

detach(readmission)

testDF <- function(dat) {
    ## `dat` comes from simSC
    ## `fit1` is from reSurv
    ## `fit2` is from Recur 
    fit1 <- with(dat, reSurv(Time, id, event, status))
    fit2 <- with(dat, Recur(Time, id, event, status))
    reDF2 <- as.data.frame(fit2@.Data)
    all(identical(fit1$reDF$Time, reDF2$time2),
        ## all.equal(fit1$reDF$id, reDF2$id),
        all.equal(fit1$reDF$id, as.numeric(levels(fit2@ID))[fit2@ID]),
        identical(fit1$reDF$recType, reDF2$event),
        identical(fit1$reDF$status, reDF2$terminal))
}

testDF(simSC(200, c(-1, 1), c(-1, 1)))
table(replicate(1000, testDF(simSC(200, c(-1, 1), c(-1, 1)))))

## ID in different orders
test1 <- function() {
    dat <- simSC(200, c(-1, 1), c(-1, 1))
    dat <- subset(dat, id %in% sample(1:200, 50))
    testDF(dat)
}

table(replicate(1000, test1()))

test2 <- function() {
    dat <- simSC(200, c(-1, 1), c(-1, 1))
    dat$id <- rep(sample(1:200, 200), table(dat$id))
    testDF(dat)
}

table(replicate(1000, test2()))

## Different event types

dat <- simSC(200, c(-1, 1), c(-1, 1))
dat$event <- dat$event * sample(1:3, nrow(dat), TRUE)
testDF(dat)

head(dat)

debug(testDF)


## #########################################################
## Testing plots
## #########################################################
set.seed(1)
n <- 20
dat <- simSC(n, c(-1, 1), c(-1, 1))
dat$x3 <- rep(sample(0:1, n, TRUE), table(dat$id))
dat

fit1 <- with(dat, reSurv(Time, id, event, status))
fit2 <- with(dat, Recur(Time, id, event, status))
str(fit2)

plotCSM(with(dat, Recur(Time, id, event, status)))
plotCSM(Recur(Time, id, event, status) ~ 1, data = dat)
plotCSM(Recur(Time, id, event, status) ~ x1, data = dat)
plotCSM(Recur(Time, id, event, status) ~ x1 + x3, data = dat)

plotCSM(with(dat, Recur(Time, id, event, status)), adjrisk = FALSE)
plotCSM(Recur(Time, id, event, status) ~ 1, data = dat, adjrisk = FALSE)
plotCSM(Recur(Time, id, event, status) ~ x1, data = dat, adjrisk = FALSE)
plotCSM(Recur(Time, id, event, status) ~ x1 + x3, data = dat, adjrisk = FALSE)

plotCSM(with(dat, Recur(Time, id, event, status)), smooth = TRUE)
plotCSM(Recur(Time, id, event, status) ~ 1, data = dat, smooth = TRUE)
plotCSM(Recur(Time, id, event, status) ~ x1, data = dat, smooth = TRUE)
plotCSM(Recur(Time, id, event, status) ~ x1 + x3, data = dat, smooth = TRUE)

plotCSM(Recur(Time, id, event, status) ~ x1, data = dat, onePanel = TRUE)
plotCSM(Recur(Time, id, event, status) ~ x1 + x3, data = dat, onePanel = TRUE)

plotEvents(with(dat, Recur(Time, id, event, status)))
plotEvents(Recur(Time, id, event, status) ~ 1, data = dat)
plotEvents(Recur(Time, id, event, status) ~ x1, data = dat)
plotEvents(Recur(Time, id, event, status) ~ x1 + x3, data = dat)

