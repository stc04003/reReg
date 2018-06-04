library(reReg)

## plotEvents examples
data(readmission)

reObj <- with(readmission, reSurv(t.stop, event, death, id))
plot(reObj)

plotEvents(reSurv(t.stop, event, death, id) ~ 1, data = readmission)
plotEvents(reSurv(t.stop, event, death, id) ~ sex, data = readmission)
plotEvents(reSurv(t.stop, event, death, id) ~ sex + chemo, data = readmission)
plotEvents(reSurv(t.stop, event * sample(1:3, 861, TRUE), death, id) ~
               sex + chemo, data = readmission)

## plotCMF examples

plot(reObj, CMF = TRUE)

plotCMF(reObj)
plotCMF(reObj ~ 1, data = readmission)
plotCMF(reObj ~ sex, data = readmission)
plotCMF(reObj ~ sex + chemo, data = readmission)


plotCMF(reSurv(t.stop, event, death, id) ~ 1, data = readmission)
plotCMF(reSurv(t.stop, event, death, id) ~ sex, data = readmission)
plotCMF(reSurv(t.stop, event, death, id) ~ sex + chemo, data = readmission)
plotCMF(reSurv(t.stop, event * sample(1:3, 861, TRUE), death, id) ~
            sex + chemo, data = readmission)
plotCMF(reSurv(t.stop, event * sample(1:3, 861, TRUE), death, id) ~
            sex + chemo, data = readmission,
        control = list(recurrent.type = letters[1:3]))

plotCMF(reSurv(t.stop, event, death, id) ~ sex + chemo, data = readmission,
        control = list(title = "Some title"))

## Same plot with collapse = TRUE
plotCMF(reSurv(t.stop, event, death, id) ~ sex, data = readmission, onePanel = TRUE)
plotCMF(reSurv(t.stop, event * sample(1:3, 861, TRUE), death, id) ~ sex,
        data = readmission, onePanel = TRUE)
plotCMF(reSurv(t.stop, event, death, id) ~ sex + chemo, data = readmission, onePanel = TRUE)

## reReg exapmles

## readmission data
data(readmission)
set.seed(123)
## Acceralted Mean Model
(fit <- reReg(reSurv(t.stop, event, death, id) ~ sex + chemo,
              data = subset(readmission, id < 50),
              method = "am.XCHWY", se = "resampling", B = 20))
summary(fit)

## Generalized Scale-Change Model
set.seed(123)
(fit <- reReg(reSurv(t.stop, event, death, id) ~ sex + chemo,
              data = subset(readmission, id < 50),
              method = "sc.XCYH", se = "resampling", B = 20))
summary(fit)

## Not run:

## simulation data
simDat <- function(n, a, b, latent = FALSE) {
    ## setting rate function
    Lam.f <- function(t, z, x, w) .5 * z * exp(-x + w) * log(1 + t * exp(x))
    Lam.f0 <- function(t) .5 * log(1 + t)
    invLam.f  <- function(t, z, x, w) (exp((2 / z) * exp(x - w) * t )- 1) / exp(x)
    ## setting hazard funciton
    ## Haz.f0 <- function(t) .5 * log(1 + t) # assume constant hazard for now
    dat <- NULL
    for (id in 1:n) {
        z <- ifelse(latent, rgamma(1, 4, 4), 1)
        x1 <- rnorm(1)
        x2 <- rnorm(1)
        x <- c(x1, x2)
        cen <- rexp(1, z * exp(x \%*\% b) / 60) ## this gives constant hazard of 1/60
        y <- min(cen, 60)
        D <- 1 * (cen == y)
        tmpt <- NULL
        while(sum(tmpt) < Lam.f(y, z, c(x \%*\% a), c(x \%*\% b))) {
            tmpt <- c(tmpt, rexp(1))
        }
        m <- length(tmpt) - 1
        if (m > 0) {
            tt <- invLam.f(cumsum(tmpt[1:m]), z, c(x \%*\% a), c(x \%*\% b))
            dat <- rbind(dat, cbind(ID = id, Time = c(tt[order(tt)], y),
                                    event = c(rep(1, m), 0), status = c(rep(0, m), D),
                                    Z = z, M = m, X1 = x1, X2 = x2))
        } else {
            dat <- rbind(dat, cbind(ID = id, Time = y, event = 0, status = D,
                                    Z = z, M = m, X1 = x1, X2 = x2))
        }
    }
    return(data.frame(dat))
}
set.seed(2017)
dat <- simDat(200, c(0, 0), c(0, 0), TRUE) ## generate data under informative censoring
fm <- reSurv(Time, event, status, ID) ~ X1 + X2
fit.HW <- reReg(fm, data = dat, method = "cox.HW", B = 5)
## End(Not run)
