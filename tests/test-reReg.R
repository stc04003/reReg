library(reReg)

## -----------------------------------------------------------------------------------
## With simulated data
## -----------------------------------------------------------------------------------

args(simDat)
fm <- reSurv(Time, id, event, status) ~ x1 + x2

simDat(200, c(1, 1), c(1, 1))
simDat(200, c(1, 1), c(1, 1), type = "am")
simDat(200, c(1, 1), c(1, 1), type = "sc")

simDat(200, c(1, 1), c(1, 1), indCen = FALSE)
simDat(200, c(1, 1), c(1, 1), type = "am", indCen = FALSE)
simDat(200, c(1, 1), c(1, 1), type = "sc", indCen = FALSE)


dat <- simDat(200, c(1, 1), c(1, 1))
coef(reReg(fm, data = dat))


do <- function(n = 200, a = c(1, 1), b = c(1, 1), type = "cox", indCen = TRUE) {
    dat <- simDat(n, a, b, indCen, type)
    c(coef(reReg(fm, data = dat)),
      coef(reReg(fm, data = dat, method = "cox.HW")),
      ## coef(reReg(fm, data = dat, method = "am.GL")),
      coef(reReg(fm, data = dat, method = "am.XCHWY")),
      coef(reReg(fm, data = dat, method = "sc.XCYH")))    
}

set.seed(1527) 
do(a = c(-1, -1), b = c(1, 1), indCen = TRUE, type = "sc")
set.seed(1526)  ## and 1527
dat <- simDat(200, a = c(-1, -1), b = c(1, 1), indCen = TRUE, type = "sc")
coef(reReg(fm, data = dat)) ## ok
coef(reReg(fm, data = dat, method = "cox.HW")) ## ok
## coef(reReg(fm, data = dat, method = "am.GL"))
coef(reReg(fm, data = dat, method = "am.XCHWY")) ## errors
coef(reReg(fm, data = dat, method = "sc.XCYH")) ## ok




do(type = "am")

foo <- matrix(NA, 500, 16)
for (i in 1:500) {
    set.seed(i + 1500)
    ## foo[i,] <- do(indCen = TRUE, type = "am")
    foo[i,] <- do(a = c(-1, -1), b = c(1, 1), indCen = TRUE, type = "sc")
    print(i)
}

matrix(colMeans(foo), 4)

## Cox, both censoring: coxph, cox.HW, sc.xcYH ok!
## am, both censoring: am.XCHWY and sc.XCYH ok!
e

## -----------------------------------------------------------------------------------
## With readmission data
## -----------------------------------------------------------------------------------
data(readmission, package = "frailtypack")
fm <- reSurv(t.stop, id, event, death) ~ sex + chemo
