library(reReg)

## -----------------------------------------------------------------------------------
## With simulated data
## -----------------------------------------------------------------------------------

args(simDat)
simDat(200, c(1, 1), c(1, 1))
simDat(200, c(1, 1), c(1, 1), type = "am")
simDat(200, c(1, 1), c(1, 1), type = "sc")

simDat(200, c(1, 1), c(1, 1), indCen = FALSE)
simDat(200, c(1, 1), c(1, 1), type = "am", indCen = FALSE)
simDat(200, c(1, 1), c(1, 1), type = "sc", indCen = FALSE)

dat <- simDat(200, c(1, 1), c(1, 1))
fm <- reSurv(Time, id, event, status) ~ x1 + x2
reReg(fm, data = dat)

e

## -----------------------------------------------------------------------------------
## With readmission data
## -----------------------------------------------------------------------------------
data(readmission, package = "frailtypack")
fm <- reSurv(t.stop, id, event, death) ~ sex + chemo
