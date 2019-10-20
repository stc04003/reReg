set.seed(1)
dat <- simSC(80, c(-1, 1), c(-1, 1), type = "am")
## Accelerated Mean Model
(fit <- reReg(reSurv(Time, id, event, status) ~ x1 + x2, 
              data = dat, method = "am.XCHWY", se = "resampling", B = 20))
summary(fit)

## Generalized Scale-Change Model
set.seed(1)
dat <- simSC(100, c(-1, 1), c(-1, 1), type = "sc")
(fit <- reReg(reSurv(Time, id, event, status) ~ x1 + x2, 
              data = dat, method = "sc.XCYH", se = "resampling", B = 20))
summary(fit)