fm <- Recur(Time, id, event, status) ~ x1 + x2

## Joint Accelerated Mean Model
set.seed(1)
dat <- simSC(100, c(-1, 1), c(-1, 1), type = "cox")
(fit <- reReg(Recur(Time, id, event, status) ~ x1 + x2, 
              data = dat, method = "cox", B = 50))
summary(fit)

## Generalized Scale-Change Model
set.seed(1)
dat <- simSC(100, c(-1, 1), c(-1, 1), type = "sc")
(fit <- reReg(Recur(Time, id, event, status) ~ x1 + x2, 
              data = dat, method = "sc|.", B = 50))
summary(fit)
