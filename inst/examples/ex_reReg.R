fm <- Recur(Time, id, event, status) ~ x1 + x2

## Joint Cox/Accelerated Mean Model
set.seed(1)
dat <- simSC(50, c(-1, 1), c(-1, 1), type = "cox|am")
(fit <- reReg(Recur(Time, id, event, status) ~ x1 + x2, 
              data = dat, method = "cox|am", B = 50))
summary(fit)

## Generalized Scale-Change Model
set.seed(1)
dat <- simSC(50, c(-1, 1), c(-1, 1), type = "sc")
(fit <- reReg(Recur(Time, id, event, status) ~ x1 + x2, 
              data = dat, method = "sc|.", B = 50))
summary(fit)
