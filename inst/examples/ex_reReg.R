data(simDat)

## Nonparametric estimate
plot(reReg(Recur(t.start %to% t.stop, id, event, status) ~ 1, data = simDat, B = 50))

fm <- Recur(t.start %to% t.stop, id, event, status) ~ x1 + x2
## Fit the Cox rate model
summary(reReg(fm, data = simDat, model = "cox", B = 50))
## Fit the joint Cox/Cox model
summary(reReg(fm, data = simDat, model = "cox|cox", B = 50))
## Fit the scale-change rate model
summary(reReg(fm, data = simDat, model = "sc", B = 50, se = "sand"))
