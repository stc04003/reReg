library(reReg)

args(Recur)
?Recur
Recur(3:5)
Recur(6:1, id = rep(1:2, 3))

time1 <- c(1, 5, 7)
time2 <- c(3, 7, 9)
ex3 <- Recur(time1 %to% time2, id = c("A1", "A1", "A2"))
head(ex3)

## Testing readmission data
data(readmission, package = "frailtypack")
head(readmission)
attach(readmission)

Recur(t.start %2% t.stop, id, event, death)
Recur(t.stop, id, event, death)

reSurv(t.start, t.stop, id, event, death)
reSurv(t.stop, id, event, death)

detach(readmission)


dat <- simSC(200, c(-1, 1), c(-1, 1))
head(dat)

with(dat, reSurv(Time, id, event, status))
with(dat, Recur(Time, id, event, status))
