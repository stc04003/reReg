library(reReg)

args(Recur)
?Recur
Recur(3:5)
Recur(6:1, id = rep(1:2, 3))

time1 <- c(1, 5, 7)
time2 <- c(3, 7, 9)
Recur(time1 %to% time2, id = c("A1", "A1", "A2"))
