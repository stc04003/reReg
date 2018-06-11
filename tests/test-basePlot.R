library(reReg)

fm <- reSurv(Time, id, event, status) ~ x1 + x2

## Use simulated data
R0 <- function(x) log(1 + x)
H0 <- function(x) log(1 + x) / 4

## Cox type

dat <- simDat(200, c(1, 1), c(1, 1), type = "cox", indCen = TRUE)
fit <- reReg(fm, data = dat)

summary(fit)
fit

dat <- as.tibble(fit$DF) %>%
    mutate(R0 = fit$rate0(Time),
           R0.upper = fit$rate0.upper(Time),
           R0.lower = fit$rate0.lower(Time),
           H0 = fit$haz0(Time),
           H0.upper = fit$haz0.upper(Time),
           H0.lower = fit$haz0.lower(Time))
dat

ggplot(data = dat, aes(x = Time, y = R0)) + geom_line() +
    geom_line(aes(x = Time, y = R0.lower)) +
    geom_line(aes(x = Time, y = R0.upper))


ggplot(data = dat, aes(x = Time, y = H0)) + geom_line() +
    geom_line(aes(x = Time, y = H0.lower)) +
    geom_line(aes(x = Time, y = H0.upper))

ggplot(data = dat, aes(x = Time, y = H0)) + geom_line() 
ggplot(data = dat, aes(x = Time, y = H0)) + geom_smooth(se = FALSE, color = 1) 

str(fit)

library(tidyverse)

ggplot(aes(x = fit$t0.rate, y = fit$rate0(fit$t0.rate))) + geom_line()
