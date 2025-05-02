#Lab 13
#question 1
#part a
library(tidyverse)
library(e1071)
library(ggplot2)
data = read_csv("zebrafinches.csv")
far.vec = data$further
mu0 = 0
n = 25
far.t.stat = (mean(far.vec) - mu0)/(sd(far.vec)/sqrt(n)) #far t statistic
far.pdf = dnorm(far.t.stat)
far.error = (skewness(far.vec)/sqrt(n))*((2*far.t.stat^2 + 1)/6)*far.pdf

#part b
t.stat = seq(-10,10, by = 0.01)
error.vec = (skewness(far.vec)/sqrt(n))*((2*t.stat^2 + 1)/6)*dnorm(t.stat)
error.dat = data.frame(t = t.stat, error = error.vec)
error.plot = ggplot(data = error.dat, aes(x = t, y = error))+
  geom_line(color = "red") +
  labs(title = "Error with respect to t statistic for further data")
error.plot

#part c
alpha = 0.05
t = qnorm(alpha) #t statistic at alpha = 0.05
#required sample size to get a tail probability within 10% of alpha
n.required = ((skewness(far.vec)/(6*.1*alpha))*((2*t^2 + 1))*dnorm(t))^2

#question 2

#part a
far.vec = data$further
close.vec = data$closer
diff.vec = data$diff
R = 10000 #simulations
n = 25
far.resamples = tibble(t.stats = rep(NA, R))
close.resamples = tibble(t.stats = rep(NA, R))
diff.resamples = tibble(t.stats = rep(NA, R))
#resamples for all 3 types of data
for(i in 1:R){
  curr.far.resample = sample(far.vec,
                          size = n,
                          replace = T)
  curr.close.resample = sample(close.vec,
                             size = n,
                             replace = T)
  curr.diff.resample = sample(diff.vec,
                             size = n,
                             replace = T)
  #resampled t statistics for close, far, and diff data
  far.resamples$t.stats[i] = mean(curr.far.resample)/(sd(far.vec)/sqrt(n))
  close.resamples$t.stats[i] = mean(curr.close.resample)/(sd(close.vec)/sqrt(n))
  diff.resamples$t.stats[i] = mean(curr.diff.resample)/(sd(diff.vec)/sqrt(n))
}
#shifted resamples
resamples.null.far = far.resamples$t.stats - mean(far.resamples$t.stats)
resamples.null.close = close.resamples$t.stats - mean(close.resamples$t.stats)
resamples.null.diff = diff.resamples$t.stats - mean(diff.resamples$t.stats)
#mean of shifted resamples
mean(resamples.null.far)
mean(resamples.null.close)
mean(resamples.null.diff)

#part b
far.t.stat = (mean(far.vec) - mu0)/(sd(far.vec)/sqrt(n)) #far t statistic
close.t.stat = (mean(close.vec) - mu0)/(sd(close.vec)/sqrt(n)) #far t statistic
diff.t.stat = (mean(diff.vec) - mu0)/(sd(diff.vec)/sqrt(n)) #far t statistic
#p values for each test
(far.p = mean(resamples.null.far <= far.t.stat))
(close.p = mean(resamples.null.close >= close.t.stat))
(diff.p = mean(resamples.null.diff >= abs(diff.t.stat)) + mean(resamples.null.diff <= -abs(diff.t.stat)))

#part c
#5th percentiles of shifted resamples
(far.5th.percentile = quantile(resamples.null.far, 0.05))
(close.5th.percentile = quantile(resamples.null.close, 0.05))
(diff.5th.percentile = quantile(resamples.null.diff, 0.05))

#part d
library(boot)
boot.mean <- function(d, i){
  mean(d[i])
}
#bootstrap confidence intervals
boots.far <- boot(data = far.vec,
              statistic = boot.mean,
              R = R)
boot.far.ci = boot.ci(boots.far, type="bca")

boot.far.ci$bca[4:5]

boots.close <- boot(data = close.vec,
                  statistic = boot.mean,
                  R = R)
boot.close.ci = boot.ci(boots.close, type="bca")

boot.close.ci$bca[4:5]

boots.diff <- boot(data = diff.vec,
                  statistic = boot.mean,
                  R = R)
boot.diff.ci = boot.ci(boots.diff, type="bca")

boot.diff.ci$bca[4:5]
