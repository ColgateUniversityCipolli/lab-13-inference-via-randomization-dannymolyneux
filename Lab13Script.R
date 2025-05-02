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

#Question 3: randomization

#part a
R = 10000
n = 25
far.rand = tibble(xbars = rep(NA, R))
close.rand = tibble(xbars = rep(NA, R))
diff.rand = tibble(xbars = rep(NA, R))
# RANDOMIZE / SHUFFLE
for(i in 1:R){
  curr.far.rand <- far.vec *
    sample(x = c(-1, 1),
           size = n,
           replace = T)
  curr.close.rand <- close.vec *
    sample(x = c(-1, 1),
           size = n,
           replace = T)
  curr.diff.rand <- diff.vec *
    sample(x = c(-1, 1),
           size = n,
           replace = T)
  
  far.rand$xbars[i] <- mean(curr.far.rand)
  close.rand$xbars[i] <- mean(curr.close.rand)
  diff.rand$xbars[i] <- mean(curr.diff.rand)
}

#part b: p-value

mu0 = 0
delta.diff <- abs(mean(diff.vec) - mu0)
low.diff <- mu0 - delta.diff
high.diff <- mu0 + delta.diff
#p-values
(far.p = mean(far.rand$xbars <= mean(far.vec)))
(close.p = mean(close.rand$xbars >= mean(close.vec)))
(diff.p = mean(diff.rand$xbars <= low.diff) + mean(diff.rand$xbars >= high.diff))

#part c: confidence interval
mu0 = 0
R <- 1000
mu0.iterate <- 0.0001

far.mu.lower <- quantile(far.vec, 0.25)
close.mu.lower <- quantile(close.vec, 0.25)
diff.mu.lower <- quantile(diff.vec, 0.25)
repeat{
  far.rand <- tibble(xbars = rep(NA, R))
  close.rand <- tibble(xbars = rep(NA, R))
  diff.rand <- tibble(xbars = rep(NA, R))
  
  # PREPROCESSING: shift the data to be mean 0 under H0
  far.shift <- far.vec - far.mu.lower
  close.shift <- close.vec - close.mu.lower
  diff.shift <- diff.vec - diff.mu.lower
  # RANDOMIZE / SHUFFLE
  for(i in 1:R){
    far.curr.rand <- far.shift *
      sample(x = c(-1, 1),
             size = length(far.shift),
             replace = T)
    close.curr.rand <- close.shift *
      sample(x = c(-1, 1),
             size = length(close.shift),
             replace = T)
    diff.curr.rand <- diff.shift *
      sample(x = c(-1, 1),
             size = length(diff.shift),
             replace = T)
    
    far.rand$xbars[i] <- mean(far.curr.rand)
    close.rand$xbars[i] <- mean(close.curr.rand)
    diff.rand$xbars[i] <- mean(diff.curr.rand)
  }
  far.rand <- far.rand |>
    mutate(xbars = xbars + far.mu.lower) # shifting back
  close.rand <- close.rand |>
    mutate(xbars = xbars + close.mu.lower) # shifting back
  diff.rand <- diff.rand |>
    mutate(xbars = xbars + diff.mu.lower) # shifting back
  # p-value 
  far.delta <- abs(mean(far.vec) - far.mu.lower)
  close.delta <- abs(mean(close.vec) - close.mu.lower)
  diff.delta <- abs(mean(diff.vec) - diff.mu.lower)
  far.low <- far.mu.lower - far.delta
  close.low <- close.mu.lower - close.delta
  diff.low <- diff.mu.lower - diff.delta
  far.high<- far.mu.lower + far.delta
  close.high<- close.mu.lower + close.delta
  diff.high<- diff.mu.lower + diff.delta
  far.p.val <- mean(far.rand$xbars <= far.low) +
      mean(far.rand$xbars >= far.high)
  close.p.val <- mean(close.rand$xbars <= close.low) +
    mean(close.rand$xbars >= close.high)
  diff.p.val <- mean(diff.rand$xbars <= diff.low) +
    mean(diff.rand$xbars >= diff.high)
  if(diff.p.val < 0.05){
    diff.mu.lower <- diff.mu.lower + mu0.iterate
  }
  if(close.p.val < 0.05){
    close.mu.lower <- close.mu.lower + mu0.iterate
  }
  if(far.p.val < 0.05){
    far.mu.lower <- far.mu.lower + mu0.iterate
  }
  if((diff.p.val>0.05) & (close.p.val>0.05) & (far.p.val>0.05)){
    break
  }
}

far.mu.upper <- quantile(far.vec, 0.75)
close.mu.upper <- quantile(close.vec, 0.75)
diff.mu.upper <- quantile(diff.vec, 0.75)
repeat{
  far.rand <- tibble(xbars = rep(NA, R))
  close.rand <- tibble(xbars = rep(NA, R))
  diff.rand <- tibble(xbars = rep(NA, R))
  
  # PREPROCESSING: shift the data to be mean 0 under H0
  far.shift <- far.vec - far.mu.upper
  close.shift <- close.vec - close.mu.upper
  diff.shift <- diff.vec - diff.mu.upper
  # RANDOMIZE / SHUFFLE
  for(i in 1:R){
    far.curr.rand <- far.shift *
      sample(x = c(-1, 1),
             size = length(far.shift),
             replace = T)
    close.curr.rand <- close.shift *
      sample(x = c(-1, 1),
             size = length(close.shift),
             replace = T)
    diff.curr.rand <- diff.shift *
      sample(x = c(-1, 1),
             size = length(diff.shift),
             replace = T)
    
    far.rand$xbars[i] <- mean(far.curr.rand)
    close.rand$xbars[i] <- mean(close.curr.rand)
    diff.rand$xbars[i] <- mean(diff.curr.rand)
  }
  far.rand <- far.rand |>
    mutate(xbars = xbars + far.mu.upper) # shifting back
  close.rand <- close.rand |>
    mutate(xbars = xbars + close.mu.upper) # shifting back
  diff.rand <- diff.rand |>
    mutate(xbars = xbars + diff.mu.upper) # shifting back
  # p-value 
  far.delta <- abs(mean(far.vec) - far.mu.upper)
  close.delta <- abs(mean(close.vec) - close.mu.upper)
  diff.delta <- abs(mean(diff.vec) - diff.mu.upper)
  far.low <- far.mu.upper - far.delta
  close.low <- close.mu.upper - close.delta
  diff.low <- diff.mu.upper - diff.delta
  far.high<- far.mu.upper + far.delta
  close.high<- close.mu.upper + close.delta
  diff.high<- diff.mu.upper + diff.delta
  far.p.val <- mean(far.rand$xbars <= far.low) +
    mean(far.rand$xbars >= far.high)
  close.p.val <- mean(close.rand$xbars <= close.low) +
    mean(close.rand$xbars >= close.high)
  diff.p.val <- mean(diff.rand$xbars <= diff.low) +
    mean(diff.rand$xbars >= diff.high)
  if(diff.p.val < 0.05){
    diff.mu.upper <- diff.mu.upper - mu0.iterate
  }
  if(close.p.val < 0.05){
    close.mu.upper <- close.mu.upper - mu0.iterate
  }
  if(far.p.val < 0.05){
    far.mu.upper <- far.mu.upper - mu0.iterate
  }
  if((diff.p.val>0.05) & (close.p.val>0.05) & (far.p.val>0.05)){
    break
  }
}
(far.ci = c(far.mu.lower, far.mu.upper))
(close.ci = c(close.mu.lower, close.mu.upper))
(diff.ci = c(diff.mu.lower, diff.mu.upper))



