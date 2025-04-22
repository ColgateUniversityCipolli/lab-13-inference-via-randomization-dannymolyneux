#Lab 13
#question 1
#part a
library(tidyverse)
library(e1071)
data = read_csv("zebrafinches.csv")
far.vec = data$further
mu0 = 0
n = 25
far.t.stat = (mean(far.vec) - mu0)/(sd(far.vec)/sqrt(n)) #far t statistic
far.pdf = dnorm(far.t.stat)
far.error = (skewness(far.vec)/sqrt(n))*((2*far.t.stat^2 + 1)/6)*far.pdf
