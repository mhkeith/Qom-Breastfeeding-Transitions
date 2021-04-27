#### Qom breastfeeding transitions ####
#### Time-to-event analysis: base hazard function screening ####

## load packages
library(tidyverse)
library(loo)
library(rstanarm)
library(cowplot)

## load data
dat.CF <- read.csv(file="~/QomTTCF.dat.csv") #time-to-complementary feeding, repeated observations 
dat.wean <- read.csv(file="~/QomTTWean.dat.csv") #time-to-weaning, single observations

##------------------- time-to-complementary feeding models -------------------##

## code events with Surv object notation
## 0=right censored, 1=event observed, 2=left censored, 3=interval censored
dat.CF$event[dat.CF$leftcensored==1] <- 2
## event times must be >0 since log(t) is used in parametric equations
dat.CF$agewks.st[dat.CF$agewks.st==0] <- 0.5

## fit exponential, Weibull, Gompertz, m-spline, and Weibull-aft base hazard functions
ttcf.exp <- stan_surv(
  formula = Surv(agewks.st, agewks.end, event, type=c("interval")) ~ sex + (1 | SubID), 
  data = dat.CF, basehaz = "exp")
ttcf.wei <- update(ttcf.exp, basehaz = "weibull")
ttcf.weiaft <- update(ttcf.exp, basehaz = "weibull-aft")
ttcf.gomp <- update(ttcf.exp, basehaz = "gompertz")
ttcf.ms <- update(ttcf.exp, basehaz = "ms") #default degree is cubic spline

## plot base hazards
exp.haz <- plotfun(ttcf.exp, "Exponential", "Time (weeks)") #exponential has constant hazard rate
wei.haz <- plotfun(ttcf.wei, "Weibull", "Time (weeks)")
weiaft.haz <- plotfun(ttcf.weiaft, "Weibull-aft", "Time (weeks)")
gomp.haz <- plotfun(ttcf.gomp, "Gompertz", "Time (weeks)")
ms.haz <- plotfun(ttcf.ms, "M-spline", "Time (weeks)")
FigureS1 <- plot_grid(exp.haz, wei.haz, weiaft.haz, gomp.haz, ms.haz, ncol = 3)

## leave-one-out cross-validation by individual (instead of by observation)
ids <- dat.CF$SubID
ll.1 <- log_lik(ttcf.wei)
ll.1 <- apply(ll.1, 1L, function(row) tapply(row, ids, sum))
ll.1 <- t(ll.1)
Iloo.wei <- loo(ll.1)

ll.2 <- log_lik(ttcf.weiaft)
ll.2 <- apply(ll.2, 1L, function(row) tapply(row, ids, sum))
ll.2 <- t(ll.2)
Iloo.weiaft <- loo(ll.2)

ll.3 <- log_lik(ttcf.gomp)
ll.3 <- apply(ll.3, 1L, function(row) tapply(row, ids, sum))
ll.3 <- t(ll.3)
Iloo.gomp <- loo(ll.3)

ll.4 <- log_lik(ttcf.ms)
ll.4 <- apply(ll.4, 1L, function(row) tapply(row, ids, sum))
ll.4 <- t(ll.4)
Iloo.ms <- loo(ll.4)

ll.5 <- log_lik(ttcf.exp)
ll.5 <- apply(ll.5, 1L, function(row) tapply(row, ids, sum))
ll.5 <- t(ll.5)
Iloo.exp <- loo(ll.5)

loo_compare(Iloo.exp, Iloo.wei, Iloo.weiaft, Iloo.gomp, Iloo.ms) #LOO
loo_compare(waic(ttcf.exp), waic(ttcf.wei), waic(ttcf.weiaft), waic(ttcf.gomp), waic(ttcf.ms)) #WAIC

##------------------------- time-to-weaning models ---------------------------##
## event times must be >0
dat.wean$agemo.st[dat.wean$agemo.st==0] <- 0.05

## fit exponential, Weibull, Gompertz, m-spline, and Weibull-aft base hazard functions
ttw.exp <- stan_surv(
  formula = Surv(agemo.st, agemo.end, event, type=c("interval")) ~ sex, 
  data = dat.wean, basehaz = "exp")
ttw.wei <- update(ttw.exp, basehaz = "weibull")
ttw.weiaft <- update(ttw.exp, basehaz = "weibull-aft")
ttw.gomp <- update(ttw.exp, basehaz = "gompertz")
ttw.ms <- update(ttw.exp, basehaz = "ms") #default degree is cubic spline

## plot base hazards
wexp.haz <- plotfun(ttw.exp, "Exponential", "Time (months)") 
wwei.haz <- plotfun(ttw.wei, "Weibull", "Time (months)")
wweiaft.haz <- plotfun(ttw.weiaft, "Weibull-aft", "Time (months)")
wgomp.haz <- plotfun(ttw.gomp, "Gompertz", "Time (months)")
wms.haz <- plotfun(ttw.ms, "M-splines", "Time (months)")
FigureS3 <- plot_grid(wexp.haz, wwei.haz, wweiaft.haz, wgomp.haz, wms.haz, ncol = 3)

## leave-one-out cross-validation and WAIC
loo_compare(loo(ttw.exp), loo(ttw.wei), loo(ttw.weiaft), loo(ttw.gomp), loo(ttw.ms)) #LOO
loo_compare(waic(ttw.exp), waic(ttw.wei), waic(ttw.weiaft), waic(ttw.gomp), waic(ttw.ms)) #WAIC
