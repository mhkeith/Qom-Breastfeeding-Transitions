#### Qom breastfeeding transitions ####
#### Time-to-event analysis with gestational covariates ####

## load packages
library(tidyverse)
library(rstanarm)
library(loo)
library(cowplot)
library(ggfortify)
library(survival)

## load data
dat.CF <- read.csv(file="~/QomTTCF.dat.csv") #time-to-complementary feeding, repeated observations 
dat.wean <- read.csv(file="~/QomTTWean.dat.csv") #time-to-weaning, single observations

##------------------- time-to-complementary feeding models -------------------##
## covariate screening

## subset to exclude missing covariate data
sub.CF <- subset(dat.CF, 
                 !(is.na(csect) | is.na(primipar2) | is.na(matage.cat) | is.na(birthweight) | is.na(gestage))) 
## 89 infants in this subset, 253 repeated observations
sub.CF$event[sub.CF$leftcensored==1] <- 2 #code left-censored events with Surv object notation
## event times must be >0 since log(t) is used in parametric equations
sub.CF$agewks.st[sub.CF$agewks.st==0] <- 0.5

## base model with sex as only covariate  
base.cf <- stan_surv(
  formula = Surv(agewks.st, agewks.end, event, type=c("interval")) ~ sex + (1 | SubID), 
  data = sub.CF, basehaz = "gompertz", #best fitting base hazard from screening
  chains = 4, cores = 2, iter = 4000)
print(base.cf)  

## birth mode model
csect.cf <- stan_surv(
  formula = Surv(agewks.st, agewks.end, event, type=c("interval")) ~ sex + csect + (1 | SubID), 
  data = sub.CF, basehaz = "gompertz",
  chains = 4, cores = 2, iter = 4000)
print(csect.cf)

## parity model
par.cf <- stan_surv(
  formula = Surv(agewks.st, agewks.end, event, type=c("interval")) ~ sex + primipar2 + (1 | SubID),
  data = sub.CF, basehaz = "gompertz",
  chains = 4, cores = 2, iter = 4000)
print(par.cf)

## maternal age model
matage.cf <- stan_surv(
  formula = Surv(agewks.st, agewks.end, event, type=c("interval")) ~ sex + matage.cat + (1 | SubID),
  data = sub.CF, basehaz = "gompertz",
  chains = 4, cores = 2, iter = 4000)
print(matage.cf)

## birthweight, continuous model
bwcont.cf <- stan_surv(
  formula = Surv(agewks.st, agewks.end, event, type=c("interval")) ~ sex + birthweight + (1 | SubID),
  data = sub.CF, basehaz = "gompertz",
  chains = 4, cores = 2, iter = 4000)
print(bwcont.cf)

## birthweight, 3 categories model
bw3cat.cf <- stan_surv(
  formula = Surv(agewks.st, agewks.end, event, type=c("interval")) ~ sex + bw.cat3 + (1 | SubID),
  data = sub.CF, basehaz = "gompertz",
  chains = 4, cores = 2, iter = 4000)
print(bw3cat.cf)

## birthweight, 2 categories model
bw2cat.cf <- stan_surv(
  formula = Surv(agewks.st, agewks.end, event, type=c("interval")) ~ sex + bw.cat2 + (1 | SubID),
  data = sub.CF, basehaz = "gompertz",
  chains = 4, cores = 2, iter = 4000)
print(bw2cat.cf)

## gestational age, continuous model
gestcont.cf <- stan_surv(
  formula = Surv(agewks.st, agewks.end, event, type=c("interval")) ~ sex + gestage + (1 | SubID),
  data = sub.CF, basehaz = "gompertz",
  chains = 4, cores = 2, iter = 4000)
print(gestcont.cf)

## gestational age, 4 categories model
gest4.cf <- stan_surv(
  formula = Surv(agewks.st, agewks.end, event, type=c("interval")) ~ sex + gest.cat4 + (1 | SubID),
  data = sub.CF, basehaz = "gompertz",
  chains = 4, cores = 2, iter = 4000)
print(gest4.cf)

## gestational age, 2 categories model
sub.CF$gest.2early <- ifelse(sub.CF$gest.cat2=="early", 1, 0) #make full term the reference category
gest2.cf <- stan_surv(
  formula = Surv(agewks.st, agewks.end, event, type=c("interval")) ~ sex + gest.2early + (1 | SubID),
  data = sub.CF, basehaz = "gompertz",
  chains = 4, cores = 2, iter = 4000)
print(gest2.cf)

## leave-one-out cross-validation by individual

ids <- sub.CF$SubID
ll.b <- log_lik(base.cf)
ll.b <- apply(ll.b, 1L, function(row) tapply(row, ids, sum))
ll.b <- t(ll.b)
loo.base <- loo(ll.b)

ll.cs <- log_lik(csect.cf)
ll.cs <- apply(ll.cs, 1L, function(row) tapply(row, ids, sum))
ll.cs <- t(ll.cs)
loo.cs <- loo(ll.cs)

ll.par <- log_lik(par.cf)
ll.par <- apply(ll.par, 1L, function(row) tapply(row, ids, sum))
ll.par <- t(ll.par)
loo.par <- loo(ll.par) 

ll.ma <- log_lik(matage.cf)
ll.ma <- apply(ll.ma, 1L, function(row) tapply(row, ids, sum))
ll.ma <- t(ll.ma)
loo.ma <- loo(ll.ma)

ll.bw <- log_lik(bwcont.cf)
ll.bw <- apply(ll.bw, 1L, function(row) tapply(row, ids, sum))
ll.bw <- t(ll.bw)
loo.bw <- loo(ll.bw)

ll.bw3 <- log_lik(bw3cat.cf)
ll.bw3 <- apply(ll.bw3, 1L, function(row) tapply(row, ids, sum))
ll.bw3 <- t(ll.bw3)
loo.bw3 <- loo(ll.bw3)

ll.bw2 <- log_lik(bw2cat.cf)
ll.bw2 <- apply(ll.bw2, 1L, function(row) tapply(row, ids, sum))
ll.bw2 <- t(ll.bw2)
loo.bw2 <- loo(ll.bw2)

ll.gest <- log_lik(gestcont.cf)
ll.gest <- apply(ll.gest, 1L, function(row) tapply(row, ids, sum))
ll.gest <- t(ll.gest)
loo.gest <- loo(ll.gest)

ll.gest4 <- log_lik(gest4.cf)
ll.gest4 <- apply(ll.gest4, 1L, function(row) tapply(row, ids, sum))
ll.gest4 <- t(ll.gest4)
loo.gest4 <- loo(ll.gest4)

ll.gest2 <- log_lik(gest2.cf)
ll.gest2 <- apply(ll.gest2, 1L, function(row) tapply(row, ids, sum))
ll.gest2 <- t(ll.gest2)
loo.gest2 <- loo(ll.gest2)

loo_compare(loo.base, loo.cs, loo.par, loo.ma, loo.bw, loo.bw3, loo.bw2, loo.gest, loo.gest4, loo.gest2) #LOO
loo_compare(waic(base.cf), waic(csect.cf), waic(par.cf), waic(matage.cf), waic(bwcont.cf), 
            waic(bw3cat.cf), waic(bw2cat.cf), waic(gestcont.cf), waic(gest4.cf), waic(gest2.cf)) #WAIC

## full time-to-complementary feeding model with gestational covariates  
fullmod.cf <- stan_surv(
  formula = Surv(agewks.st, agewks.end, event, type=c("interval")) ~
    sex + csect + matage.cat + gestage + (1 | SubID), #best fitting covariates from screening
  data = sub.CF, basehaz = "gompertz", iter = 4000)
print(fullmod.cf)  
summary(fullmod.cf)
plot(fullmod.cf, "trace", par = "(Intercept)") #trace plots

##------------------------ time-to-weaning models ----------------------------##
#### covariate screening models 

dat.wean$agemo.st[dat.wean$agemo.st==0]<-0.05 #log(t) requires t>0
## subset to exclude missing covariate data
sub.wean <- subset(dat.wean, !(is.na(csect) | is.na(primipar2) | is.na(matage.cat) | is.na(birthweight) | is.na(gestage)))

## base model with sex as only covariate  
base.wean <- stan_surv(
  formula = Surv(agemo.st, agemo.end, event, type=c("interval")) ~ sex, 
  data = sub.wean, basehaz = "weibull", #best fitting base hazard from screening
  chains = 4, cores = 2, iter = 4000)
print(base.wean)  

## birth mode
csect.wean <- stan_surv(
  formula = Surv(agemo.st, agemo.end, event, type=c("interval")) ~ sex + csect, 
  data = sub.wean, basehaz = "weibull",
  chains = 4, cores = 2, iter = 4000)
print(csect.wean)  

## parity 
par.wean <- stan_surv(
  formula = Surv(agemo.st, agemo.end, event, type=c("interval")) ~ sex + primipar2, 
  data = sub.wean, basehaz = "weibull",
  chains = 4, cores = 2, iter = 4000)
print(par.wean)

## maternal age 
matage.wean <- stan_surv(
  formula = Surv(agemo.st, agemo.end, event, type=c("interval")) ~ sex + matage.cat, 
  data = sub.wean, basehaz = "weibull",
  chains = 4, cores = 2, iter = 4000)
print(matage.wean)

## birthweight, continuous
bwcont.wean <- stan_surv(
  formula = Surv(agemo.st, agemo.end, event, type=c("interval")) ~ sex + birthweight, 
  data = sub.wean, basehaz = "weibull",
  chains = 4, cores = 2, iter = 4000)
print(bwcont.wean)

## birthweight, 3 categories
bw3cat.wean <-stan_surv(
  formula = Surv(agemo.st, agemo.end, event, type=c("interval")) ~ sex + bw.cat3, 
  data = sub.wean, basehaz = "weibull",
  chains = 4, cores = 2, iter = 4000)
print(bw3cat.wean)

## birthweight, 2 categories
bw2cat.wean <- stan_surv(
  formula = Surv(agemo.st, agemo.end, event, type=c("interval")) ~ sex + bw.cat2, 
  data = sub.weant, basehaz = "weibull",
  chains = 4, cores = 2, iter = 4000)
print(bw2cat.wean)

## gestational age, continuous
gestcont.wean <- stan_surv(
  formula = Surv(agemo.st, agemo.end, event, type=c("interval")) ~ sex + gestage, 
  data = sub.wean, basehaz = "weibull",
  chains = 4, cores = 2, iter = 4000)
print(gestcont.wean)

## gestational age, 4 categories
gest4.wean <- stan_surv(
  formula = Surv(agemo.st, agemo.end, event, type=c("interval")) ~ sex + gest.cat4, 
  data = sub.wean, basehaz = "weibull",
  chains = 4, cores = 2, iter = 4000)
print(gest4.wean)

## gestational age, 2 categories
sub.wean$gest.2early <- ifelse(sub.wean$gest.cat2=="early", 1, 0) #make full-term the reference
gest2.wean <- stan_surv(
  formula = Surv(agemo.st, agemo.end, event, type=c("interval")) ~ sex + gest.2early, 
  data = sub.wean, basehaz = "weibull",
  chains = 4, cores = 2, iter = 4000)
print(gest2.wean)

## leave-one-out cross-validation
loo_compare(loo(base.wean), loo(csect.wean), loo(par.wean), loo(matage.wean), loo(bwcont.wean), 
            loo(bw3cat.wean), loo(bw2cat.wean), loo(gestcont.wean), loo(gest4.wean), loo(gest2.wean))
## WAIC
loo_compare(waic(base.wean), waic(csect.wean), waic(par.wean), waic(matage.wean), waic(bwcont.wean),
            waic(bw3cat.wean), waic(bw2cat.wean), waic(gestcont.wean), waic(gest4.wean), waic(gest2.wean))

## age at complementary feeding model
## down-sample to infants with known age at complementary feeding
wean.cfsub <- subset(sub.wean, !(is.na(sub.wean$ageCF.mo))) #58 infants
agecf.wean <- stan_surv(
  formula = Surv(agemo.st, agemo.end, event, type=c("interval")) ~ sex + ageCF.mo, 
  data = wean.cfsub, basehaz = "weibull",
  chains = 4, cores = 2, iter = 4000)
print(agecf.wean)

##full time-to-weaning model with gestational covariates  
fullmod.wean <- stan_surv(
  formula = Surv(agemo.st, agemo.end, event, type=c("interval")) ~
  sex + csect + primipar2 + gest.cat2, #best fitting covariates from screening
  data = sub.wean, basehaz = "weibull", iter = 4000)
print(fullmod.wean)  
summary(fullmod.wean)
plot(fullmod.wean, "trace", par = "(Intercept)") #trace plots


##-------------------------------- figures -----------------------------------##
## Figure 1
## time-to-complementary feeding curve
ps.cf <- ps_check(base.cf) 
ttcf <- ps.cf + ggplot2::theme_classic() + 
  xlab("Time (weeks)") + ggtitle("Time-to-complementary feeding") + 
  geom_hline(yintercept = 0.5, color="black", linetype="dashed")
## time-to-weaning curve
ps.w <- ps_check(base.wean)
ttw <- ps.w + ggplot2::theme_classic() + 
  xlab("Time (months)") + ggtitle("Time-to-weaning") + xlim(0,60) + 
  geom_hline(yintercept = 0.5, color="black", linetype="dashed")
## plot together
Figure1 <- plot_grid(ttcf, ttw, ncol = 2)

## Figure 2
## hazard ratio estimates from full time-to-complementary feeding model  
coef.fullcf <- data.frame(fullmod.cf$stan_summary[2:6, c("mean", "2.5%", "97.5%")])
rownames(coef.fullcf) <- c("Male", "C-section", "Maternal age <20", "Maternal age 30+", "Gestational age")
coef.fullcf[ ,4:6] <- exp(coef.fullcf[ ,1:3]) #transform coefficients to hazard ratios 
colnames(coef.fullcf) <- c("post.mean", "2.5%", "97.5%", "haz.mean", "haz.low", "haz.high")
## forest plot of covariate hazard ratios (means and 95% credible intervals) 
coef.fullcf$Covariate <- factor(row.names(coef.fullcf))
hazplot.fullcf <- ggplot(data = coef.fullcf, aes(x = Covariate, y = haz.mean, ymin = haz.low, ymax = haz.high)) 
fullcf <- hazplot.fullcf +  geom_point(size=3,position = position_dodge(width = 0.7)) + 
  scale_y_continuous(trans="log10") +
  geom_errorbar(position = position_dodge(width = 0.5), width = 0.3) + ylab("Hazard Ratio (log scale)") + 
  xlab(NULL) + ggtitle("Time-to-complementary feeding covariates") + 
  theme_classic() + geom_hline(yintercept=1.0) + coord_flip()

## hazard ratio estimates from the full time-to-weaning model  
coef.fullw <- data.frame(fullmod.wean$stan_summary[2:5, c("mean", "2.5%", "97.5%")])
rownames(coef.fullw) <- c("Male", "C-section", "Primiparous", "Early term")
coef.fullw[ ,4:6] <- exp(coef.fullw[ ,1:3]) #transform coefficients to hazard ratios 
colnames(coef.fullw) <- c("post.mean", "2.5%", "97.5%", "haz.mean", "haz.low", "haz.high")
## forest plot of covariate hazard ratios (means and 95% credible intervals) 
coef.fullw$Covariate <- factor(row.names(coef.fullw))
hazplot.fullw <- ggplot(data = coef.fullw, aes(x = Covariate, y = haz.mean, ymin = haz.low, ymax = haz.high)) 
fullw <- hazplot.fullw +  geom_point(size=3,position = position_dodge(width = 0.7)) + 
  scale_y_continuous(trans="log10", limits = c(0.2, 3.5)) +
  geom_errorbar(position = position_dodge(width = 0.5), width = 0.3) + ylab("Hazard Ratio (log scale)") + 
  xlab(NULL) + ggtitle("Time-to-weaning covariates") + theme_classic() + geom_hline(yintercept=1.0) + coord_flip()

Figure2 <- plot_grid(fullcf, fullw, ncol = 1, scale=0.90)

## Figure 3
## birth mode-adjusted time-to-weaning curves
sub.wean$csect[sub.wean$csect==0] <- "vaginal"
sub.wean$csect[sub.wean$csect==1] <- "c-section"
fit.csect <- survfit(Surv(agemo.st, agemo.end, event, type="interval") ~ csect, data = sub.wean)
csect.plot <- autoplot(fit.csect, conf.int = FALSE, xlab="Time (months)", ylab="Percent breastfeeding", censor = FALSE) 
Figure3 <- csect.plot + theme_classic() + guides(fill=FALSE) + 
  labs(colour = "Birth Mode") + scale_color_manual(values=c("black", "gray")) +
  ggtitle("Birth mode-adjusted time-to-weaning") + geom_hline(yintercept = 0.5, linetype="dashed") + xlim(0,60) +  
  theme(legend.position = c(0.85,0.85), legend.background = element_rect(fill = "white", linetype="solid", color="black"))

## Figure S2
## sex-adjusted time-to-complementary feeding curves
fit.cf.sex <- survfit(Surv(agewks.st, agewks.end, event, type="interval") ~ sex, data = sub.CF)
sex.plot.cf <- autoplot(fit.cf.sex, conf.int = FALSE, xlab="Time (weeks)", ylab="Percent excl. breastfeeding", censor = FALSE) 
FigureS2 <- sex.plot.cf + theme_classic() + guides(fill=FALSE) + labs(colour = "Sex") +
  ggtitle("Sex-adjusted time-to-complementary feeding") + geom_hline(yintercept = 0.5, linetype="dashed") +
  scale_color_manual(values=c("gray", "black")) +  
  theme(legend.position = c(0.85,0.85), legend.background = element_rect(fill = "white", linetype="solid", color="black"))

## Figure S4
## gestational age-adjusted time-to-weaning curves
fit.gest <- survfit(Surv(agemo.st, agemo.end, event, type="interval") ~ gest.cat2, data = sub.wean)
gest.plot <- autoplot(fit.gest, conf.int = FALSE, xlab="Time (months)", ylab="Percent breastfeeding", censor = FALSE) 
FigureS4 <- gest.plot + theme_classic() + guides(fill=FALSE) + 
  xlim(0,60) + labs(colour = "Gestational Age") + ggtitle("Early/Full term-adjusted time-to-weaning") + 
  geom_hline(yintercept = 0.5, linetype="dashed") + scale_color_manual(values=c("gray", "black")) +  
  theme(legend.position = c(0.85,0.85), legend.background = element_rect(fill = "white", linetype="solid", color="black"))

## Figure S5
## sex-adjusted time-to-weaning curves
fit.wsex <- survfit(Surv(agemo.st, agemo.end, event, type="interval") ~ sex, data = sub.wean)
wsex.plot <- autoplot(fit.wsex, conf.int = FALSE, xlab="Time (months)", ylab="Percent breastfeeding", censor = FALSE) 
FigureS5 <- wsex.plot + theme_classic() + guides(fill=FALSE) + labs(colour = "Sex") + 
  ggtitle("Sex-adjusted time-to-weaning") + xlim(0,60) + geom_hline(yintercept = 0.5, linetype="dashed") + 
  scale_color_manual(values=c("gray", "black")) +  
  theme(legend.position = c(0.85,0.85), legend.background = element_rect(fill = "white", linetype="solid", color="black"))
