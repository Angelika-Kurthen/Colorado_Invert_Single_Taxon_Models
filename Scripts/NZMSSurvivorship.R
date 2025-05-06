############################################################
## Code to fit curve to P. anitpodarum mortality rates to data
#############################################################

library(readxl)
library(minpack.lm)
library(tidyverse)
#library(car)
library(boot)
NZMSVitalRates <- read_excel("VitalRates.xlsx", sheet = "NZMS Mortality Rates")
NZMSVitalRates <- as.data.frame(NZMSVitalRates)

flow.surv.fit <- function(magnitude, mortality, Qmin){
  x <- magnitude # rename x and y
  y <- 1 - mortality # convert to survival
  df <- as.data.frame(cbind(x,y))
  df <- rbind(df, c(Qmin, 1.375))  # make sure we specify that we have 100% survival at Qmin (0.25)
  nls.fit <- nlsLM(formula = (y ~ k* exp(-h*x)), data = df, start = c(k = 1, h = 1)) # fit to negative exponential function
  return(nls.fit)
}

flow.surv.rate <- function(h, k, max, min, interval, Qmin) {
  Q <- seq(min, max, by = interval)
  surv <- k*(exp(-h*Q))
  surv.df <- as.data.frame(cbind(Q, surv))
  surv.df$surv[which(surv.df$Q <= Qmin)] <- 1
  return(surv.df)
}

surv.fit.NZMS <- flow.surv.fit(NZMSVitalRates$`Max Event Discharge/Bankfull Discharge`, NZMSVitalRates$Mortality, 0.25)
surv.df.NZMS <- flow.surv.rate(surv.fit.NZMS$m$getPars()[2] , surv.fit.NZMS$m$getPars()[1], 2, 0.001, 0.001, 0.25)

ggplot(surv.df.NZMS, aes(x = Q, y = surv))+
  geom_line()+
  geom_point(data = NZMSVitalRates, aes(x = `Max Event Discharge/Bankfull Discharge` , y = 1-(Mortality), color = Citation))+
  coord_cartesian(ylim = c(0,1)) +
  ylab('Immediate Post-Disturbance Survival') +
  theme_bw()+
  xlab('`Max Event Discharge/Bankfull Discharge`')

# calculating temperature fecundity relationship for NZMS using a power function
# log(y)=log(b*x^z)=log(b)+z*log(x)
# In this case lm(log(y)~log(x) ) solve your problem

# max recruitment between 16 and 19 C (Dybahl and Kane)
# recruitment stops below 9C (Bennett )
# above 27 C everything stops working well (Dybahl and Kane)
# x <- c(9, 16, 17.5, 19, 27)
# y <- c(0.001, 1, 1, 1, 0.001)
# data <- as.data.frame(cbind(x,y))
# fit <- nlsLM(y~ -b*(x - 17.5)^4 + 1,start = list(b = 0),data=data)

# eq = y = -0.0001427(x - 17.5)^4 + 1 * F 

# Calculating Temperature Dependent Mortality
NZMSSurvRates <- read_excel("VitalRates.xlsx", sheet = "NZMS Survival Rates")
NZMSSurvRates <- as.data.frame(NZMSSurvRates)
fit <- nlsLM(logit(Survival) ~ a*Temperature^2+ b*Temperature + c, data = NZMSSurvRates, start = c(a = 1, b = 1, c = 1))
TempSurv_NZMS <- function(x){
  y = -0.08814*x^2  +3.09981*x -9.18655 
  return(inv.logit(y))
}


tem <- seq(0, 40, by = 1)
tempNZMS <- as.data.frame(cbind(TempSurv_NZMS(tem)))

# based on the two fits we have neg binom and inv logit 2nd deg I think the neg inv. binom fits the best, based on 
# 
# min.RSS <- function(par){
#   mod <- dnbinom(as.integer(-NZMSSurvRates$Temperature + 34), size = par[2], prob = par[1])
#   a <- sum(NZMSSurvRates$Survival - (mod*(max(NZMSSurvRates$Survival)/max(mod))))^2
# }
#  params <- optim(par = c(0.2, 2), fn = min.RSS, method = "BFGS")
# # 
# TempSurv <- function(x){
#   a <-  dnbinom(-x + 34, size = params$par[2] , prob = params$par[1])*(max(NZMSSurvRates$Survival)/max(dnbinom(as.integer(-NZMSSurvRates$Temperature + 34), size =params$par[2], prob = params$par[1])))
#   if (x <= 0){
#     a <- 0
#   }
#   return((a))
# }
# 
# TempSurv <- function(x){
#   a <-  dnbinom(-x + 32, size = 3, prob = 0.2)
#   return((a))
# }
# 
# tem <- seq(0, 40, by = 1)
# plot(NZMSSurvRates$Temperature, NZMSSurvRates$Survival, col = "red", pch = 16, xlab = "Temperature", ylab = "Survival", xlim = c(0,40), ylim = c(0, 1))
# lines(tem, TempSurv_NZMS(tem))

# 
