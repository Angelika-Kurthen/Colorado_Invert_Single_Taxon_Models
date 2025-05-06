############################################################
## Code to fit curve to G. lacustris mortality rates to data
#############################################################

library(readxl)
library(minpack.lm)
library(tidyverse)
#library(car)
library(boot)
GAMMVitalRates <- read_excel("VitalRates.xlsx", sheet = "Gammarus Mortality Rates")
GAMMVitalRates <- as.data.frame(GAMMVitalRates)

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

surv.fit.GAMM <- flow.surv.fit(GAMMVitalRates$`Max/Bankful`, GAMMVitalRates$Mortality, 0.25)
surv.df.GAMM <- flow.surv.rate(surv.fit.GAMM$m$getPars()[2] , surv.fit.GAMM$m$getPars()[1], 4, 0.001, 0.001, 0.25)

# ggplot(surv.df.GAMM, aes(x = Q, y = surv))+
#   geom_line()+
#   geom_point(data = GAMMVitalRates, aes(x = `Max/Bankful` , y = 1-(Mortality), color = Citation))+
#   coord_cartesian(ylim = c(0,1)) +
#   ylab('Immediate Post-Disturbance Survival') +
#   theme_bw()+
#   xlab('`Max Event Discharge/Bankfull Discharge`')

# Calculating Temperature Dependent Mortality
GAMMSurvRates <- read_excel("VitalRates.xlsx", sheet = "Gammarus Survival Rates")
GAMMSurvRates <- as.data.frame(GAMMSurvRates)
fit <- nlsLM(logit(Survival) ~ a*Temperature^4+ b*Temperature^3 + c*Temperature^2 + d*Temperature + e, data = GAMMSurvRates, start = c(a = 1, b = 1, c = 1, d= 1, e = 1))
#fit1 <- nlsLM(logit(Survival)~ a*Temperature^2 + b*Temperature + c, data = GAMMSurvRates, start = c(a =1, b=1, c=1))
TempSurv_GAMM <- function(x){
  #y = -0.05241*x^2 + 1.75693*x -10.20119 
  y = -2.489e-04*x^4 +1.689e-02*x^3 -4.211e-01*x^2 +4.5*x -1.449e+01 
  return(inv.logit(y))
}
tem <- seq(0, 40, by = 1)
tempGamm <- as.data.frame(cbind(tem, TempSurv_GAMM(tem)))

# based on the two fits we have neg binom and inv logit 2nd deg I think the neg inv. binom fits the best, based on 
# 


# tem <- seq(0, 40, by = 1)
# plot(GAMMSurvRates$Temperature, GAMMSurvRates$Survival, col = "red", pch = 16, xlab = "Temperature", ylab = "Survival", xlim = c(0,40), ylim = c(0, 1))
# lines(tem, TempSurv_GAMM(tem))
# 
# 



timestep_to_mat <- function(x){
  # get growth rate in mm/day  - data for Gammarus fossarum, slightly smaller 
  # multiply by 14 to get growth in 2 week intervals
  rate <- (0.0014*x)-0.0024
  r <- rate * 14 
  # at a given temperature, it will take a certain number of time step to 
  stagedur1 <- 4.5/r
  stagedur2 <- 2/r
  return(list(round(stagedur1), round(stagedur2)))
}

  