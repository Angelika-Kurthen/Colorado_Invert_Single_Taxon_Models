#####################################
## Code to fit curve to Hydropsyche spp mortality rates to data
####################################
#Code for HPC
# library(readxl, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
# library(minpack.lm, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
# library(readxl, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
# library(car, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
# library(boot, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
# library(data.table, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
# library(tibble, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
# library(tidyr, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
# library(dplyr, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
# library(purrr, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
# library(forcats, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
# library(stringr, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
# library(readr, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
# 



library(readxl)
library(minpack.lm)
library(tidyverse)
#library(car)
library(boot)
library(data.table)
HYOSVitalRates <- read_excel("VitalRates.xlsx", sheet = "Hydropsyche Mortality Rates")
HYOSVitalRates <- as.data.frame(HYOSVitalRates)

flow.surv.fit <- function(magnitude, mortality, Qmin){
  x <- magnitude # rename x and y
  y <- 1 - mortality # convert to survival
  df <- as.data.frame(cbind(x,y))
  df <- rbind(df, c(Qmin, 1))  # make sure we specify that we have 100% survival at Qmin (0.25)
  nls.fit <- nlsLM(formula = (y ~ k* exp(-h*x)), data = df, start = c(k = 1, h = 1)) # fit to negative exponential function
  return(nls.fit)
}

flow.surv.rate <- function(h, k, max, min, interval, Qmin) {
  Q <- seq(min, max, by = interval)
  surv <- k*exp(-h*Q)
  surv.df <- as.data.frame(cbind(Q, surv))
  surv.df$surv[which(surv.df$Q <= Qmin)] <- 1
  return(surv.df)
}

surv.fit.HYOS <- flow.surv.fit(HYOSVitalRates$`Max Event Discharge/Bankfull Discharge`, HYOSVitalRates$Mortality, 0.15)
surv.df.HYOS <- flow.surv.rate(surv.fit.HYOS$m$getPars()[2] , surv.fit.HYOS$m$getPars()[1], 2, 0.001, 0.001, 0.2)

ggplot(surv.df.HYOS, aes(x = Q, y = surv))+
  geom_line()+
  geom_point(data = HYOSVitalRates, aes(x = `Max Event Discharge/Bankfull Discharge` , y = 1-(Mortality), color = Citation))+
  #coord_cartesian(ylim = c(0,1)) +
  ylab('Immediate Post-Disturbance Survival') +
  theme_bw()+
  xlab('`Max Event Discharge/Bankfull Discharge`')

# Calculate temperature dependent development time
HYOSDevRates <- read_excel("VitalRates.xlsx", sheet = "Hydropsyche Development Rates")
HYOSDevRates <- as.data.frame(HYOSDevRates)

polyfit <- nlsLM(logit(MaturationRate) ~ a*Temperature^2 + b*Temperature + c, data = HYOSDevRates, start = c(a = 1, b = 1, c = 1))

devtime <- function(x){
  y = -0.01385  *x^2+ 0.30973*x -5.72982
  return(1/inv.logit(y))
}
# 
HYOSSurvRates <- read_excel("VitalRates.xlsx", sheet = "Hydropsyche Survival Rates ")
HYOSSurvRates <- as.data.frame(HYOSSurvRates)
# fit <- nlsLM(logit(Survival)~ a*Temperature^2 + b*Temperature + c, data = HYOSSurvRates, start = c(a= 1, b=1, c = 1))
# 
TempSurv_HYOS <- function(x){
  a <- -0.09934*x^2 +3.44127*x -15.47038
  return(inv.logit(a))
}
tem <- seq(0, 40, by = 1)
tempHyos <- as.data.frame(cbind(tem,TempSurv_HYOS(tem)))
# min.RSS <- function(par){
#   mod <- dnbinom(as.integer(-HYOSSurvRates$Temperature + 31), size = par[2], prob = par[1])
#   a <- sum(HYOSSurvRates$Survival - (mod*(max(HYOSSurvRates$Survival)/max(mod))))^2
# }
# params <- optim(par = c(0.23, 4.5), fn = min.RSS)
# 
# TempSurv_HYOS <- function(n){
#   if (n <= 0){
#     a <- 0
#   }else{
#   a <-  dnbinom(as.integer(-n + 32), size = 2.8835371 , prob = 0.1932115)*(max(HYOSSurvRates$Survival)/max(dnbinom(as.integer(-HYOSSurvRates$Temperature + 32), size =2.8835371, prob = 0.1932115 )))
#   }
#   return((a))
# }


# 
tem <- seq(0, 40, by = 1)
# temSurv <- unlist(lapply(tem, TempSurv_HYOS))
# plot(HYOSSurvRates$Temperature, HYOSSurvRates$Survival, col = "red", pch = 16, xlab = "Temperature", ylab = "Survival", xlim = c(0,40), ylim = c(0, 1))
# lines(tem,  unlist(lapply(tem, TempSurv_HYOS)))
# 
