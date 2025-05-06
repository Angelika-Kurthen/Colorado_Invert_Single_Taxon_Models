#############################################################
## Code to fit curve to Baetidae spp mortality rates to data
#############################################################

library(readxl)
library(minpack.lm)
library(tidyverse)
# library(car)
library(boot)

BAETVitalRates <- read_excel("VitalRates.xlsx", sheet = "Baetid Mortality Rates")
BAETVitalRates <- as.data.frame(BAETVitalRates)

flow.surv.fit <- function(magnitude, mortality, Qmin){
  x <- magnitude # rename x and y
  y <- 1 - mortality # convert to survival
  df <- as.data.frame(cbind(x,y))
  df <- rbind(df, c(Qmin, 1.11))  # make sure we specify that we have 100% survival at Qmin (0.25)
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

surv.fit.BAET <- flow.surv.fit(BAETVitalRates$`Max Event Discharge/Bankfull Discharge`, BAETVitalRates$Mortality, 0.25)
surv.df.BAET <- flow.surv.rate(surv.fit.BAET$m$getPars()[2] , surv.fit.BAET$m$getPars()[1], 4, 0.001, 0.001, 0.25)

# Calculate Development Time Based on Literature
BAETDevRates <- read_excel("VitalRates.xlsx", sheet = "Baetid Development")
BAETDevRates <- as.data.frame(BAETDevRates)
polyfit <- nlsLM(DevelopmentTime ~ a*Temperature^4 + b* Temperature ^3 + c*Temperature ^2 +d*Temperature  + e, data = BAETDevRates, start = c(a = 1, b =1, c = 1, d = 1, e = 1))

devtime <- function(x){
  a <- 2.525e-03*x^4 -2.508e-01*x^3+  9.379e+00*x^2 -1.580e+02*x +  1.040e+03 
  return(a)
}

MaturationRate <- function(x){
  a <- 1/x
  return(a)
}

# Calculate Temperature Dependent Survival
BAETSurvRate <- read_excel("VitalRates.xlsx", sheet = "Baetid Survival Rates")
BAETSurvRate <- as.data.frame(BAETSurvRate)
#fit <- nlsLM(logit(Survival) ~ a*Temperature^4 + b*Temperature^3 + c*Temperature^2 + d*Temperature + e, data = BAETSurvRate, start = c(a = 1, b = 1, c = 1, d = 1, e = 1))
#inv.logit(predict(fit))
#-0.02837*x^2 + 1.21299*x  -10.92723

#fit <- nlsLM(logit(Survival) ~ a*Temperature^2 + b*Temperature+ c, data = BAETSurvRate, start = c(a = 1, b = 1, c = 1))
# # 
TempSurv_BAET <- function(n){
  a <-  -0.02429*n^2 +0.73459 *n -2.59225 
  return(inv.logit(a))
}

tem <- seq(0, 40, by = 1)
tempBaet <- as.data.frame(cbind(tem, TempSurv_BAET(tem)))
# TempSurv_BAET <- function(n){
#   a <- -5.782e-04*n^4+ 4.425e-02*n^3 -1.200*n^2+ 1.348e+01*n -4.991e+01 
#   return(inv.logit(a))
# }
# # 
# min.RSS <- function(par){
#   mod <- dnbinom(as.integer(-BAETSurvRate$Temperature + 34), size = par[2], prob = par[1])
#   a <- sum(BAETSurvRate$Survival - (mod*(max(BAETSurvRate$Survival)/max(mod))))^2
# }
# params <- optim(par = c(0.2, 4), fn = min.RSS, method = "BFGS")
# 
# TempSurv_BAET <- function(n){
#   if (n <= 0){
#     a <- 0
#   }else {
#   a <-  dnbinom(as.integer(-n + 34), size = params$par[2] , prob = params$par[1])*(max(BAETSurvRate$Survival)/max(dnbinom(as.integer(-BAETSurvRate$Temperature + 34), size =params$par[2], prob = params$par[1])))
#   }
#    return((a))
# }


# ggplot(surv.df.BAET, aes(x = Q, y = surv))+
#   geom_line()+
#   geom_point(data = BAETVitalRates, aes(x = `Max Event Discharge/Bankfull Discharge` , y = 1-(Mortality), color = Citation))+
#   #coord_cartesian(ylim = c(0,1)) +
#   ylab('Immediate Post-Disturbance Survival') +
#   theme_bw()+
#   xlab('`Max Event Discharge/Bankfull Discharge`')
# 
# 
# tem <- seq(0, 40, by = 1)
# plot(BAETSurvRate$Temperature, BAETSurvRate$Survival, col = "red", pch = 16, xlab = "Temperature", ylab = "Survival", xlim = c(0,40), ylim = c(0, 1))
# lines(tem, TempSurv_BAET(tem))
#lines(tem,  dnbinom(as.integer(-tem + 33), size = params$par[2] , prob = params$par[1])*(max(BAETSurvRate$Survival)/max(dnbinom(as.integer(-BAETSurvRate$Temperature + 33), size =params$par[2], prob = params$par[1]))))
# plot(temp, s, xlab = "Temperature C", ylab = "Survival", col = "red", pch = 16, cex = 1.5, xlim = c(0,40), ylim = c(0,1))
# points(temp, predict(fit.betalogit), col = "blue", pch = 1)
# points(temp, inv.logit(predict(fit4)), col = "hotpink", pch = 2)
# points(temp, inv.logit(predict(fit2)), col = "green", pch = 5)     
# legend(-0.5, 1, legend=c("Data", "Beta Regression 2D Poly", "Logit 4D Poly", "Logit 2D Poly"),
#        col=c("red", "blue", "hotpink", "green"), pch=c(16, 1, 2, 5), cex=0.8)     
