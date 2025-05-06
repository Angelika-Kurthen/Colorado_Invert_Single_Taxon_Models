#############################################################
## Code to fit curve to Chirnomus spp mortality rates to data
#############################################################

library(readxl)
library(minpack.lm)
library(tidyverse)
#library(car)
library(boot)
CHIRVitalRates <- read_excel("VitalRates.xlsx", sheet = "Chiro Mortality Rates")
CHIRVitalRates <- as.data.frame(CHIRVitalRates)

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

surv.fit.CHIR <- flow.surv.fit(CHIRVitalRates$`Max/Bankful`, CHIRVitalRates$Mortality, 0.25)
surv.df.CHIR <- flow.surv.rate(surv.fit.CHIR$m$getPars()[2] , surv.fit.CHIR$m$getPars()[1], 4, 0.001, 0.001, 0.25)

# Calculate Development Time Based on Literature
# BAETDevRates <- read_excel("VitalRates.xlsx", sheet = "Baetid Development")
# BAETDevRates <- as.data.frame(BAETDevRates)
# polyfit <- nlsLM(DevelopmentTime ~ a*Temperature^4 + b* Temperature ^3 + c*Temperature ^2 +d*Temperature  + e, data = BAETDevRates, start = c(a = 1, b =1, c = 1, d = 1, e = 1))
# 
# devtime <- function(x){
#   a <- 2.525e-03*x^4 -2.508e-01*x^3+  9.379e+00*x^2 -1.580e+02*x +  1.040e+03 
#   return(a)
# }
# 
# MaturationRate <- function(x){
#   a <- 1/x
#   return(a)
# }

# Calculate Temperature Dependent Survival
CHIRSurvRate <- read_excel("VitalRates.xlsx", sheet = "Chiro Survival")
CHIRSurvRate <- as.data.frame(CHIRSurvRate)
#fit <- nlsLM(logit(Survival) ~ a*Temp^4 + b*Temp^3 + c*Temp^2 + d*Temp + e, data = CHIRSurvRate, start = c(a = 1, b = 1, c = 1, d = 1, e = 1))
#fit1 <- nlsLM(logit(Survival) ~ a*Temp^2 + b*Temp + c, data = CHIRSurvRate, start = c(a = 1, b = 1, c = 1))
# inv.logit(predict(fit))
# fit <- nlsLM(Survival ~ a*Temp^2 + b*Temp + c, data = CHIRSurvRate, start = c(a=1, b=1, c=1))
TempSurv_CHIR <- function(n){
    a <- -0.03178*n^2+  1.20308*n -9.25551 
#    #a <- -0.001785*n^2+ 0.074341*n -0.154283
#    #with 0,0
# #    #a <- -0.02184*n^2 +0.92739*n -8.55032
# #  #with 40,0 and 0,0
# #     #a <- -0.02469*n^2 +  1.02303*n -9.10223
# #   # with 40,0
# #    #a <- -0.02187*n^2 +  0.88293*n -7.58722
# #    #wo sankarperumal
# #    #a <- -0.02009*n^2 +0.81130*n -7.02633
# #    #wo reyes malndonald
# #    #a <- -0.01361*n^2 + 0.60686*n -5.79246
# #    #wo eggermont
# #    #a <- -0.01216*n^2 +  0.45588*n -3.53420
# #     # wo stevens
# #    #a <- -0.0241*n^2 + 0.9536*n -7.7931
# #    # wo stratmont
# #    a <- -0.01954*n^2  + 0.81415*n -7.54678 
    return(inv.logit(a))
 }


#a = -1.016e-04  b = 9.412e-03 c = -3.121e-01 d = 4.317e+00 e =-2.032e+01 
# 

#THIS IS THE ONE
TempSurv_CHIR <- function(n){
       a <- -2.164e-04*n^4 + 1.582e-02*n^3 -4.134e-01*n^2+  4.657*n -1.864e+01
#    #a <-  4.011e-05*n^4 -2.260e-03*n^3 -6.379e-04*n^2 + 1.231e+00*n -1.096e+01   #with 0,0
#    #a <- -2.481e-05*n^4 + 2.180e-03*n^3 -8.318e-02*n^2 +  1.494e+00*n -9.197e+00
#    # with 0,0
#    #a <- -4.233e-05*n^4 +3.728e-03*n^3 -1.302e-01*n^2 + 2.055e+00*n -1.127e+01
#    #with 40, 0 and 0,0
#    #a <- -6.733e-05*n^4 + 5.475e-03*n^3 -1.687e-01*n^2 +2.339e+00*n -1.166e+01
#     #with 40,0
#    #a <- -7.434e-05*n^4 +  6.127e-03*n^3 -1.895e-01*n^2  + 2.599e+00*n -1.267e+01
#    #wo sankarperumal
#    #a <- -3.678e-05*n^4  +2.981e-03*n^3 -1.014e-01*n^2 +1.660*n -9.679e+00
#    #wo reyes maldonado
#    #a<- -5.786e-05*n^4 + 5.156e-03*n^3 -1.704e-01*n^2 + 2.451e+00*n -1.230e+01
#    # wo stevens
#   #a <- -1.735e-05*n^4 +  1.254e-03*n^3 -5.321e-02*n^2 + 1.194*n -8.343
#   # wo stratment
#   a <- -3.742e-05*n^4 + 3.902e-03*n^3 -1.566e-01*n^2  + 2.655*n -1.495e+01
  return(inv.logit(a))
       }

tem <- seq(0, 40, by = 1)
tempChir <- as.data.frame(cbind(tem, TempSurv_CHIR(tem)))
# TempSurv <- function(n){
#   #a <- -0.02884*n^2+  1.12797*n -9.38446
#   a <- -1.016e-04*n^4 +  9.412e-03*n^3 -3.121e-01*n^2 + 4.317*n -2.032e+01
#   return(inv.logit(a))
# }
#
# min.RSS <- function(par){
#   mod <- dnbinom(as.integer(CHIRSurvRate$Temp), size = par[2], prob = par[1])
#   a <- sum(CHIRSurvRate$Survival - (mod*(max(CHIRSurvRate$Survival)/max(mod))))^2
# }
# params <- optim(par = c(0.01, 6), fn = min.RSS, method = "BFGS")
# 


#params <- optim(par = c(0, 0), fn = min.RSS, method = "BFGS")
# #
# TempSurv_CHIR <- function(n){
#   if (n <= 0){
#     a <- 0
#   }else {
#     a <-  dnbinom(as.integer(-n + 36), size = params$par[2] , prob = params$par[1])*(max(CHIRSurvRate$Survival)/max(dnbinom(as.integer(-CHIRSurvRate$Temp + 36), size =params$par[2], prob = params$par[1])))
#   }
#   return((a))
# }


# survs <- vector()
# for(i in 1:length(tem)){
#   survs[i] <- TempSurv(tem[i])
# }

# ggplot(surv.df.CHIR, aes(x = Q, y = surv))+
#   geom_line()+
#   geom_point(data = CHIRVitalRates, aes(x = `Max/Bankful` , y = 1-(Mortality), color = Citation))+
#   coord_cartesian(xlim = c(0,4)) +
#   ylab('Immediate Post-Disturbance Survival') +
#   theme_bw()+
#   xlab('`Max Event Discharge/Bankfull Discharge`')

# 
#tem <- seq(0, 40, by = 1)#http://127.0.0.1:23003/graphics/ed0483c1-b811-4223-bfde-c3d3a3dceb62.png
#plot(CHIRSurvRate$Temp, CHIRSurvRate$Survival, col = "red", pch = 16, xlab = "Temperature", ylab = "Survival", xlim = c(0,40), ylim = c(0, 1))
# lines(tem, survs, col = "green")
#lines(tem, TempSurv_CHIR(tem), col = "green")
# #lines(tem,  dnbinom(as.integer(-tem + 37.5), size = params$par[2] , prob = params$par[1])*(max(CHIRSurvRate$Survival)/max(dnbinom(as.integer(-CHIRSurvRate$Temp + 37.5), size =params$par[2], prob = params$par[1]))))
# plot(temp, s, xlab = "Temperature C", ylab = "Survival", col = "red", pch = 16, cex = 1.5, xlim = c(0,40), ylim = c(0,1))
# points(temp, predict(fit.betalogit), col = "blue", pch = 1)
# points(temp, inv.logit(predict(fit4)), col = "hotpink", pch = 2)points(temp, inv.logit(predict(fit2)), col = "green", pch = 5)     
# legend(-0.5, 1, legend=c("Data", "Beta Regression 2D Poly", "Logit 4D Poly", "Logit 2D Poly"),
#        col=c("red", "blue", "hotpink", "green"), pch=c(16, 1, 2, 5), cex=0.8)     
