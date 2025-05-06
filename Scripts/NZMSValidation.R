#########################
## Code to Validate NZMS model
###########################
library(purrr)
library(tidyverse)
library(lubridate)
library(plyr)
library(dplyr)
library(ggplot2)
# data retrieval tool from USGS
library(dataRetrieval)
#install.packages("devtools")
library(devtools)
#install_github(repo = "jmuehlbauer-usgs/R-packages", subdir = "foodbase")
library(foodbase)


source("NZMS_1sp_Model.R")

discharge <- readNWISdv("09380000", "00060", "2007-10-01", "2023-05-01")
flow.magnitude <- TimestepDischarge(discharge, 85000)
temp <- readNWISdv("09380000", "00010", "2007-10-01", "2023-05-01")
temps <- TimestepTemperature(temp)

out <- NZMSmodel(flow.data = flow.magnitude$Discharge, temp.data = temps, disturbanceK = 1000000, baselineK = 5000, Qmin = 0.25, extinct = 50, iteration = 9, peaklist = 0.17, peakeach = length(temps$Temperature))


# adults<-as.data.frame(cbind(as.Date(temps$dts), out[1:length(temps$dts),2:3,1]))
# colnames(adults) <- c("Time","Adult")
# adults$Time <- as.Date(adults$Time, origin = "1970-01-01")

means.list.NZMS <- mean.data.frame(out,burnin = 50, iteration= 9)
means.list.NZMS <- cbind(means.list.NZMS, temps$dts[1:length(means.list.NZMS$timesteps)])
means.list.NZMS$`temps$dts` <- as.Date(means.list.NZMS$`temps$dts`)

# culling data by lags - Temporal autocorrelation: a neglected factor in the study of behavioral repeatability and plasticity  

source(NMix_HPC.R)


acf(estimate$adjust)

cor.df <- left_join(estimate, means.list.NZMS, by=c('date'="temps$dts"), copy = T)
cor.lm <- lm((cor.df$mean.abund) ~ (cor.df$est_flow))
cor.test((cor.df$est_flow), (cor.df$mean.abund+1), method = "spearman")

NZMS.samp.sum1 <- cor.df %>% slice(which(row_number() %% 3 == 0))
NZMS.samp.sum2 <- cor.df %>%  slice(which(row_number() %% 3 == 1))
NZMS.samp.sum3 <- cor.df %>%  slice(which(row_number() %%  3 == 2))
  
cor.df1 <- left_join(NZMS.samp.sum1, means.list.NZMS, by=c('date'="temps$dts"), copy = T)
cor.lm1 <- lm(cor.df1$mean.abund.x ~ cor.df1$est_flow)
summary(cor.lm1)
plot(cor.df1$V2, cor.df1$mean.abund)
rho1 <- cor.test((cor.df1$est_flow), (cor.df1$mean.abund.x), method = "spearman")

cor.df2 <- left_join(NZMS.samp.sum2, means.list.NZMS, by=c('date'="temps$dts"), copy = T)
cor.lm2 <- lm(cor.df2$mean.abund.x ~ cor.df2$est_flow)
rho2 <- cor.test((cor.df2$est_flow+1), (cor.df2$mean.abund.x+1), method = "spearman")

cor.df3 <- left_join(NZMS.samp.sum3, means.list.NZMS, by=c('date'="temps$dts"), copy = T)
cor.lm3 <- lm(cor.df3$mean.abund.x ~ cor.df3$est_flow)
summary(cor.lm3)
rho3 <- cor.test((cor.df3$est_flow+1), (cor.df3$mean.abund.x+1), method = "spearman")

rho <- mean(c(rho1$estimate, rho2$estimate, rho3$estimate))
sd <- sd(c(rho1$estimate, rho2$estimate, rho3$estimate))

summary(cor.df)
ggplot(data = cor.df, aes(x = (est) , y = (mean.abund)))+
  geom_point()+
  stat_smooth(method = "lm",
              formula = y ~ x,
              geom = "smooth")+
  geom_text(x = 3e+05, y = 6100, label = "y = 0.0047x, R^2 = 0.05")+
  labs(y = "NZMS Model Output", x = "NZMS Emprical Data")
  

#flow
cor.df <- left_join(estimate, means.list.NZMS, by=c('date'="temps$dts"), copy = T)
cor.lm <- lm((cor.df$mean.abund) ~ (cor.df$adjust))
cor.test((cor.df$adjust), (cor.df$mean.abund+1), method = "spearman")

NZMS.samp.sum1 <- cor.df %>% slice(which(row_number() %% 3 == 0))
NZMS.samp.sum2 <- cor.df %>%  slice(which(row_number() %% 3 == 1))
NZMS.samp.sum3 <- cor.df %>%  slice(which(row_number() %%  3 == 2))

cor.df1 <- left_join(NZMS.samp.sum1, means.list.NZMS, by=c('date'="temps$dts"), copy = T)
cor.lm1 <- lm(cor.df1$mean.abund.x ~ cor.df1$adjust)
summary(cor.lm1)
plot(cor.df1$V2, cor.df1$mean.abund)
rho1 <- cor.test((cor.df1$adjust), (cor.df1$mean.abund.x), method = "spearman")

cor.df2 <- left_join(NZMS.samp.sum2, means.list.NZMS, by=c('date'="temps$dts"), copy = T)
cor.lm2 <- lm(cor.df2$mean.abund.x ~ cor.df2$adjust)
rho2 <- cor.test((cor.df2$adjust+1), (cor.df2$mean.abund.x+1), method = "spearman")

cor.df3 <- left_join(NZMS.samp.sum3, means.list.NZMS, by=c('date'="temps$dts"), copy = T)
cor.lm3 <- lm(cor.df3$mean.abund.x ~ cor.df3$adjust)
summary(cor.lm3)
rho3 <- cor.test((cor.df3$adjust+1), (cor.df3$mean.abund.x+1), method = "spearman")

rho <- mean(c(rho1$estimate, rho2$estimate, rho3$estimate))
sd <- sd(c(rho1$estimate, rho2$estimate, rho3$estimate))



abund.trends.NZMS <- ggplot(data = means.list.NZMS, aes(x = `temps$dts`,
                                                        y = (mean.abund), group = 1, color = "Model")) +
  geom_ribbon(aes(ymin = mean.abund - 1.96 * se.abund,
                   ymax = mean.abund + 1.96 * se.abund),
               colour = 'transparent',
               alpha = .15,
               show.legend = T) +
  geom_line(show.legend = T, linewidth = 0.7) +
  geom_line(data = NZMS.samp.sum[-125,], aes(x =V1, y = (log(V2)+1)*600, color = "Empirical"), show.legend = T)+
  geom_point(data = NZMS.samp.sum[125,], aes(x = V1, y = (log(V2)+1)*600, color = "Empirical"), show.legend = T)+
  #geom_line(data = flow.magnitude, aes(x = as.Date(dts), y = X_00060_00003), color = "blue") +
  #geom_line(data = temps, aes(x = as.Date(dts), y = Temperature*1000), color = "green")+
  #coord_cartesian(ylim = c(0,2000000)) +
  ylab('New Zealand Mudsnail Abundance Density (m2)') +
  xlab("")+
  labs(colour=" ")+
  theme(text = element_text(size = 14), axis.text.x = element_text(angle=45, hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), )+
  scale_x_date(date_labels="%B", date_breaks  ="6 months")




# ggplot(data = adults, aes(x = Time,
#                                    y = Adult, group = 1)) +
#   # #geom_ribbon(aes(ymin = mean.abund - 1.96 * se.abund,
#   #                 ymax = mean.abund + 1.96 * se.abund),
#   #             colour = 'transparent',
#   #             alpha = .5,
#   #             show.legend = FALSE) +
#   geom_line(show.legend = FALSE, linewidth = 0.7) +
#   geom_point()+
#   geom_line(data = NZMS.samp, aes(x = V1, y = exp(x)), color = "red")+
#   geom_line(data = flow.magnitude, aes(x = as.Date(dts), y = X_00060_00003), color = "blue") +
#   geom_line(data = temps, aes(x = as.Date(dts), y = Temperature*1000), color = "green")+
#   coord_cartesian(ylim = c(0,70000)) +
#   ylab('New Zealand Mudsnail Abundance') +
#   xlab("")+
#   theme(text = element_text(size = 14), axis.text.x = element_text(angle=45, hjust = 1, size = 12.5), 
#         axis.text.y = element_text(size = 13))+
#   scale_x_date(date_labels="%B", date_breaks  ="6 months")



#NZMS.samp.sum <- NZMS.samp.sum[which(NZMS.samp.sum$V1 %in% sample(NZMS.samp.sum$V1, size = 30)),]
#forecast::auto.arima(NZMS.samp.sum$V2, ic = "bic")
#NZMS.samp.sum$V2[2:125] <- diff(NZMS.samp.sum$V2, differences=1)




install.packages("ubms")
library(ubms)

NZMS.samp <- read.csv("NZMSsamp.csv", header = T)

# need to fit our data to an unmarked frame

# observations need to be in MxJ matrix, where M is # of sites and J is max number of obs per site

# we measured over temps$dts (those are our timesteps)

max_visits <- vector()
for (i in 1:length(temps$dts)){
  d <- NZMS.samp[which(NZMS.samp$Date >= temps$dts[i] & NZMS.samp$Date < temps$dts[i+1]), ]
  max_visits[i] <- length(d$Date)
  }

R <- length(temps$dts)
J <- max(max_visits)

site_mat <- matrix(data = NA, nrow = R, ncol = J)
dens_mat <- matrix(data = NA, nrow = R, ncol = J)
flows <- vector()
volumes <- matrix(data = NA, nrow = R, ncol = J)
for (i in 1:length(temps$dts)){
  d <- NZMS.samp[which(NZMS.samp$Date >= temps$dts[i] & NZMS.samp$Date < temps$dts[i+1]), ]
  flows[i] <- mean(d$X_00060_00003)
  site_mat[i, ] <- c(d$CountTotal, rep(NA, times = (J- length(d$CountTotal))))
  dens_mat[i, ] <- c(as.integer(d$Density), rep (NA, times = (J - length(d$Density))))
  volumes[i, ] <- c(scale(d$Volume),rep(NA, times = (J- length(d$CountTotal))))
}

# we need to remove all timesteps that are just NAs
# first identify all the timsteps that don't have data (so we can match them up later)
nodata <- which(is.na(site_mat[,1]))
length(temps$dts) - length(nodata)
site_mat <- as.matrix(site_mat[-nodata,])
dens_mat <- as.matrix(dens_mat[-nodata, ])
flows <- as.data.frame(flows[-nodata])
volumes <- as.matrix(volumes[-nodata, ])
dimnames(volumes) <- list(temps$dts[-nodata], seq(1:48))


volumes <- list(volumes)
names(volumes) <- c("vol")

dat <- unmarkedFramePCount(y = site_mat, siteCovs = flows, obsCovs = volumes)
dat_dens <- unmarkedFramePCount(y = dens_mat, siteCovs = flows, obsCovs = volumes)
# Model with no covariates or random effect on abundance
mod_null <- stan_pcount(~ 1 ~1, dat, K=11000,
                        chains=3, iter=2000, seed=123, prior_coef_det = logistic(0.69, 0.97), log_lik = F)
mod_vol <- stan_pcount(~1 + offset(vol)~1, dat, K = 11000, chains = 3, iter = 2000, seed= 123,prior_coef_det = logistic(0, 1), log_lik= F)

fmList <- fitList(Null=mod_null, vol = mod_vol)
# Model selection
modSel(fmList, nullmod="Null")

# Function returning three fit-statistics.
fitstats <- function(fm, na.rm=TRUE) {
  observed <- getY(fm@data)
  expected <- fitted(fm)
  resids <- residuals(fm)
  sse <- sum(resids^2, na.rm=na.rm)
  chisq <- sum((observed - expected)^2 / expected, na.rm=na.rm)
  freeTuke <- sum((sqrt(observed) - sqrt(expected))^2, na.rm=na.rm)
  out <- c(SSE=sse, Chisq=chisq, freemanTukey=freeTuke)
  return(out)
}
(pb <- parboot(mod_um_dens, fitstats, nsim=25, report=1))
plot(pb, main="")
# Finite-sample inference for a derived parameter.
# Population size in sampled area
Nhat <- function(fm) {
  sum(bup(ranef(fm, K=50)))
}
set.seed(345)
(pb.N <- parboot(fm, Nhat, nsim=25, report=5))
# Compare to empirical Bayes confidence intervals
colSums(confint(ranef(fm, K=50)))

# Estimates of conditional abundance distribution at each site
(re <- ranef(fm))
# Best Unbiased Predictors
bup(re, stat="mean") # Posterior mean
bup(re, stat="mode") # Posterior mode
confint(re, level=0.9) # 90% CI
# Plots
plot(re, subset=site %in% c(1:10), layout=c(5, 2), xlim=c(-1,20))

mod_um_dens <- pcount(~1 ~1, dat_dens, K = 1176244)
mod_um_dens_vol <- pcount(~scale(vol) ~1, dat_dens, K = 1176244)


dens_mod_mull <- stan_pcount(~1 ~1, dat_dens, K = 1176244, chains = 3, iter = 2000, seed = 123, log_lik = T)

# occupancy to get global p 
# site_mat[site_mat >0] <- 1
# occu_dat <- unmarkedFramePCount(y = site_mat, siteCovs = flows, obsCovs = volumes)
# 
# occu_mod_null <- stan_occu(~1 ~1, occu_dat, chains = 3, iter = 2000)
# occu_mod_vol <- stan_occu(~scale(vol)~1, occu_dat, chains = 3, iter = 2000)

#posterior
names(mod_um)
## [1] "beta_state[(Intercept)]" "beta_det[(Intercept)]"
occ_intercept <- extract(mod_um, "det")[[1]]
hist(occ_intercept, freq=FALSE)
lines(density(occ_intercept), col='red', lwd=2)

#Compare the models
#First we combine the models into a fitList:
mods <- fitList(mod_um, mod_um_vol)
#Then we generate a model selection table:
round(modSel(mods), 3)
#the model with the largest elpd performed best
plot_residuals(mod_um_vol, submodel="det")
fit_top_gof <- gof(mod_um_vol, raws=100, quiet=TRUE)

prob_det <- predict(mod_um_vol, submodel="det")
# from predicted detection probabilities, get range
range(prob_det$Predicted, na.rm = T)


