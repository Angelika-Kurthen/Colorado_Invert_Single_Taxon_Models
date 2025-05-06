##################
# GAMM N Mixture
#################
# data retrieval tool from USGS
#install.packages("dataRetrieval")
library(dataRetrieval)
#install.packages("devtools")
#install_github(repo = "jmuehlbauer-usgs/R-packages", subdir = "foodbase")
library(devtools)
library(foodbase)
library(lubridate)


library(purrr)
library(tidyverse)
library(lubridate)
library(plyr)
library(dplyr)
library(ggplot2)
# data retrieval tool from USGS
library(dataRetrieval)
library(foodbase)

source("1spFunctions.R")
source("GAMM_1sp_model.R")

temp <- readNWISdv("09380000", "00010", "2007-10-01", "2023-05-01")
temps <- TimestepTemperature(temp)
discharge <- readNWISdv("09380000", "00060", "2007-10-01", "2023-05-01")
flow.magnitude <- TimestepDischarge(discharge, 85000)

drift.data.total <- readDB(gear = "Drift", type = "Sample", updater = F)
# get drift data from between Lees Ferry and RM -6
drift.LF <- drift.data.total[which(drift.data.total$RiverMile >= -6 & drift.data.total$RiverMile <= 0),]
#specify
GAMM.samp.LF <- sampspec(samp = drift.LF, species = c("GAMM"), stats = T)
# pull stats and merge with sample info
GAMM.samp <- merge(GAMM.samp.LF$Statistics, GAMM.samp.LF$Samples, by = "BarcodeID", all = T)
# make sure we are using the same gear
GAMM.samp <- GAMM.samp[which(GAMM.samp$GearID == 4),] 
GAMM.samp <- GAMM.samp[which(GAMM.samp$FlagStrange == 0), ]
# calculate density
GAMM.samp$Density <- GAMM.samp$CountTotal/GAMM.samp$Volume
GAMM.samp <- merge(GAMM.samp, discharge[, 3:4], by = "Date")
GAMM.samp <- aggregate(GAMM.samp$Density, list(GAMM.samp$Date), FUN = mean)

means <- vector()
for (i in 1:length(temps$dts)){
  d <- GAMM.samp[which(GAMM.samp$Group.1 >= temps$dts[i] & GAMM.samp$Group.1 < temps$dts[i+1]),]
  if (any(is.nan(mean(d$x))) == T || any(is.na((d$x) == T))) {
    s = NA
  } else {
    s <- mean(d$x)}
  means <- append(means, s)
  # we know the last value doesn't fit into interval and needs to be manually added
}
means <- as.data.frame(cbind(means, as.Date(temps$dts)))
means$V2 <- as.Date(means$V2, origin = "1970-01-01")
#acf(na.omit(means$means))

set.seed(111)
out <- GAMMmodel(flow.data = flow.magnitude$Discharge, temp.data = temps, disturbanceK = 40000, baselineK = 10000, Qmin = 0.25, extinct = 50, iteration = 1000, peaklist = 0.17, peakeach = length(temps$Temperature))

repdf <- plyr::adply(out, c(1,2,3))
names(repdf) <- c('timesteps', 'stage', 'rep', 'abund')
repdf$timesteps <- as.numeric(as.character(repdf$timesteps))
repdf$timesteps <- as.factor(repdf$timesteps)
means.list.GAMM <- repdf %>%
  dplyr::group_by(timesteps, rep) %>% # combining stages
  dplyr::summarise(abund = sum(abund)) %>%
  ungroup() %>%
  dplyr::group_by(timesteps) %>%
  dplyr::summarise(mean.abund = mean(abund),
                   sd.abund = sd(abund),
                   se.abund = sd(abund)/sqrt(1000)) %>%
  ungroup()

means.list.GAMM <- as.data.frame(cbind(means.list.GAMM[1:404,], temps$dts))
colnames(means.list.GAMM) <- c("timesteps", "mean.abund", "sd.abund", "se.abund", "Date")
means.list.GAMM$Date <- as.Date(as.POSIXct(means.list.GAMM$Date, origin = "1970-01-01"))
means.list.GAMM <- as.data.frame(cbind(means.list.GAMM, means))
means.list.GAMM <- means.list.GAMM[-(1:199), ] # use first 150 as burn in 
#means.list.GAMM <- means.list.GAMM[which(means.list.GAMM$Date < "2019-02-06"),] # we have some one off data points
cor.df <- na.omit(means.list.GAMM)

# checking for temporal autocorrelation
#acf(cor.df$means) # there is some

#splitting data by lags
cor.df1 <- cor.df %>% slice(which(row_number() %% 2 == 0))
rho1 <- cor.test(cor.df1$mean.abund, cor.df1$means, method = "spearman")
cor.df2 <- cor.df %>%  slice(which(row_number() %% 2 == 1))
rho2 <- cor.test(cor.df2$mean.abund, cor.df2$means, method = "spearman")

rho <- mean(c(rho1$estimate, rho2$estimate))
#cortests <- append(cortests, rho)
#}
#}
  #cor.test((cor.df$means), (cor.df$mean.abund), method = "spearman")

means.list.GAMMnonas <- na.omit(means.list.GAMM)

rmse.gamm<- sqrt(mean((means.list.GAMMnonas$means - means.list.GAMMnonas$mean.abund)^2))

rmse.gamm.scale <- sqrt(mean((scale(means.list.GAMMnonas$means) - scale(means.list.GAMMnonas$mean.abund))^2))
coverage <- mean(scale(means.list.GAMMnonas$means) >= (scale(means.list.GAMMnonas$mean.abund) - (1.96*rmse.gamm.scale)) & scale(means.list.GAMMnonas$means) <= (scale(means.list.GAMMnonas$mean.abund) + (1.96*rmse.gamm.scale)))
colors <- c("#CCBB44", "black" )


GAMMts <- ggplot(data = means.list.GAMMnonas[-nrow(means.list.GAMMnonas), ], aes(x = Date)) +
  geom_ribbon(aes(ymin = scale(mean.abund) - 1.96 * rmse.gamm.scale,
                  ymax = scale(mean.abund) + 1.96 * rmse.gamm.scale),
              colour = 'transparent',
              alpha = .1,
              fill = "black",
              show.legend = FALSE) +
  
  geom_line(aes(y = scale(mean.abund), color = "Model"), linewidth = 1, alpha = 0.8) +
  geom_line(aes(y = scale(means), color = "Empirical"), linewidth = 1, alpha = 0.8) +
  #geom_point(aes(x = last(means.list.GAMMnonas$Date), y = last(scale(means.list.GAMMnonas$mean.abund)), color = "Model"), size = 2, show.legend = FALSE) +
  #geom_point(aes(x = last(means.list.GAMMnonas$Date), y = last(scale(means.list.GAMMnonas$means)), color = "Empirical"), size = 2, show.legend = FALSE) +
  geom_text(mapping = aes(x = as.Date("2017-12-01"), y = 5, label = paste('rho', "==", 0.41)), parse = TRUE, color = "black", size = 5.5) +
  geom_text(mapping = aes(x = as.Date("2017-12-01"), y = 5.75, label = 'C = 98%'), color = "black", size = 5.5) +
  geom_text(mapping = aes(x = as.Date("2017-12-01"), y = 6.5, label = 'Scaled RMSE = 1.21'), color = "black", size = 5.5) +
  
  labs(y = "Scaled Abundance", title = expression(italic("G. lacustris"))) +
  xlab("") +
  labs(colour = " ") +
  theme_bw() +
  ylim(c(-4, 7)) +
  scale_color_manual(values = colors) +
  theme(text = element_text(size = 15), 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
        axis.text.y = element_text(size = 15)) +
  scale_x_date(date_labels = "%Y")

##################
# N mix models
##################
library(rjags)
library(jagsUI)
library(MCMCvis)

drift.data.total <- readDB(gear = "Drift", type = "Sample", updater = F)
# get drift data from between Lees Ferry and RM -6
drift.LF <- drift.data.total[which(drift.data.total$RiverMile >= -6 & drift.data.total$RiverMile <= 0),]
#specify
GAMM.samp.LF <- sampspec(samp = drift.LF, species = c("GAMM"), stats = T)
# pull stats and merge with sample info
GAMM.samp <- merge(GAMM.samp.LF$Statistics, GAMM.samp.LF$Samples, by = "BarcodeID", all = T)
# make sure we are using the same gear
GAMM.samp <- GAMM.samp[which(GAMM.samp$GearID == 4),] 
GAMM.samp <- GAMM.samp[which(GAMM.samp$FlagStrange == 0), ]

vals <- vector()
for (i in 1:length(temps$dts)){
  d <- GAMM.samp[which(GAMM.samp$Date %within% interval(temps$dts[i], temps$dts[i+1]-1) == T),]
  if (length(d$CountTotal) > 0) {
    s<- rep(i, times = length(d$CountTotal))
    vals <- append(vals, s)}
}

# add to data frame
GAMM.samp <- cbind(GAMM.samp, vals)

# now we need into include mean water temperature and discharge for each 
GAMM.samp <- cbind(GAMM.samp, temps$Temperature[GAMM.samp$vals])
GAMM.samp <- cbind(GAMM.samp, flow.magnitude$Discharge[GAMM.samp$vals])
max_visits <- vector() # vector to put max obs per site
means <- vector()
for (i in 1:length(temps$dts)){  # cycle through all the 14 day timesteps that we have model output for
  # pull abundances between each each timestep
  d <- GAMM.samp[which(GAMM.samp$Date >= temps$dts[i] & GAMM.samp$Date < temps$dts[i+1]), ]
  max_visits[i] <- length(d$Date) # number of observations
  means[i] <- mean(d$CountTotal) # means - this is just for checking NAs
}
# make data frame with the timestep, the mean
df <- cbind.data.frame(temps$dts, means)
# remove anythig where there are NA values of Density (means that volume or count is missing)
df <- df[!is.na(df$means), ]
# phenology may also play into dynamics, so include month column as well
#month <- month(GAMM.samp$Date)
#install.packages("aspace")
library(aspace)
df$circdate <- sin(as_radians((lubridate::yday(df$`temps$dts`)/365)*360))

# define our RxJ matrix
R <- length(temps$dts)
J <- max(max_visits)
site_mat <- matrix(data = NA, nrow = R, ncol = J)
dens_mat <- matrix(data = NA, nrow = R, ncol = J)
obs_intercept <- matrix(data = 1, nrow = R, ncol = J)
# make vector for flows at each timestep
# make RxJ matrix full of densities
# make RxJ matrix full of raw counts
# make RxJ matrix full of volumes sampled for each abundance
flows <- vector()
temperature <- vector()
volumes <- matrix(data = NA, nrow = R, ncol = J)
for (i in 1:length(temps$dts)){
  d <- GAMM.samp[which(GAMM.samp$Date >= temps$dts[i] & GAMM.samp$Date < temps$dts[i+1]), ]
  flows[i] <- mean(d$`flow.magnitude$Discharge[GAMM.samp$vals]`)
  #date <- df[]
  temperature[i] <- mean(d$`temps$Temperature[GAMM.samp$vals]`)
  #windspeed[i, ] <- c(d$WindSpeed, rep(NA, times = (J- length(d$CountTotal))))
  site_mat[i, ] <- c(d$CountTotal, rep(NA, times = (J- length(d$CountTotal))))
  #habitat[i, ] <- c(d$Habitat, rep(NA, times = (J- length(d$CountTotal))))
  #time[i, ] <- c(d$TimeElapsed,rep(NA, times = (J- length(d$CountTotal))))
  #weather[i, ] <- c(d$Weather, rep(NA, times = (J- length(d$CountTotal))))
  volumes[i, ] <- c((d$Volume),rep(NA, times = (J- length(d$CountTotal))))
}

# we need to remove all timesteps that are just NAs
nodata <- which(is.na(site_mat[,1]))
# first identify all the timsteps that don't have data (so we can match them up later)
site_mat <- as.matrix(site_mat[-nodata,]) # count data
#dens_mat <- as.matrix(dens_mat[-nodata, ]) # density data
obs_intercept <- as.matrix(obs_intercept[-nodata,]) # intercept for obs cov

flows <- as.data.frame(scale(flows[-nodata])) # site cov flow 
temperature <- as.data.frame(scale(temperature[-nodata])) # site cov temp
circdate <- as.data.frame(df$circdate)
#windspeed <- as.matrix((scale(windspeed[-nodata,])))
##windspeed[is.na(windspeed)] <- mean(windspeed, na.rm = TRUE) # replace NAs with mean duration time since NAs not allowed in predictors or offsets
#habitat <- as.matrix(habitat[-nodata,])
#Mode <- function(x) {
#  ux <- unique(x)
#  ux[which.max(tabulate(match(x, ux)))]
#}
#habitat[is.na(habitat)] <- Mode(na.omit(habitat)) # replace NAs with mean duration time since NAs not allowed in predictors or offsets
#weather <- as.matrix(weather[-nodata,])
#weather[is.na(weather)] <- Mode(na.omit(weather))

site_intercept <- rep(1, times = length(flows$V1)) 
site_covs<- as.matrix(cbind(site_intercept, temperature)) #flows,temperature, circdate)
obs_covs <- array(data= NA, dim = c(length(flows$V1),J,1))
obs_covs[,,1] <- obs_intercept                                  

#offset
offset <- as.matrix(scale(log(volumes[-nodata, ])))
offset[is.na(offset)] <- mean(offset, na.rm = TRUE) # replace NAs with mean duration time since NAs not allowed in predictors or offsets
# 

sink("N-mixturePoisGAMM.jags")
cat("
model{
    # Priors
    for(i in 1:nAlpha){ # nAlpha is the number of site predictor variables
      alpha[i] ~ dnorm(0, 0.1) # alphas are the site covariates
    }

    for(i in 1:nBeta){ # nBeta is number of site x observation predictors
      beta[i] ~ dnorm(0, 0.1) # betas are the observation covariates
    }

    # Likelihood
    for(r in 1:R){
     N[r] ~ dpois(lambda[r]) #start with pulling from Poisson
     log(lambda[r]) <- sum(alpha * XN[r, ]) #XN is matrix of site covariates

     for(j in 1:J){
      y[r, j] ~ dbinom(p[r, j], N[r]) #binomial for observed counts
      logit(p[r, j]) <- sum(off[r, j] + beta * Xp[r,j,]) #XP(obs covariates)

        ## Expected count at site r, sample j
        exp[r,j] <- N[r] * p[r, j]

        ## Discrepancy
        ## (note small value added to denominator to avoid potential divide by zero)
        ## This is the X2 + descrepancy
        E[r, j] <- pow((y[r, j] - exp[r, j]), 2) / (exp[r, j] + 0.5)

        ## Simulate new count from model
        y.rep[r, j] ~ dbinom(p[r, j], N[r])

        ## X2
        E.rep[r, j] <- pow((y.rep[r, j] - exp[r, j]), 2) / (exp[r, j] + 0.5)
     }
    }
     # chi-squared test statistics
    fit <- sum(E[,])
    fit.rep <- sum(E.rep[ , ])
    } # End model
", fill = TRUE)
sink()

jags_data <- list(y = (site_mat),
                  XN = site_covs,
                  Xp = (obs_covs),
                  J = dim(site_mat)[2], #visits
                  R = dim(site_mat)[1], #sites
                  off = (offset), #obs offset
                  nAlpha = dim(site_covs)[2],
                  nBeta = dim(obs_covs)[3])

nAlpha <- dim(site_covs)[2]
nBeta <- dim(obs_covs)[3]
jags_inits <- function(){
  list(
    N = apply(jags_data$y, 1, max, na.rm=TRUE),
    alpha=runif(nAlpha,-1,1),
    beta=runif(nBeta,-1,1))}

parameters <- c("alpha", "beta", "lambda", "p", "N")


nc <- 3
ni <- 10000
nb <- 2500
nt <- 1

Nmix_fit <- jags.model("N-mixturePoisGAMM.jags",data = jags_data, inits = jags_inits, n.chains = nc, n.adapt = 1000)

update(Nmix_fit, n.iter = 1000)

Nmix_fit_UI1 <- jagsUI::jags(data = jags_data, inits = jags_inits, parameters.to.save = parameters, model.file = "N-mixturePoisGAMM.jags",  n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, parallel = TRUE)
#
print(Nmix_fit_UI1)

zm = coda.samples(Nmix_fit, variable.names = c("lambda"), n.iter = ni, n.thin = nt)


lam <- MCMCpstr(zm, "lambda")
# N <- MCMCpstr(zm, "N")
# N <- as.data.frame(cbind(as.Date(temps$dts[-nodata]), unlist(N)))
# N$V1 <- as.Date(N$V1, origin = "1970-01-01")


lam <- as.data.frame(cbind(as.Date(temps$dts[-nodata]), unlist(lam)))
lam$V1 <- as.Date(lam$V1, origin = "1970-01-01")

# y.rep <- MCMCpstr(zm, "y.rep")
# exp <- MCMCpstr(zm, "exp")
#
# plot(unlist(y.rep), unlist(site_mat))
#
# plot(unlist(y.rep), unlist(site_mat))
# abline(0, 1)
#
# fit <- MCMCchains(zm, "fit")
# fit.rep <- MCMCchains(zm, "fit.rep")
# # mean(fit > fit.rep) # close to 1 so bad fit?
# plot(fit.rep ~ fit)
# abline(0, 1) # 1 to 1 line not even there


# fit_df <- data.frame(y = c(c(unlist(site_mat)), c(unlist(y.rep))),
#                      data = rep(c("Observed", "Simulated"), each = length(site_mat)))
# library(ggplot2)
# ggplot(fit_df, aes(x = y, fill = data)) + geom_histogram() + facet_grid(.~data)# not getting all the 0s and missing the really high #s

#cor.df <- left_join(N, means.list.GAMM, by=c('V1'="Date"), copy = T)
cor.df <- na.omit(left_join(lam, means.list.GAMM, by=c('V1'="Date"), copy = T))

cor.df1 <- cor.df %>% slice(which(row_number() %% 2 == 0))
rho1 <- cor.test((cor.df1$V2.x), (cor.df1$mean.abund), method = "spearman")

cor.df2 <- cor.df %>%  slice(which(row_number() %% 2 == 1))
rho2 <- cor.test(cor.df2$V2.x, cor.df2$mean.abund, method = "spearman")

Poislam  <- mean(c(rho1$estimate, rho2$estimate))
#cor.lm <- lm(cor.df$mean.abund ~ cor.df$V2)
#PoisN <- cor.test((cor.df$V2.x), (cor.df$mean.abund), method = "spearman")

#cor.lm <- lm(cor.df$mean.abund ~ cor.df$V2)
#Poislam <- cor.test((cor.df$V2.x), (cor.df$mean.abund), method = "spearman")

rm(zm)
#rm(Nmix_fit_UI1)
rm(Nmix_fit)

sink("N-mixtureZIPGAMM.jags")
cat("
model{
    # Priors
    omega ~ dbeta(1,1)

    for(i in 1:nAlpha){ # nAlpha is the number of site predictor         variables
      alpha[i] ~ dnorm(0, 0.1) # alphas are the site covariates
    }

    for(i in 1:nBeta){ # nBeta is number of site x observation predictors
      beta[i] ~ dnorm(0, 0.1) # betas are the observation covariates
    }

    # Likelihood
    for(r in 1:R){
     z[r] ~ dbern(omega) # either there or not
     N[r] ~ dpois(lambda[r] * z[r]) #start with pulling from Poisson with z variable
     log(lambda[r]) <- sum(alpha * XN[r, ]) #XN is matrix of site covariates

     for(j in 1:J){
      y[r, j] ~ dbinom(p[r, j], N[r]) #binomial for observed counts
      logit(p[r, j]) <- sum(off[r, j] + beta * Xp[r,j,]) #XP(obs covariates)

        ## Expected count at site r, sample j
        exp[r,j] <- N[r] * p[r, j]

        ## Discrepancy
        ## (note small value added to denominator to avoid potential divide by zero)
        ## This is the X2 + descrepancy
        E[r, j] <- pow((y[r, j] - exp[r, j]), 2) / (exp[r, j] + 0.5)

        ## Simulate new count from model
        y.rep[r, j] ~ dbinom(p[r, j], N[r])

        ## X2
        E.rep[r, j] <- pow((y.rep[r, j] - exp[r, j]), 2) / (exp[r, j] + 0.5)
     }
    }
     # chi-squared test statistics
    fit <- sum(E[,])
    fit.rep <- sum(E.rep[,])
    } # End model
", fill = TRUE)
sink()

jags_data <- list(y = (site_mat),
                  XN = site_covs,
                  Xp = (obs_covs),
                  J = dim(site_mat)[2], #visits
                  R = dim(site_mat)[1], #sites
                  off = (offset), #obs offset
                  nAlpha = dim(site_covs)[2],
                  nBeta = dim(obs_covs)[3])

nAlpha <- dim(site_covs)[2]
nBeta <- dim(obs_covs)[3]
jags_inits <- function(){
  list(
    omega = runif(1, 0.5, 0.7),
    N = apply(jags_data$y, 1, max, na.rm=TRUE),
    alpha=runif(nAlpha,-1,1),
    beta=runif(nBeta,-1,1))}

parameters <- c("omega", "alpha", "beta", "lambda", "p", "N")

nc <- 3
ni <- 10000
nb <- 2500
nt <- 1

Nmix_fit <- jags.model("N-mixtureZIPGAMM.jags",data = jags_data, inits = jags_inits, n.chains = nc, n.adapt = 1000)

update(Nmix_fit, n.iter = 1000)

Nmix_fit_UI2 <- jagsUI::jags(data = jags_data, inits = jags_inits, parameters.to.save = parameters, model.file = "N-mixtureZIPGAMM.jags",  n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, parallel = TRUE)
#
print(Nmix_fit_UI2)

zm = coda.samples(Nmix_fit, variable.names = c("lambda"), n.iter = ni, n.thin = nt)


lam <- MCMCpstr(zm, "lambda")
# N <- MCMCpstr(zm, "N")
# N <- as.data.frame(cbind(as.Date(temps$dts[-nodata]), unlist(N)))
# N$V1 <- as.Date(N$V1, origin = "1970-01-01")


lam <- as.data.frame(cbind(as.Date(temps$dts[-nodata]), unlist(lam)))
lam$V1 <- as.Date(lam$V1, origin = "1970-01-01")

# cor.df <- left_join(N, means.list.GAMM, by=c('V1'="Date"), copy = T)
# cor.lm <- lm(cor.df$mean.abund ~ cor.df$V2.x)
# ZipN <-  cor.test((cor.df$V2.x), (cor.df$mean.abund), method = "spearman")

cor.df <- na.omit(left_join(lam, means.list.GAMM, by=c('V1'="Date"), copy = T))

cor.df1 <- cor.df %>% slice(which(row_number() %% 2 == 0))
rho1 <- cor.test((cor.df1$V2.x), (cor.df1$mean.abund), method = "spearman")

cor.df2 <- cor.df %>%  slice(which(row_number() %% 2 == 1))
rho2 <- cor.test(cor.df2$V2.x, cor.df2$mean.abund, method = "spearman")

Ziplam  <- mean(c(rho1$estimate, rho2$estimate))

rm(zm)
rm(Nmix_fit)
rm(Nmix_fit_UI2)


sink("N-mixtureZIPoverdispGAMM.jags")
cat("
model{
    # Priors
    omega ~ dbeta(1,1)

    tau.p <- pow(sd.p, -2)
    sd.p ~ dunif(0,3)

    for(i in 1:nAlpha){ # nAlpha is the number of site predictor         variables
      alpha[i] ~ dnorm(0, 0.1) # alphas are the site covariates
    }

    for(i in 1:nBeta){ # nBeta is number of site x observation predictors
      beta[i] ~ dnorm(0, 0.1) # betas are the observation covariates
    }

    # Likelihood
    for(r in 1:R){
     z[r] ~ dbern(omega) # either there or not
     N[r] ~ dpois(lambda[r] * z[r]) #start with pulling from Poisson with z variable
     log(lambda[r]) <- sum(alpha * XN[r, ]) #XN is matrix of site covariates

     for(j in 1:J){
      y[r, j] ~ dbinom(p[r, j], N[r]) #binomial for observed counts
      logit(p[r, j]) <- lp[r,j]
      mu.lp[r, j] <- sum(off[r, j] + beta * Xp[r,j,]) #XP(obs covariates)
      lp[r,j] ~ dnorm(mu.lp[r,j], tau.p) #sample effect based on mean p

        ## Expected count at site r, sample j
        exp[r,j] <- N[r] * p[r, j]

        ## Discrepancy
        ## (note small value added to denominator to avoid potential divide by zero)
        ## This is the X2 + descrepancy
        E[r, j] <- pow((y[r, j] - exp[r, j]), 2) / (exp[r, j] + 0.5)

        ## Simulate new count from model
        y.rep[r, j] ~ dbinom(p[r, j], N[r])

        ## X2
        E.rep[r, j] <- pow((y.rep[r, j] - exp[r, j]), 2) / (exp[r, j] + 0.5)
     }
    }
     # chi-squared test statistics
    fit <- sum(E[,])
    fit.rep <- sum(E.rep[,])
    } # End model
", fill = TRUE)
sink()

jags_data <- list(y = (site_mat),
                  XN = site_covs,
                  Xp = (obs_covs),
                  J = dim(site_mat)[2], #visits
                  R = dim(site_mat)[1], #sites
                  off = (offset), #obs offset
                  nAlpha = dim(site_covs)[2],
                  nBeta = dim(obs_covs)[3])

nAlpha <- dim(site_covs)[2]
nBeta <- dim(obs_covs)[3]
jags_inits <- function(){
  list(
    omega = runif(1, 0.5, 0.7),
    sd.p = runif(1, 0.3, 0.7),
    N = apply(jags_data$y, 1, max, na.rm=TRUE),
    alpha=runif(nAlpha,-1,1),
    beta=runif(nBeta,-1,1))}

parameters <- c("omega", "alpha", "beta", "lambda", "p", "N")

nc <- 3
ni <- 10000
nb <- 2500
nt <- 1

Nmix_fit <- jags.model("N-mixtureZIPoverdispGAMM.jags",data = jags_data, inits = jags_inits, n.chains = nc, n.adapt = 1000)

update(Nmix_fit, n.iter = 1000)

Nmix_fit_UI3 <- jagsUI::jags(data = jags_data, inits = jags_inits, parameters.to.save = parameters, model.file = "N-mixtureZIPoverdispGAMM.jags",  n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, parallel = TRUE)

print(Nmix_fit_UI3)

zm = coda.samples(Nmix_fit, variable.names = c("lambda"), n.iter = ni, n.thin = nt)


lam <- MCMCpstr(zm, "lambda")
# p <- MCMCpstr(zm, "p")
# N <- MCMCpstr(zm, "N")
# N <- as.data.frame(cbind(as.Date(temps$dts[-nodata]), unlist(N)))
# N$V1 <- as.Date(N$V1, origin = "1970-01-01")
# 

lam <- as.data.frame(cbind(as.Date(temps$dts[-nodata]), unlist(lam)))
lam$V1 <- as.Date(lam$V1, origin = "1970-01-01")

#plot(unlist(y.rep), unlist(site_mat))
#abline(0, 1)

#
# mean(fit > fit.rep) # close to 1 so bad fit?
# plot(fit.rep ~ fit)
# abline(0, 1) # 1 to 1 line not even there


# fit_df <- data.frame(y = c(c(unlist(site_mat)), c(unlist(y.rep))),
#                     data = rep(c("Observed", "Simulated"), each = length(site_mat)))
# library(ggplot2)
# ggplot(fit_df, aes(x = y, fill = data)) + geom_histogram() + facet_grid(.~data) #still not getting all the 0s and missing the really high #s

cor.df <- na.omit(left_join(lam, means.list.GAMM, by=c('V1'="Date"), copy = T))

cor.df1 <- cor.df %>% slice(which(row_number() %% 2 == 0))
rho1 <- cor.test((cor.df1$V2.x), (cor.df1$mean.abund), method = "spearman")

cor.df2 <- cor.df %>%  slice(which(row_number() %% 2 == 1))
rho2 <- cor.test(cor.df2$V2.x, cor.df2$mean.abund, method = "spearman")

Zip_ovd_lam  <- mean(c(rho1$estimate, rho2$estimate))
# cor.df <- left_join(N, means.list.GAMM, by=c('V1'="V2"), copy = T)
# #cor.lm <- lm(cor.df$mean.abund ~ cor.df$V2.x)
# Zip_ovdN <- cor.test((cor.df$V2), (cor.df$mean.abund), method = "spearman")
# 
# cor.df <- left_join(lam, means.list.GAMM, by=c('V1'="V2"), copy = T)
# #cor.lm <- lm(cor.df$mean.abund ~ cor.df$V2.x)
# Zip_ovdlam <- cor.test((cor.df$V2), (cor.df$mean.abund), method = "spearman")

rm(zm)
rm(Nmix_fit_UI3)
rm(Nmix_fit)

sink("N-mixturePoisoverdispGAMM.jags")
cat("
model{
    # Priors
    for(i in 1:nAlpha){ # nAlpha is the number of site predictor variables
      alpha[i] ~ dnorm(0, 0.1) # alphas are the site covariates
    }

    for(i in 1:nBeta){ # nBeta is number of site x observation predictors
      beta[i] ~ dnorm(0, 0.1) # betas are the observation covariates
    }

    tau.p <- pow(sd.p, -2)
    sd.p ~ dunif(0,3)

    # Likelihood
    for(r in 1:R){
     N[r] ~ dpois(lambda[r]) #start with pulling from Poisson
     log(lambda[r]) <- sum(alpha * XN[r, ]) #XN is matrix of site covariates

     for(j in 1:J){
      y[r, j] ~ dbinom(p[r, j], N[r]) #binomial for observed counts
      logit(p[r, j]) <- lp[r,j]
      mu.lp[r, j] <- sum(off[r, j] + beta * Xp[r,j,]) #XP(obs covariates)
      lp[r,j] ~ dnorm(mu.lp[r,j], tau.p) #sample effect based on mean p

        ## Expected count at site r, sample j
        exp[r,j] <- N[r] * p[r, j]

        ## Discrepancy
        ## (note small value added to denominator to avoid potential divide by zero)
        ## This is the X2 + descrepancy
        E[r, j] <- pow((y[r, j] - exp[r, j]), 2) / (exp[r, j] + 0.5)

        ## Simulate new count from model
        y.rep[r, j] ~ dbinom(p[r, j], N[r])

        ## X2
        E.rep[r, j] <- pow((y.rep[r, j] - exp[r, j]), 2) / (exp[r, j] + 0.5)
     }
    }
     # chi-squared test statistics
    fit <- sum(E[,])
    fit.rep <- sum(E.rep[,])
    } # End model
", fill = TRUE)
sink()

jags_data <- list(y = (site_mat),
                  XN = site_covs,
                  Xp = (obs_covs),
                  J = dim(site_mat)[2], #visits
                  R = dim(site_mat)[1], #sites
                  off = (offset), #obs offset
                  nAlpha = dim(site_covs)[2],
                  nBeta = dim(obs_covs)[3])

nAlpha <- dim(site_covs)[2]
nBeta <- dim(obs_covs)[3]
jags_inits <- function(){
  list(
    sd.p = runif(1, 0.5, 0.8),
    N = apply(jags_data$y, 1, max, na.rm=TRUE),
    alpha=runif(nAlpha,-1,1),
    beta=runif(nBeta,-1,1))}

parameters <- c("omega", "alpha", "beta", "lambda", "p", "N")

nc <- 3
ni <- 10000
nb <- 2500
nt <- 1

Nmix_fit <- jags.model("N-mixturePoisoverdispGAMM.jags",data = jags_data, inits = jags_inits, n.chains = nc, n.adapt = 1000)

update(Nmix_fit, n.iter = 1000)

Nmix_fit_UI4 <- jagsUI::jags(data = jags_data, inits = jags_inits, parameters.to.save = parameters, model.file = "N-mixturePoisoverdispGAMM.jags",  n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, parallel = TRUE)

print(Nmix_fit_UI4)
rm(Nmix_fit_UI4)
zm = coda.samples(Nmix_fit, variable.names = c("lambda"), n.iter = ni, n.thin = nt)

lam <- MCMCpstr(zm, "lambda")
# p <- MCMCpstr(zm, "p")
# N <- MCMCpstr(zm, "N")
# N <- as.data.frame(cbind(as.Date(temps$dts[-nodata]), unlist(N)))
# N$V1 <- as.Date(N$V1, origin = "1970-01-01")
# 

lam <- as.data.frame(cbind(as.Date(temps$dts[-nodata]), unlist(lam)))
lam$V1 <- as.Date(lam$V1, origin = "1970-01-01")

# close to 1 so bad fit?
# plot(fit.rep ~ fit)
# abline(0, 1) # 1 to 1 line not even there

#cor.df <- left_join(N, means.list.GAMM, by=c('V1'="V2"), copy = T)
#cor.lm <- lm(cor.df$mean.abund ~ cor.df$V2.x)
#Pois_ovdN <- cor.test((cor.df$V2), (cor.df$mean.abund), method = "spearman")


#cor.df <- left_join(lam, means.list.GAMM, by=c('V1'="V2"), copy = T)
#cor.lm <- lm(cor.df$mean.abund ~ cor.df$V2.x)
#Pois_ovdlam <- cor.test((cor.df$V2), (cor.df$mean.abund), method = "spearman")


cor.df <- na.omit(left_join(lam, means.list.GAMM, by=c('V1'="Date"), copy = T))

cor.df1 <- cor.df %>% slice(which(row_number() %% 2 == 0))
rho1 <- cor.test((cor.df1$V2.x), (cor.df1$mean.abund), method = "spearman")

cor.df2 <- cor.df %>%  slice(which(row_number() %% 2 == 1))
rho2 <- cor.test(cor.df2$V2.x, cor.df2$mean.abund, method = "spearman")

Pois_ovd_lam  <- mean(c(rho1$estimate, rho2$estimate))


rm(zm)
rm(Nmix_fit)
#rm(Nmix_fit_UI4)

sink("N-mixtureNBGAMM.jags")
cat("
model{

# State model
for (r in 1:R){
 N[r] <- n[r] * z[r]
  n[r] ~ dnegbin(s[r], th)
  s[r] <- th / (th + lambda[r])
  log(lambda[r]) <- sum(alpha * XN[r, ]) #XN is matrix of site covariates
   z[r] ~ dbern(omega)
}

omega ~ dbeta(1,1)
th ~ dgamma(0.01, 0.01)
phi <- 1/th
theta <- th

  for(i in 1:nAlpha){ # nAlpha is the number of site predictor variables
      alpha[i] ~ dnorm(0, 0.01) # alphas are the site covariates
    }


# Detection model
for (r in 1:R){
  for (j in 1:J){
    logit(p[r,j]) <- max(1e-5, min (0.999999, sum(off[r, j] + (beta * Xp[r,j,]))))
    y[r,j] ~ dbinom(p[r,j], N[r])
  }
}

 for(i in 1:nBeta){ # nBeta is number of site x observation predictors
      beta[i] ~ dnorm(0, 0.01) # betas are the observation covariates
    }


# Fit statistic for real data
for (r in 1:R){
  for (j in 1:J){
    yhat[r,j] <- N[r] * p[r,j] + 0.001 # add small value to avoid divide by zero
    chi2[r,j] <- (y[r,j] - yhat[r,j])^2 / yhat[r,j]
  }
}
fit <- sum(chi2)

# Fit statistic for simulated data
for (r in 1:R){
  for (j in 1:J){
    y_new[r,j] ~ dbinom(p[r,j], N[r]) # simulate new datapoint
    chi2_new[r,j] <- (y_new[r,j] - yhat[r,j])^2 / yhat[r,j]
  }
}
fit_new <- sum(chi2_new)

sumN <- sum(N[])

}
", fill = TRUE)
sink()

jags_data <- list(y = (site_mat),
                  XN = site_covs,
                  Xp = (obs_covs),
                  J = dim(site_mat)[2], #visits
                  R = dim(site_mat)[1], #sites
                  off = (offset), #obs offset
                  nAlpha = dim(site_covs)[2],
                  nBeta = dim(obs_covs)[3])



nAlpha <- dim(site_covs)[2]
nBeta <- dim(obs_covs)[3]
jags_inits <- function(){
  list(
    omega = runif(1, 0.5, 0.7),
    n = apply(jags_data$y, 1, max, na.rm=TRUE),
    alpha=runif(nAlpha,-1,1),
    beta=runif(nBeta,-1,1))}

parameters <- c("alpha", "beta", "lambda", "p", "N", "theta", "phi", "fit", "fit_new")

nc <- 3
ni <- 10000
nb <- 2500
nt <- 1

Nmix_fit <- jags.model("N-mixtureNBGAMM.jags",data = jags_data, inits = jags_inits, n.chains = nc, n.adapt = 1000)

update(Nmix_fit, n.iter = 1000)



zm = coda.samples(Nmix_fit, variable.names = c("lambda", "N"), n.iter = ni, n.thin = nt)


lam <- MCMCpstr(zm, "lambda")
# N <- MCMCpstr(zm, "N")
# N <- as.data.frame(cbind(as.Date(temps$dts[-nodata]), unlist(N)))
# N$V1 <- as.Date(N$V1, origin = "1970-01-01")
# 

lam <- as.data.frame(cbind(as.Date(temps$dts[-nodata]), unlist(lam)))
lam$V1 <- as.Date(lam$V1, origin = "1970-01-01")


Nmix_fit_UI5 <- jagsUI::jags(data = jags_data, inits = jags_inits, parameters.to.save = parameters, model.file = "N-mixtureNBGAMM.jags",  n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, parallel = TRUE)

# cor.df <- left_join(lam, means.list.GAMM, by=c('V1'="V2"), copy = T)
# #cor.lm <- lm(cor.df$mean.abund ~ cor.df$V2.x)
# nblam <- cor.test((cor.df$V2), (cor.df$mean.abund), method = "spearman")
# 
# 
# cor.df <- left_join(N, means.list.GAMM, by=c('V1'="V2"), copy = T)
# #cor.lm <- lm(cor.df$mean.abund ~ cor.df$V2.x)
# nbN <- cor.test((cor.df$V2), (cor.df$mean.abund), method = "spearman")

cor.df <- na.omit(left_join(lam, means.list.GAMM, by=c('V1'="Date"), copy = T))

cor.df1 <- cor.df %>% slice(which(row_number() %% 2 == 0))
rho1 <- cor.test((cor.df1$V2.x), (cor.df1$mean.abund), method = "spearman")

cor.df2 <- cor.df %>%  slice(which(row_number() %% 2 == 1))
rho2 <- cor.test(cor.df2$V2.x, cor.df2$mean.abund, method = "spearman")

NB_lam <- mean(c(rho1$estimate, rho2$estimate))

print(Nmix_fit_UI5)


rm(zm)
rm(Nmix_fit)
rm(Nmix_fit_UI5)
