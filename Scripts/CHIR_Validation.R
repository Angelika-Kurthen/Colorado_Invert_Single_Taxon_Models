##################
# CHIR N Mixture
#################
# data retrieval tool from USGS
#install.packages("dataRetrieval")
library(dataRetrieval)
#install.packages("devtools")
#install_github(repo = "jmuehlbauer-usgs/R-packages", subdir = "foodbase")
library(devtools)
library(foodbase)
library(lubridate)
library(parallel)
library(doParallel)
library(foreach)
library(purrr)
library(tidyverse)
#library(car)
library(lubridate)
library(plyr)
library(dplyr)
library(ggplot2)
# data retrieval tool from USGS
library(dataRetrieval)
library(foodbase)

getwd()

source("1spFunctions.R")
source("CHIR_1sp_Model.R")


temp <- readNWISdv("09380000", "00010", "2007-10-01", "2023-05-01")
temps <- TimestepTemperature(temp)
discharge <- readNWISdv("09380000", "00060", "2007-10-01", "2023-05-01")
flow.magnitude <- TimestepDischarge(discharge, 85000)

drift.data.total <- foodbase::readDB(gear = "Drift", type = "Sample", updater = F)
# get drift data from between Lees Ferry and RM -6
drift.LF <- drift.data.total[which(drift.data.total$RiverMile >= -6 & drift.data.total$RiverMile <= 0),]
#specify
CHIR.samp.LF <- sampspec(samp = drift.LF, species = c("CHIL"), stats = T)
# pull stats and merge with sample info
CHIR.samp <- merge(CHIR.samp.LF$Statistics, CHIR.samp.LF$Samples, by = "BarcodeID", all = T)
# make sure we are using the same gear
CHIR.samp <- CHIR.samp[which(CHIR.samp$GearID == 4),]
CHIR.samp <- CHIR.samp[which(CHIR.samp$FlagStrange == 0), ]
# calculate density
CHIR.samp$Density <- CHIR.samp$CountTotal/CHIR.samp$Volume
CHIR.samp <- merge(CHIR.samp, discharge[, 3:4], by = "Date")

CHIR.samp$Density <- (CHIR.samp$Density)/(1.3 *(CHIR.samp$X_00060_00003 * 0.0283168)^4.1)
CHIR.samp <- aggregate(CHIR.samp$Density, list(CHIR.samp$Date), FUN = mean)

# write.csv2(CHIR.samp, file = "CHIR.samp.csv")
# CHIR.samp <- read.csv2("CHIR.samp.csv")
means <- vector()
for (i in 1:length(temps$dts)){
  d <- CHIR.samp[which(CHIR.samp$Group.1 >= temps$dts[i] & CHIR.samp$Group.1 < temps$dts[i+1]),]
  if (any(is.nan(mean(d$x))) == T || any(is.na((d$x) == T))) {
    s = NA
  } else {
    s <- mean(d$x)}
  means <- append(means, s)
  # we know the last value doesn't fit into interval and needs to be manually added
}
means <- as.data.frame(cbind(means, as.Date(temps$dts)))
means$V2 <- as.Date(means$V2, origin = "1970-01-01")

  # Run CHIR model
  set.seed(111)
  out <- CHIRmodel(flow.data = flow.magnitude$Discharge, temp.data = temps, disturbanceK = 40000, 
                   baselineK = 10000, Qmin = 0.25, extinct = 50, iteration = 1000, 
                   peaklist = 0.17, peakeach = length(temps$Temperature) )#, g1s = g1s, g2s = g2s)
  
  # Process the model output
  repdf <- plyr::adply(out, c(1, 2, 3))
  names(repdf) <- c('timesteps', 'stage', 'rep', 'abund')
  repdf$timesteps <- as.numeric(as.character(repdf$timesteps))
  repdf <- subset(repdf, stage == "S1" | stage == "S2")
  means.list.CHIR <- repdf %>%
    dplyr::group_by(timesteps, rep) %>%
    dplyr::summarise(abund = sum(abund)) %>%
    ungroup() %>%
    dplyr::group_by(timesteps) %>%
    dplyr::summarise(mean.abund = mean(abund),
                     sd.abund = sd(abund),
                     se.abund = sd(abund) / sqrt(1000)) %>%
    ungroup()
  
  means.list.CHIR <- as.data.frame(cbind(means.list.CHIR[1:404, ], temps$dts))
  colnames(means.list.CHIR)[5] <- "Date"
  means.list.CHIR <- as.data.frame(cbind(means.list.CHIR, means))
  means.list.CHIR <- means.list.CHIR[-(1:199), ]  # Use first 200 as burn-in
  
  # Calculate Spearman's rho
  cor.df <- na.omit(means.list.CHIR)
  cor.df1 <- cor.df %>% slice(which(row_number() %% 2 == 0))
  rho1 <- cor.test(cor.df1$mean.abund, cor.df1$means, method = "spearman")
  cor.df2 <- cor.df %>% slice(which(row_number() %% 2 == 1))
  rho2 <- cor.test(cor.df2$mean.abund, cor.df2$means, method = "spearman")
  rho <- mean(c(rho1$estimate, rho2$estimate))


rmse.chir<- sqrt(mean((cor.df$means - cor.df$mean.abund)^2))

rmse.chir.scale <- sqrt(mean((scale(cor.df$means) - scale(cor.df$mean.abund))^2))
# coverage 
coverage <- mean(scale(cor.df$means) >= (scale(cor.df$mean.abund) - (1.96*rmse.chir.scale)) & scale(cor.df$means) <= (scale(cor.df$mean.abund) + (1.96*rmse.chir.scale)))



colors <- c("#228833", "black" )
linetypes <- c("solid", "twodash")

means.list.CHIR <- means.list.CHIR[which(means.list.CHIR$Date < "2021-07-07"),]
means.list.CHIR$Date <- as.Date(means.list.CHIR$Date)

CHIRts <- ggplot(data = subset(means.list.CHIR[-nrow(means.list.CHIR), ], !is.na(means)), aes(x = Date,  y = scale(mean.abund), group = 1, color = "Model")) +
  geom_ribbon(aes(ymin = scale(mean.abund) - 1.96 * rmse.chir.scale,
                  ymax = scale(mean.abund) + 1.96 * rmse.chir.scale),
              colour = 'transparent',
              fill = "black",
              alpha = .15,
              show.legend = F) +
  geom_line(show.legend = T, linewidth = 1, alpha = 0.8) +
  geom_line(data = subset(means.list.CHIR[-nrow(means.list.CHIR), ], !is.na(means)),
            aes(x = Date, y = scale(means), color = "Empirical"), linewidth = 1,  show.legend = T, alpha = 0.8)+
  #geom_point(aes(x = last(means.list.CHIR$Date), y = last(scale(means.list.CHIR$mean.abund)), color = "Model"), size = 2, show.legend = FALSE) +
  #geom_point(aes(x = last(means.list.CHIR$Date), y = last(scale(means.list.CHIR$means)), color = "Empirical"), size = 2, show.legend = FALSE) +
  geom_text(mapping = aes(x = as.Date("2018-01-01"), y =5, label = paste('rho', "==", 0.18)), parse = T, color = "black", size = 5.5)+
  geom_text(mapping = aes(x = as.Date("2018-01-01"), y =5.75, label = paste("C = 92%")), color = "black", size = 5.5)+
  geom_text(mapping = aes(x = as.Date("2018-01-01"), y =6.5,label = paste("Scaled RMSE = 1.21")), color = "black", size = 5.5)+
  ylab("Scaled Abundance")+
  xlab("")+
  ylim(c(-4, 7))+
  labs(colour=" ", title = expression(paste(italic("Chironomidae"), " spp.")))+
  theme_bw()+
  scale_color_manual(values = colors)+

theme(text = element_text(size = 15), axis.text.x = element_text(angle=45, hjust = 1, size = 15),
        axis.text.y = element_text(size = 15), )+
  scale_x_date(date_labels="%Y")

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
CHIR.samp.LF <- sampspec(samp = drift.LF, species = c("CHIL"), stats = T)
# pull stats and merge with sample info
CHIR.samp <- merge(CHIR.samp.LF$Statistics, CHIR.samp.LF$Samples, by = "BarcodeID", all = T)
# make sure we are using the same gear
CHIR.samp <- CHIR.samp[which(CHIR.samp$GearID == 4),]
CHIR.samp <- CHIR.samp[which(CHIR.samp$FlagStrange == 0), ]

vals <- vector()
for (i in 1:length(temps$dts)){
  d <- CHIR.samp[which(CHIR.samp$Date %within% interval(temps$dts[i], temps$dts[i+1]-1) == T),]
  if (length(d$CountTotal) > 0) {
    s<- rep(i, times = length(d$CountTotal))
    vals <- append(vals, s)}
}

# add to data frame
CHIR.samp <- cbind(CHIR.samp, vals)

# now we need into include mean water temperature and discharge for each
CHIR.samp <- cbind(CHIR.samp, temps$Temperature[CHIR.samp$vals])
CHIR.samp <- cbind(CHIR.samp, flow.magnitude$Discharge[CHIR.samp$vals])

max_visits <- vector() # vector to put max obs per site
means <- vector()
for (i in 1:length(temps$dts)){  # cycle through all the 14 day timesteps that we have model output for
  # pull abundances between each each timestep
  d <- CHIR.samp[which(CHIR.samp$Date >= temps$dts[i] & CHIR.samp$Date < temps$dts[i+1]), ]
  max_visits[i] <- length(d$Date) # number of observations
  means[i] <- mean(d$CountTotal) # means - this is just for checking NAs
}
# make data frame with the timestep, the mean
df <- cbind.data.frame(temps$dts, means)
# remove anythig where there are NA values of Density (means that volume or count is missing)
df <- df[!is.na(df$means), ]
# phenology may also play into dynamics, so include month column as well
#month <- month(CHIR.samp$Date)
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
  d <- CHIR.samp[which(CHIR.samp$Date >= temps$dts[i] & CHIR.samp$Date < temps$dts[i+1]), ]
  flows[i] <- mean(d$`flow.magnitude$Discharge[CHIR.samp$vals]`)
  #date <- df[]
  temperature[i] <- mean(d$`temps$Temperature[CHIR.samp$vals]`)
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
obs_intercept <- as.matrix(obs_intercept[-nodata,]) # intercept for obs cov

flows <- as.data.frame(scale(flows[-nodata])) # site cov flow
temperature <- as.data.frame(scale(temperature[-nodata])) # site cov temp
circdate <- as.data.frame(df$circdate)


site_intercept <- rep(1, times = length(flows$V1))
site_covs<- as.matrix(cbind(site_intercept, temperature, circdate)) #flows,temperature, circdate)
obs_covs <- array(data= NA, dim = c(length(flows$V1),J,1))
obs_covs[,,1] <- obs_intercept

#offset
offset <- as.matrix(scale(log(volumes[-nodata, ])))
offset[is.na(offset)] <- mean(offset, na.rm = TRUE) # replace NAs with mean duration time since NAs not allowed in predictors or offsets
#
sink("N-mixturePoisCHIR.jags")
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

Nmix_fit <- jags.model("N-mixturePoisCHIR.jags",data = jags_data, inits = jags_inits, n.chains = nc, n.adapt = 1000)

update(Nmix_fit, n.iter = 1000)

Nmix_fit_UI1 <- jagsUI::jags(data = jags_data, inits = jags_inits, parameters.to.save = parameters, model.file = "N-mixturePoisCHIR.jags",  n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, parallel = TRUE)
#
#print(Nmix_fit_UI1)

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

# cor.df <- left_join(N, means.list.CHIR, by=c('V1'="Date"), copy = T)
# #cor.lm <- lm(cor.df$mean.abund ~ cor.df$V2)
# PoisN <- cor.test((cor.df$V2.x), (cor.df$mean.abund), method = "spearman")

cor.df <- left_join(lam, means.list.CHIR, by=c('V1'="Date"), copy = T)

cor.df1 <- cor.df %>% slice(which(row_number() %% 2 == 0))
rho1 <- cor.test(cor.df1$V2.x, cor.df1$means, method = "spearman")
cor.df2 <- cor.df %>% slice(which(row_number() %% 2 == 1))
rho2 <- cor.test(cor.df2$V2.x, cor.df2$means, method = "spearman")
Poislam <- mean(c(rho1$estimate, rho2$estimate))


#cor.lm <- lm(cor.df$mean.abund ~ cor.df$V2)

rm(zm)
#rm(Nmix_fit_UI1)
rm(Nmix_fit)


sink("N-mixtureZIPCHIR.jags")
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

Nmix_fit <- jags.model("N-mixtureZIPCHIR.jags",data = jags_data, inits = jags_inits, n.chains = nc, n.adapt = 1000)

update(Nmix_fit, n.iter = 1000)

Nmix_fit_UI2 <- jagsUI::jags(data = jags_data, inits = jags_inits, parameters.to.save = parameters, model.file = "N-mixtureZIPCHIR.jags",  n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, parallel = TRUE)
#
print(Nmix_fit_UI2)

zm = coda.samples(Nmix_fit, variable.names = c("lambda"), n.iter = ni, n.thin = nt)


lam <- MCMCpstr(zm, "lambda")
# N <- MCMCpstr(zm, "N")
# N <- as.data.frame(cbind(as.Date(temps$dts[-nodata]), unlist(N)))
# N$V1 <- as.Date(N$V1, origin = "1970-01-01")
# 

lam <- as.data.frame(cbind(as.Date(temps$dts[-nodata]), unlist(lam)))
lam$V1 <- as.Date(lam$V1, origin = "1970-01-01")


# cor.df <- left_join(N, means.list.CHIR, by=c('V1'="Date"), copy = T)
# cor.lm <- lm(cor.df$mean.abund ~ cor.df$V2.x)
# ZipN <-  cor.test((cor.df$V2.x), (cor.df$mean.abund), method = "spearman")

cor.df <- left_join(lam, means.list.CHIR, by=c('V1'="Date"), copy = T)

cor.df1 <- cor.df %>% slice(which(row_number() %% 2 == 0))
rho1 <- cor.test(cor.df1$V2.x, cor.df1$means, method = "spearman")
cor.df2 <- cor.df %>% slice(which(row_number() %% 2 == 1))
rho2 <- cor.test(cor.df2$V2.x, cor.df2$means, method = "spearman")
Ziplam <- mean(c(rho1$estimate, rho2$estimate))

# cor.lm <- lm(cor.df$mean.abund ~ cor.df$V2.x)
# Ziplam <- cor.test((cor.df$V2.x), (cor.df$mean.abund), method = "spearman")

rm(zm)
rm(Nmix_fit)
#rm(Nmix_fit_UI2)


sink("N-mixtureZIPoverdispCHIR.jags")
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

Nmix_fit <- jags.model("N-mixtureZIPoverdispCHIR.jags",data = jags_data, inits = jags_inits, n.chains = nc, n.adapt = 1000)

update(Nmix_fit, n.iter = 1000)

Nmix_fit_UI3 <- jagsUI::jags(data = jags_data, inits = jags_inits, parameters.to.save = parameters, model.file = "N-mixtureZIPoverdispCHIR.jags",  n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, parallel = TRUE)

print(Nmix_fit_UI3)

zm = coda.samples(Nmix_fit, variable.names = c("lambda"), n.iter = ni, n.thin = nt)


lam <- MCMCpstr(zm, "lambda")
# p <- MCMCpstr(zm, "p")
# N <- MCMCpstr(zm, "N")
# N <- as.data.frame(cbind(as.Date(temps$dts[-nodata]), unlist(N)))
# N$V1 <- as.Date(N$V1, origin = "1970-01-01")


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

# cor.df <- left_join(N, means.list.CHIR, by=c('V1'="Date"), copy = T)
# cor.lm <- lm(cor.df$mean.abund ~ cor.df$V2.x)
# Zip_ovdN <- cor.test((cor.df$V2.x), (cor.df$mean.abund), method = "spearman")

cor.df <- left_join(lam, means.list.CHIR, by=c('V1'="Date"), copy = T)

cor.df1 <- cor.df %>% slice(which(row_number() %% 2 == 0))
rho1 <- cor.test(cor.df1$V2.x, cor.df1$means, method = "spearman")
cor.df2 <- cor.df %>% slice(which(row_number() %% 2 == 1))
rho2 <- cor.test(cor.df2$V2.x, cor.df2$means, method = "spearman")
Zip_ovdlam <- mean(c(rho1$estimate, rho2$estimate))


rm(zm)
#rm(Nmix_fit_UI3)
rm(Nmix_fit)

sink("N-mixturePoisoverdispCHIR.jags")
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

Nmix_fit <- jags.model("N-mixturePoisoverdispCHIR.jags",data = jags_data, inits = jags_inits, n.chains = nc, n.adapt = 1000)

update(Nmix_fit, n.iter = 1000)

Nmix_fit_UI4 <- jagsUI::jags(data = jags_data, inits = jags_inits, parameters.to.save = parameters, model.file = "N-mixturePoisoverdispCHIR.jags",  n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, parallel = TRUE)

print(Nmix_fit_UI4)

zm = coda.samples(Nmix_fit, variable.names = c("lambda"), n.iter = ni, n.thin = nt)

lam <- MCMCpstr(zm, "lambda")
# p <- MCMCpstr(zm, "p")
# N <- MCMCpstr(zm, "N")
# N <- as.data.frame(cbind(as.Date(temps$dts[-nodata]), unlist(N)))
# N$V1 <- as.Date(N$V1, origin = "1970-01-01")


lam <- as.data.frame(cbind(as.Date(temps$dts[-nodata]), unlist(lam)))
lam$V1 <- as.Date(lam$V1, origin = "1970-01-01")

# close to 1 so bad fit?
# plot(fit.rep ~ fit)
# abline(0, 1) # 1 to 1 line not even there
cor.df1 <- cor.df %>% slice(which(row_number() %% 2 == 0))
rho1 <- cor.test(cor.df1$V2.x, cor.df1$means, method = "spearman")
cor.df2 <- cor.df %>% slice(which(row_number() %% 2 == 1))
rho2 <- cor.test(cor.df2$V2.x, cor.df2$means, method = "spearman")
Pois_ovdlam <- mean(c(rho1$estimate, rho2$estimate))

rm(zm)
rm(Nmix_fit)
#rm(Nmix_fit_UI4)

sink("N-mixtureNBCHIR.jags")
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

Nmix_fit <- jags.model("N-mixtureNBCHIR.jags",data = jags_data, inits = jags_inits, n.chains = nc, n.adapt = 1000)

update(Nmix_fit, n.iter = 1000)



zm = coda.samples(Nmix_fit, variable.names = c("lambda"), n.iter = ni, n.thin = nt)


lam <- MCMCpstr(zm, "lambda")
# N <- MCMCpstr(zm, "N")
# N <- as.data.frame(cbind(as.Date(temps$dts[-nodata]), unlist(N)))
# N$V1 <- as.Date(N$V1, origin = "1970-01-01")


lam <- as.data.frame(cbind(as.Date(temps$dts[-nodata]), unlist(lam)))
lam$V1 <- as.Date(lam$V1, origin = "1970-01-01")


Nmix_fit_UI5 <- jagsUI::jags(data = jags_data, inits = jags_inits, parameters.to.save = parameters, model.file = "N-mixtureNBCHIR.jags",  n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, parallel = TRUE)

print(Nmix_fit_UI5)

cor.df <- left_join(lam, means.list.CHIR, by=c('V1'="Date"), copy = T)
cor.df1 <- cor.df %>% slice(which(row_number() %% 2 == 0))
rho1 <- cor.test(cor.df1$V2.x, cor.df1$means, method = "spearman")
cor.df2 <- cor.df %>% slice(which(row_number() %% 2 == 1))
rho2 <- cor.test(cor.df2$V2.x, cor.df2$means, method = "spearman")
nblam <- mean(c(rho1$estimate, rho2$estimate))


# cor.df <- left_join(N, means.list.CHIR, by=c('V1'="Date"), copy = T)
# cor.lm <- lm(cor.df$mean.abund ~ cor.df$V2.x)
# nbN <- cor.test((cor.df$V2.x), (cor.df$mean.abund), method = "spearman")

rm(zm)
rm(Nmix_fit)
#rm(Nmix_fit_UI5)
