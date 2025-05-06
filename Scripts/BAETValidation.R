##############################
## Code to validate BAET model
##############################
library(purrr)
library(tidyverse)
library(lubridate)
library(plyr)
library(dplyr)
library(ggplot2)
# data retrieval tool from USGS
library(dataRetrieval)
library(rjags)
library(MCMCvis)

# load Baet model
source("BAET_1sp_Model.R")
# pull discharge and temps from below flaming gorge dam
discharge <- readNWISdv("09234500", "00060", "1986-10-01", "1999-10-06")
# Bankfull discharge for Green River from https://digitalcommons.usu.edu/cgi/viewcontent.cgi?article=7604&context=etd
flow.magnitude <- TimestepDischarge(discharge, 22424.813)
temp <- readNWISdv("09234500", "00010", "2004-02-05", "2023-05-01")
temps <- average.yearly.temp(temp, "X_00010_00003","Date")
temps <- rep.avg.year(temps, 15, change.in.temp = 0, years.at.temp = 15)
# align dates
temps <- temps[20:359,2:3]
temps$dts <- flow.magnitude$dts

out <- BAETmodel(flow.data = flow.magnitude$Discharge, temp.data = temps, disturbanceK = 40000, baselineK = 5000, Qmin = 0.1, extinct = 50, iteration = 1000, peaklist = 0.13, peakeach = length(temps$Temperature), stage_output = "larvae")

# upload larval baet data from Flaming Gorge Dam 
bugdata <- read_delim("bugdata.txt", delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
bugdata <- as.data.frame(bugdata[-c(1:6, 3732:3740),])
names(bugdata) <- c("Sample", "Location", "Date", "Citation", "Method", "Area", "Density", "Phylum", "Class", "Order", "Family", "Subfamily", "Genus", "Species")
bugdata <- bugdata[which(bugdata$Location == "0.8KDD" | bugdata$Location == "6KDD" | bugdata$Location == "12KDD"),]
bugdata$Date <- as.Date(bugdata$Date, "%m/%d/%Y")
bugdata$Density <- as.numeric(bugdata$Density)
BAETdata <- bugdata[which(bugdata$Date >= "1986-10-01" & bugdata$Family == "Baetidae"), ]

BAET.samp <- aggregate(BAETdata$Density, list(BAETdata$Date), FUN = sum)

means <- vector()
for (i in 1:length(temps$dts)){
  d <- BAET.samp[which(BAET.samp$Group.1 >= temps$dts[i] & BAET.samp$Group.1 < temps$dts[i+1]),]
  if (is.nan(mean(d$x)) == T || is.na(mean(d$x)) == T) {
    s = NA
  } else {
    s<- mean(d$x)}
  means <- append(means, s)
  # we know the last value doesn't fit into interval and needs to be manually added
}
means[length(temps$dts)] <- BAET.samp$x[length(BAET.samp$x)]

means.list.BAET <- mean.data.frame(out, burnin = 200, iteration= 1000)
means.list.BAET <- cbind(means.list.BAET, temps$dts[200:341])
means.list.BAET$`temps$dts` <- as.Date(means.list.BAET$`temps$dts`)

BAET.samp.sum <- na.omit(as.data.frame(cbind(as.Date(means.list.BAET$`temps$dts`), means[200:340])))
BAET.samp.sum$V1 <- as.Date(BAET.samp.sum$V1, origin = "1970-01-01")


cor.df <- left_join(BAET.samp.sum, means.list.BAET, by=c('V1'="temps$dts"), copy = T)


rho <- cor.test((cor.df$V2), (cor.df$mean.abund), method = "spearman")

colors <- c("#66CCEE", "black")
# ggplot(data = cor.df, aes(x = (V2) , y = (mean.abund)))+
#   geom_point()+
#   stat_smooth(method = "lm",
#               formula = y ~ x,
#               geom = "smooth")+
#   geom_text(x = 1000, y = 3250, label = " ")+
#   labs(y = "Baetidae Model Output", x = "Baetidae Emprical Data")



rmse.baet<- sqrt(mean((cor.df$V2 - cor.df$mean.abund)^2))

rmse.baet.scale <- sqrt(mean((scale(cor.df$V2) - scale(cor.df$mean.abund))^2))
# coverage 
coverage <- mean(scale(cor.df$V2) >= (scale(cor.df$mean.abund) - (1.96*rmse.baet.scale)) & scale(cor.df$V2) <= (scale(cor.df$mean.abund) + (1.96*rmse.baet.scale)))

BAETts <- ggplot(data = means.list.BAET, aes(x = `temps$dts`,  y = scale(mean.abund), group = 1, color = "Model")) +
  geom_ribbon(aes(ymin = scale(mean.abund) - 1.96 * rmse.baet.scale,
                  ymax = scale(mean.abund) + 1.96 * rmse.baet.scale),
              colour = 'transparent',
              fill = "black",
              alpha = .1,
              show.legend = F) +
  geom_line(show.legend = T, linewidth = 1, alpha = 0.8) +
  geom_line(data =BAET.samp.sum, aes(x = as.Date(V1, origin = "1970-01-01"), y = scale(V2), color = "Empirical"), linewidth = 1,  show.legend = T, alpha = 0.8)+
  labs(y= "Scaled Abundance", title = expression(paste(italic("Baetidae"), " spp.") ))+
  geom_text(mapping = aes(x = as.Date("1998-01-01"), y =5, label = paste('rho', "==", 0.63)), parse = T, color = "black", size = 5.5)+
  geom_text(mapping = aes(x = as.Date("1998-01-01"), y =5.75, label = paste("C = 89%")), color = "black", size = 5.5)+
  geom_text(mapping = aes(x = as.Date("1998-01-01"), y =6.5, label = paste("Scaled RMSE = 0.96")), color = "black", size = 5.5)+
    xlab("")+
  ylim(c(-4,7))+
  labs(colour=" ")+
  theme_bw()+
  scale_color_manual(values = colors)+
  #scale_y_continuous(
    # sec.axis = sec_axis(~., name="Baetidae Larvae (inds/m2)"
    # ))+
  theme(text = element_text(size = 15), axis.text.x = element_text(angle=45, hjust = 1, size = 15), 
        axis.text.y = element_text(size = 15), )+
  scale_x_date(date_labels="%Y")

# Thats nice but maybe n mix will be better
# we have discharge from 1956 on
discharge <- readNWISdv("09234500", "00060", "1986-05-01", "1999-10-21", statCd = "00003")
# Bankfull discharge for Green River from https://digitalcommons.usu.edu/cgi/viewcontent.cgi?article=7604&context=etd
flow.magnitude <- TimestepDischarge(discharge, 22424.813)
flow.magnitude$dts <- as.Date(flow.magnitude$dts)
# missing a big chunk of temperature data between the 60s and mid-80s 
# and all throughout so we need to start in the 'mid 80s so we can have paired temp, flow, and invert dat
temp <- readNWISdv("09234500", "00010", "1986-05-01", "1999-10-21", statCd = "00003")
temp$Date <- as.Date(temp$Date)
# all dates
all_dates <- as.data.frame(seq.Date(from = as.Date("1986-05-01"), to = as.Date("1999-10-21"), by = "days"))
names(all_dates) <- "Date"
# join them
temp <- full_join(temp, all_dates)
# organize by date
temp <- temp[order(as.Date(temp$Date)),]
# get biweekly avgs - we still have some NAs left over and we can't have NAs in predictors 
temps <- TimestepTemperature(temp)
temps$dts <- as.Date(temps$dts)

# now for the issues of NAs, need to add in averages 
simtemp <- readNWISdv("09234500", "00010", "1986-07-01", "1999-10-21")
simtemp$Date <- as_datetime(simtemp$Date)
simtemp$Date <- yday(simtemp$Date)
simtemp$Temperature <- simtemp$X_00010_00003
simtemp <- simtemp %>% group_by(Date) %>% dplyr::summarise(Temperature = mean(Temperature))

temps$Temperature[which(is.na(temps$Temperature)==T)] <- simtemp$Temperature[yday(temps$dts[which(is.na(temps$Temperature)==T)])]
# that is our filled in temperature data

out <- BAETmodel(flow.data = flow.magnitude$Discharge, temp.data = temps, disturbanceK = 40000, baselineK = 5000, Qmin = 0.1, extinct = 50, iteration = 9, peaklist = 0.13, peakeach = length(temps$Temperature))


means.list.BAET <- mean.data.frame(out, burnin = 200, iteration= 9)
means.list.BAET <- cbind(means.list.BAET, temps$dts[199:352])
means.list.BAET$`temps$dts` <- as.Date(means.list.BAET$`temps$dts`)


BAETdata$Count <- round(as.numeric(BAETdata$Density)*as.numeric(BAETdata$Area))
vals <- vector()
for (i in 1:length(temps$dts)){
  d <- BAETdata[which(BAETdata$Date %within% interval(as.Date(temps$dts[i]), as.Date(temps$dts[i+1])-1) == T),]
  if (length(d$Count) > 0) {
    s<- rep(i, times = length(d$Count))
    vals <- append(vals, s)}
}
# add to data frame
BAETdata <- cbind(BAETdata, vals)



max_visits <- vector() # vector to put max obs per site
means <- vector()
for (i in 1:length(temps$dts)){  # cycle through all the 14 day timesteps that we have model output for
  # pull abundances between each each timestep
  d <- BAETdata[which(BAETdata$Date >= temps$dts[i] & BAETdata$Date < temps$dts[i+1]), ]
  max_visits[i] <- length(d$Date) # number of observations
  means[i] <- mean(d$Count, na.rm = T) # means - this is just for checking NAs
}

# make data frame with the timestep, the mean
df <- cbind.data.frame(temps$dts, means)
# remove anythig where there are NA values of Density (means that volume or count is missing)
df <- df[!is.na(df$means), ]
# phenology may also play into dynamics, so include month column as well
#month <- month(BAETdata$Date)
library(aspace)
df$circdate <- sin(as_radians((lubridate::yday(df$`temps$dts`)/365)*360))


# define our RxJ matrix
R <- length(temps$dts)
J <- max(max_visits)
site_mat <- matrix(data = NA, nrow = R, ncol = J)
dens_mat <- matrix(data = NA, nrow = R, ncol = J)
obs_intercept <- matrix(data = 1, nrow = R, ncol = J)
areas <- matrix(data = NA, nrow = R, ncol = J)
# make vector for flows at each timestep
# make RxJ matrix full of densities
# make RxJ matrix full of raw counts
# make RxJ matrix full of volumes sampled for each abundance
#time <- matrix(data = NA, nrow = R, ncol = J)
for (i in 1:length(temps$dts)){
  d <- BAETdata[which(BAETdata$Date >= temps$dts[i] & BAETdata$Date < temps$dts[i+1]), ]
  site_mat[i, ] <- c(d$Count, rep(NA, times = (J- length(d$Count))))
  areas[i, ] <- c((d$Area),rep(NA, times = (J- length(d$Area))))
}


# we need to remove all timesteps that are just NAs
nodata <- which(is.na(site_mat[,1]))
# first identify all the timsteps that don't have data (so we can match them up later)
site_mat <- as.matrix(site_mat[-nodata,]) # count data
#dens_mat <- as.matrix(dens_mat[-nodata, ]) # density data
obs_intercept <- as.matrix(obs_intercept[-nodata,]) # intercept for obs cov
#time <- as.matrix(scale(time[-nodata,])) # duration in H20 obs cov
#time[is.na(time)] <- mean(time, na.rm = TRUE) # replace NAs with mean duration time since NAs not allowed in predictors
#offset
class(areas) <- "numeric"
offset <- as.matrix(scale(log(areas[-nodata, ])))
offset[is.na(offset)] <- mean(offset, na.rm = TRUE) # replace NAs with mean 


flows <- as.data.frame(scale(flow.magnitude$Discharge[-nodata])) # site cov flow 
temperature <- as.data.frame(scale(temps$Temperature[-nodata])) # site cov temp
circdate <- as.data.frame(df$circdate[2:54])
# dimnames(time) <- list(temps$dts[-nodata], seq(1:48))
# time <- list(time)
# names(time) <- c("time")
site_intercept <- rep(1, times = length(flows$V1)) 
site_covs<- as.matrix(cbind(site_intercept, circdate)) #flows,temperature, circdate)
obs_covs <- array(data= NA, dim = c(length(flows$V1),J,1))
obs_covs[,,1] <- obs_intercept                                  

sink("N-mixturePoisBAET.jags")
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
    fit.rep <- sum(E.rep[,])
    } # End model
", fill = TRUE)
sink()

jags_data <- list(y = site_mat,
                  XN = site_covs,
                  Xp = obs_covs,
                  J = dim(site_mat)[2], #visits
                  R = dim(site_mat)[1], #sites
                  off = offset,
                  nAlpha = dim(site_covs)[2],
                  nBeta = dim(obs_covs)[3])

nAlpha <- dim(site_covs)[2]
nBeta <- dim(obs_covs)[3]
jags_inits <- function(){
  list(
    N = (apply(jags_data$y, 1, max, na.rm=TRUE)),
    alpha=runif(nAlpha,-1,1),
    beta=runif(nBeta,-1,1))}

parameters <- c("alpha", "beta", "lambda", "p", "N")


nc <- 3
ni <- 10000
nb <- 2500
nt <- 1

Nmix_fit <- jags.model("N-mixturePoisBAET.jags",data = jags_data, inits = jags_inits, n.chains = nc, n.adapt = 1000)

update(Nmix_fit, n.iter = 1000)

Nmix_fit_UI1 <- jagsUI::jags(data = jags_data, inits = jags_inits, parameters.to.save = parameters, model.file = "N-mixturePoisBAET.jags",  n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, parallel = TRUE)
# 


zm = coda.samples(Nmix_fit, variable.names = c("lambda"), n.iter = ni, n.thin = nt)


lam <- MCMCpstr(zm, "lambda")
# N <- MCMCpstr(zm, "N")
# N <- as.data.frame(cbind(as.Date(temps$dts[-nodata]), unlist(N)))
# N$V1 <- as.Date(N$V1, origin = "1970-01-01")


lam <- as.data.frame(cbind(as.Date(temps$dts[-nodata]), unlist(lam)))
lam$V1 <- as.Date(lam$V1, origin = "1970-01-01")

# plot(unlist(y.rep), unlist(site_mat))
# y.rep <- MCMCpstr(zm, "y.rep")
# exp <- MCMCpstr(zm, "exp")
# 
# 
# plot(unlist(y.rep), unlist(site_mat))
# abline(0, 1)
# 
cor.df <- left_join(lam, means.list.BAET, by=c('V1'="temps$dts"), copy = T)
cor.lm <- lm(cor.df$mean.abund ~ cor.df$V2)
Poislam <- cor.test((cor.df$V2), (cor.df$mean.abund), method = "spearman")
# 
# cor.df <- left_join(N, means.list.BAET, by = c("V1" = "temps$dts"), copy = T)
# cor.test(cor.df$V2, cor.df$mean.abund, method = "spearman")
sink("N-mixtureZIPBAET.jags")
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

Nmix_fit <- jags.model("N-mixtureZIPBAET.jags",data = jags_data, inits = jags_inits, n.chains = nc, n.adapt = 1000)

update(Nmix_fit, n.iter = 1000)

Nmix_fit_UI2 <- jagsUI::jags(data = jags_data, inits = jags_inits, parameters.to.save = parameters, model.file = "N-mixtureZIPBAET.jags",  n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, parallel = TRUE)
#
#print(Nmix_fit_UI)

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
# 
# plot(unlist(y.rep), unlist(site_mat))
# abline(0, 1)
# 
# fit <- MCMCchains(zm, "fit")
# fit.rep <- MCMCchains(zm, "fit.rep")
# mean(fit > fit.rep) # close to 1 so bad fit?
# plot(fit.rep ~ fit)
# abline(0, 1) # 1 to 1 line not even there
# 
# 
# fit_df <- data.frame(y = c(c(unlist(site_mat)), c(unlist(y.rep))),
#                      data = rep(c("Observed", "Simulated"), each = length(site_mat)))
# library(ggplot2)
# ggplot(fit_df, aes(x = y, fill = data)) + geom_histogram() + facet_grid(.~data)#still not getting all the 0s and missing the really high #s
# 
cor.df <- left_join(lam, means.list.BAET, by=c('V1'="temps$dts"), copy = T)
cor.lm <- lm(cor.df$mean.abund ~ cor.df$V2)
Ziplam <- cor.test((cor.df$V2), (cor.df$mean.abund), method = "spearman")

sink("N-mixtureZIPoverdispBAET.jags")
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

Nmix_fit <- jags.model("N-mixtureZIPoverdispBAET.jags",data = jags_data, inits = jags_inits, n.chains = nc, n.adapt = 1000)

update(Nmix_fit, n.iter = 1000)

Nmix_fit_UI3 <- jagsUI::jags(data = jags_data, inits = jags_inits, parameters.to.save = parameters, model.file = "N-mixtureZIPoverdispBAET.jags",  n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, parallel = TRUE)


zm = coda.samples(Nmix_fit, variable.names = c("lambda"), n.iter = ni, n.thin = nt)


lam <- MCMCpstr(zm, "lambda")
# p <- MCMCpstr(zm, "p")
# N <- MCMCpstr(zm, "N")
# N <- as.data.frame(cbind(as.Date(temps$dts[-nodata]), unlist(N)))
# N$V1 <- as.Date(N$V1, origin = "1970-01-01")


lam <- as.data.frame(cbind(as.Date(temps$dts[-nodata]), unlist(lam)))
lam$V1 <- as.Date(lam$V1, origin = "1970-01-01")

# plot(unlist(y.rep), unlist(site_mat))
# abline(0, 1)
# 
# 
# mean(fit > fit.rep) # close to 1 so bad fit?
# plot(fit.rep ~ fit)
# abline(0, 1) # 1 to 1 line not even there
# 
# 
# fit_df <- data.frame(y = c(c(unlist(site_mat)), c(unlist(y.rep))),
#                      data = rep(c("Observed", "Simulated"), each = length(site_mat)))
# library(ggplot2)
# ggplot(fit_df, aes(x = y, fill = data)) + geom_histogram() + facet_grid(.~data) #still not getting all the 0s and missing the really high #s
# 
cor.df <- left_join(lam, means.list.BAET, by=c('V1'="temps$dts"), copy = T)
cor.lm <- lm(cor.df$mean.abund ~ cor.df$V2)
Zip_ovdlam <- cor.test((cor.df$V2), (cor.df$mean.abund), method = "spearman")

sink("N-mixturePoisoverdispBAET.jags")
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

Nmix_fit <- jags.model("N-mixturePoisoverdispBAET.jags",data = jags_data, inits = jags_inits, n.chains = nc, n.adapt = 1000)

update(Nmix_fit, n.iter = 1000)

Nmix_fit_UI4 <- jagsUI::jags(data = jags_data, inits = jags_inits, parameters.to.save = parameters, model.file = "N-mixturePoisoverdispBAET.jags",  n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, parallel = TRUE)

zm = coda.samples(Nmix_fit, variable.names = c("lambda"), n.iter = ni, n.thin = nt)

lam <- MCMCpstr(zm, "lambda")
# p <- MCMCpstr(zm, "p")
# N <- MCMCpstr(zm, "N")
# N <- as.data.frame(cbind(as.Date(temps$dts[-nodata]), unlist(N)))
# N$V1 <- as.Date(N$V1, origin = "1970-01-01")
# 

lam <- as.data.frame(cbind(as.Date(temps$dts[-nodata]), unlist(lam)))
lam$V1 <- as.Date(lam$V1, origin = "1970-01-01")

# plot(unlist(y.rep), unlist(site_mat))
# abline(0, 1)
# 
# plot(unlist(site_mat), unlist(exp))
# 
# fit <- MCMCchains(zm, "fit")
# fit.rep <- MCMCchains(zm, "fit.rep")
# mean(fit > fit.rep) # close to 1 so bad fit?
# plot(fit.rep ~ fit)
# abline(0, 1) # 1 to 1 line not even there

cor.df <- left_join(lam, means.list.BAET, by=c('V1'="temps$dts"), copy = T)
cor.lm <- lm(cor.df$mean.abund ~ cor.df$V2)
Pois_ovdlam <- cor.test((cor.df$V2), (cor.df$mean.abund), method = "spearman")

# fit_df <- data.frame(y = c(c(unlist(site_mat)), c(unlist(y.rep))),
#                      data = rep(c("Observed", "Simulated"), each = length(site_mat)))
# library(ggplot2)
# ggplot(fit_df, aes(x = y, fill = data)) + geom_histogram() + facet_grid(.~data) #still not getting all the 0s and missing the really high #s

sink("N-mixtureNBBAET.jags")
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

Nmix_fit <- jags.model("N-mixtureNBBAET.jags",data = jags_data, inits = jags_inits, n.chains = nc, n.adapt = 1000)

update(Nmix_fit, n.iter = 1000)



zm = coda.samples(Nmix_fit, variable.names = c("lambda"), n.iter = ni, n.thin = nt)


lam <- MCMCpstr(zm, "lambda")
# N <- MCMCpstr(zm, "N")
# N <- as.data.frame(cbind(as.Date(temps$dts[-nodata]), unlist(N)))
# N$V1 <- as.Date(N$V1, origin = "1970-01-01")
# 

lam <- as.data.frame(cbind(as.Date(temps$dts[-nodata]), unlist(lam)))
lam$V1 <- as.Date(lam$V1, origin = "1970-01-01")

# 
# 
# th <- MCMCchains(zm, "theta")
# s <- MCMCchains(zm, "s")
# fit <- MCMCchains(zm, "fit")
# fit_new <- MCMCchains(zm, "fit_new")
# yhat <- MCMCchains(zm, "yhat")
# y_new <- MCMCchains(zm, "y_new")
# dim(lam)
# dim(th)
# 
# y.rep <- MCMCpstr(zm, "yhat")
# exp <- MCMCpstr(zm, "y_new")
# 
# 
# plot(unlist(y.rep), unlist(site_mat))
# abline(0, 1)
# 
# plot(unlist(site_mat), unlist(exp))
# 
# fit <- MCMCchains(zm, "fit")
# fit.rep <- MCMCchains(zm, "fit_new")
# mean(fit > fit.rep) # close to 1 so bad fit?
# plot(fit.rep ~ fit)
# abline(0, 1) # 1 to 1 line not even there
# 
# 
# fit_df <- data.frame(y = c(c(unlist(site_mat)), c(unlist(y.rep))),
#                      data = rep(c("Observed", "Simulated"), each = length(site_mat)))
# library(ggplot2)
# ggplot(fit_df, aes(x = y, fill = data)) + geom_histogram() + facet_grid(.~data) #still not getting all the 0s and missing the really high #s


Nmix_fit_UI5 <- jagsUI::jags(data = jags_data, inits = jags_inits, parameters.to.save = parameters, model.file = "N-mixtureNBBAET.jags",  n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, parallel = TRUE)


cor.df <- left_join(lam, means.list.BAET, by=c('V1'="temps$dts"), copy = T)
cor.lm <- lm(cor.df$mean.abund ~ cor.df$V2)
nblam <- cor.test((cor.df$V2), (cor.df$mean.abund), method = "spearman")

