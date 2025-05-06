######################################
## Code to Validate 
######################################

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
library(rjags)
library(MCMCvis)

source("1spFunctions.R")
source("HYOS_1sp.R")
# load Water Temperature data from above Diamond Creek Confluence (RM226)
temp <- read.delim("CRaboveDC_Temp.tsv", header=T)
colnames(temp) <- c("Date", "Temperature")
temp <- subset(temp, Date >= "2004-05-23")
temp$Date <- as.Date(temp$Date, format = "%Y-%m-%d")
temp <- aggregate(temp$Temperature, by = list(temp$Date), FUN = mean)
colnames(temp) <- c("Date", "Temperature")
simtemp <- temp
# all dates
all_dates <- as.data.frame(seq.Date(from = as.Date("2004-05-23"), to = as.Date("2024-08-23"), by = "days"))
names(all_dates) <- "Date"
# join them
temp <- full_join(temp, all_dates)


simtemp$Date <- yday(simtemp$Date)
simtemp <- simtemp %>% group_by(Date) %>% dplyr::summarise(Temperature = mean(Temperature))


temp$Temperature[which(is.na(temp$Temperature)==T)] <- simtemp$Temperature[yday(temp$Date[which(is.na(temp$Temperature)==T)])]
temp <-temp[order(temp$Date),]
temps <- TimestepTemperature(temp) # intermittantly missing data until  2000-12-26 
# load discharge data from above Diamond Creek Confluence (RM226)
# now for the issues of NAs, need to add in averages 

# that is our filled in temperature data

discharge <- readNWISdv("09404200", "00060", "2004-05-22", "2024-08-23")
discharge <- full_join(discharge, all_dates)
flow.magnitude <- TimestepDischarge(discharge, 85000)

set.seed(111)
out <- HYOSmodel(flow.data = flow.magnitude$Discharge, temp.data = temps, disturbanceK = 40000, baselineK = 10000, Qmin = 0.3, extinct = 50, iteration = 1000, peaklist = 0.17, peakeach = length(temps$Temperature))

drift.data.total <- readDB(gear = "LightTrap", type = "Sample", updater = F)

drift.CR <- drift.data.total[which(drift.data.total$Region == "GrandCanyon" & drift.data.total$RiverMile >= 219 & drift.data.total$RiverMile <= 225),]


specieslist <- c("HYOS")
HYOS.samp.CR <- sampspec(samp = drift.CR, stats = T, species = specieslist)
HYOS.samp <- merge(HYOS.samp.CR$Statistics, HYOS.samp.CR$Samples, by = "BarcodeID", all = T)
HYOS.samp <- HYOS.samp[which(HYOS.samp$FlagStrange == F),]
HYOS.samp$Density <- HYOS.samp$CountTotal/HYOS.samp$TimeElapsed

#HYOS.samp <- HYOS.samp[-which(HYOS.samp$Habitat == "Rock"),]
HYOS.samp <- HYOS.samp[-which(HYOS.samp$Weather == "Rain"),]
HYOS.samp <- HYOS.samp[-which(HYOS.samp$Bats == T), ]
HYOS.samp <- aggregate(HYOS.samp$Density, list(HYOS.samp$Date), FUN = mean)
means <- vector()
for (i in 1:length(temps$dts)){
  d <- HYOS.samp[which(HYOS.samp$Group.1 >= temps$dts[i] & HYOS.samp$Group.1 < temps$dts[i+1]),]
  if (any(is.nan(mean(d$x))) == T || any(is.na((d$x) == T))) {
    s = NA
  } else {
    s<- mean(d$x)}
  means <- append(means, s)
  # we know the last value doesn't fit into interval and needs to be manually added
}
means[length(temps$dts)] <- HYOS.samp$x[length(HYOS.samp$x)]

repdf <- plyr::adply(out, c(1,2,3))
names(repdf) <- c('timesteps', 'stage', 'rep', 'abund')
repdf$timesteps <- as.numeric(as.character(repdf$timesteps))
repdf$timesteps <- as.factor(repdf$timesteps)
repdfhyos <- subset(repdf,stage == "S3")
means.list.HYOS<- repdfhyos %>%
  dplyr::group_by(timesteps) %>%
  dplyr::summarise(mean.abund = mean(abund),
                   sd.abund = sd(abund),
                   se.abund = sd(abund)/sqrt(1000)) %>%
  ungroup()



means.list.HYOS <- as.data.frame(cbind(means.list.HYOS[200:530,], temps$dts[199:529]))
colnames(means.list.HYOS) <- c("timesteps", "mean.abund", "sd.abund", "se.abund", "Date")
means.list.HYOS$Date <- as.Date(as.POSIXct(means.list.HYOS$Date, origin = "1970-01-01"))


HYOS.samp.sum <- na.omit(as.data.frame(cbind(as.Date(means.list.HYOS$Date), means[200:530])))
HYOS.samp.sum$V1 <- as.Date(HYOS.samp.sum$V1, origin = "1970-01-01")
HYOS.samp.sum <- HYOS.samp.sum[which(HYOS.samp.sum$V1 < "2022-01-01"),]

#test for temporal autocorrelation
#acf(na.omit(means)) # we have some

# use data splitting method
HYOS.samp.sum1 <- HYOS.samp.sum %>% slice(which(row_number() %% 5 == 0))
HYOS.samp.sum2 <- HYOS.samp.sum %>%  slice(which(row_number() %% 5 == 1))
HYOS.samp.sum3 <- HYOS.samp.sum %>%  slice(which(row_number() %% 5 == 2))
HYOS.samp.sum4 <- HYOS.samp.sum %>%  slice(which(row_number() %% 5 == 3))
HYOS.samp.sum5 <- HYOS.samp.sum %>%  slice(which(row_number() %% 5 == 4))

cor.df <- left_join(HYOS.samp.sum, means.list.HYOS, by=c('V1'="Date"), copy = T)
# cor.lm <- lm(cor.df$mean.abund ~ cor.df$V2)
# 
# summary(cor.lm)

cor.df1 <- left_join(HYOS.samp.sum1, means.list.HYOS, by=c('V1'="Date"), copy = T)
#cor.lm1 <- lm(cor.df1$mean.abund ~ cor.df1$V2)
# summary(cor.lm1)
# plot(cor.df1$V2, cor.df1$mean.abund)
rho1 <- cor.test((cor.df1$V2), (cor.df1$mean.abund), method = "spearman")

cor.df2 <- left_join(HYOS.samp.sum2, means.list.HYOS, by=c('V1'="Date"), copy = T)
#cor.lm2 <- lm(cor.df2$mean.abund ~ cor.df2$V2)
rho2 <- cor.test((cor.df2$V2), (cor.df2$mean.abund), method = "spearman")

cor.df3 <- left_join(HYOS.samp.sum3, means.list.HYOS, by=c('V1'="Date"), copy = T)
# cor.lm3 <- lm(cor.df3$mean.abund ~ cor.df3$V2)
# summary(cor.lm3)
rho3 <- cor.test((cor.df3$V2), (cor.df3$mean.abund), method = "spearman")

cor.df4 <- left_join(HYOS.samp.sum4, means.list.HYOS, by=c('V1'="Date"), copy = T)
#cor.lm2 <- lm(cor.df2$mean.abund ~ cor.df2$V2)
rho4 <- cor.test((cor.df4$V2), (cor.df4$mean.abund), method = "spearman")

cor.df5 <- left_join(HYOS.samp.sum5, means.list.HYOS, by=c('V1'="Date"), copy = T)
#cor.lm2 <- lm(cor.df2$mean.abund ~ cor.df2$V2)
rho5 <- cor.test((cor.df5$V2), (cor.df5$mean.abund), method = "spearman")



rho <- mean(c(rho1$estimate, rho2$estimate, rho3$estimate, rho4$estimate, rho5$estimate))
sd <- sd(c(rho1$estimate, rho2$estimate, rho3$estimate))
# 
# ggplot(data = cor.df, aes(x = (V2) , y = (mean.abund)))+
#   geom_point()+
#   stat_smooth(method = "lm",
#               formula = y ~ x,
#               geom = "smooth")+
#   geom_text(x = 3, y = 300, label = "")+
#   labs(y = "Hydropsychidae Model Output", x = "Hydropsychidae Empirical Data")
# 
#cor.test((cor.df$V2), (cor.df$mean.abund), method = "spearman")


#hist(cor.df$V2, xlab = "Hydropsychidae Adults (#/hour)", col = "#CADBD7")

colors <- c("#AA3377", "black" )

rmse.hyos <- sqrt(mean((cor.df$V2 - cor.df$mean.abund)^2))

rmse.hyos.scale <- sqrt(mean((scale(cor.df$V2) - scale(cor.df$mean.abund))^2))
coverage <- mean(scale(cor.df$V2) >= (scale(cor.df$mean.abund) - (1.96*rmse.hyos.scale)) & scale(cor.df$V2) <= (scale(cor.df$mean.abund) + (1.96*rmse.hyos.scale)))

HYOSts <- ggplot(data = cor.df, aes(x = V1, y = scale(mean.abund), group = 1, color = "Model")) +
  geom_ribbon(aes(ymin = scale(mean.abund) - (1.96 * rmse.hyos.scale),
                  ymax = scale(mean.abund) + (1.96 * rmse.hyos.scale)), 
  alpha = .1,
  fill = "black",
  color = "transparent",
  show.legend = F) +
  geom_line(show.legend = T, linewidth = 1, alpha = 0.8) +
  geom_line(data =cor.df, aes(x = V1, y = scale(V2), color = "Empirical"), linewidth = 1, alpha = 0.8,show.legend = T)+
  # geom_line(data = flow.magnitude, aes(x = as.Date(dts), y = Discharge), color = "blue") +
  # geom_line(data = temps, aes(x = as.Date(dts), y = Temperature*0.1), color = "green")+
  # #coord_cartesian(ylim = c(0,2000000)) +
  labs(y="Scaled Abundance")+
  geom_text(mapping = aes(x = as.Date("2015-06-01"), y =5, label = paste('rho', "==", 0.69)), parse = T, color = "black", size = 5.5)+
  geom_text(mapping = aes(x = as.Date("2015-06-01"), y =5.75, label = paste('C = 93%')), color = "black", size = 5.5)+
  geom_text(mapping = aes(x = as.Date("2015-06-01"), y =6.5, label = paste('Scaled RMSE = 1.19')), color = "black", size = 5.5)+
    #scale_y_continuous(
    # Features of the first axis
    # Add a second axis and specify its features
  #   sec.axis = sec_axis(~., name="Hydropsychidae Adults (ind/hour)")
  # ) + 
  xlab(" ")+
  ylim(c(-4,7))+
  labs(colour=" ", title = expression(paste(italic("Hydropsyche"), " spp.")))+
  theme_bw()+
  theme(text = element_text(size = 15), axis.text.x = element_text(angle=45, hjust = 1, size = 15), 
        axis.text.y = element_text(size = 15))+
  scale_x_date(date_labels="%Y")+
  scale_color_manual(values = colors)

# plot(means.list.HYOS$Date, means.list.HYOS$mean.abund, type = "l")
# lines(as.Date(temps$dts), temps$Temperature * 5, col = "red")
# lines(as.Date(flow.magnitude$dts), flow.magnitude$Discharge * 1000, col = "blue")
# lines(as.Date(HYOS.samp.sum$V1), HYOS.samp.sum$V2*2, col = "darkgreen")

drift.data.total <- readDB(gear = "LightTrap", type = "Sample", updater = F)

drift.CR <- drift.data.total[which(drift.data.total$Region == "GrandCanyon" & drift.data.total$RiverMile >= 219 & drift.data.total$RiverMile <= 225),]

specieslist <- c("HYOS")
HYOS.samp.CR <- sampspec(samp = drift.CR, stats = T, species = specieslist)
HYOS.samp <- merge(HYOS.samp.CR$Statistics, HYOS.samp.CR$Samples, by = "BarcodeID", all = T)
HYOS.samp <- HYOS.samp[which(HYOS.samp$FlagStrange == F),]

vals <- vector()
for (i in 1:length(temps$dts)){
  d <- HYOS.samp[which(HYOS.samp$Date %within% interval(temps$dts[i], temps$dts[i+1]-1) == T),]
  if (length(d$CountTotal) > 0) {
    s<- rep(i, times = length(d$CountTotal))
    vals <- append(vals, s)}
}

# add to data frame
HYOS.samp <- cbind(HYOS.samp, vals)

# now we need into include mean water temperature and discharge for each 
HYOS.samp <- cbind(HYOS.samp, temps$Temperature[HYOS.samp$vals])
HYOS.samp <- cbind(HYOS.samp, flow.magnitude$Discharge[HYOS.samp$vals])

HYOS.samp <- HYOS.samp[-which(HYOS.samp$Habitat == "Rock"), ]
HYOS.samp$Habitat <- as.character(HYOS.samp$Habitat)
HYOS.samp$Habitat[which(HYOS.samp$Habitat == "Sand")] <- 0
HYOS.samp$Habitat[which(HYOS.samp$Habitat == "Vegetation")] <- 1
HYOS.samp$Habitat <- as.numeric(HYOS.samp$Habitat)
#HYOS.samp <- HYOS.samp[-which(HYOS.samp$Weather == "Rain"), ]

HYOS.samp$Weather <- as.character(HYOS.samp$Weather)
HYOS.samp$Weather[which(HYOS.samp$Weather == "Rain")] <- 1
HYOS.samp$Weather[which(HYOS.samp$Weather == "NoRain")] <- 0
HYOS.samp$Weather <- as.numeric(HYOS.samp$Weather)

HYOS.samp <- HYOS.samp[-which(HYOS.samp$Bats == T), ]
max_visits <- vector() # vector to put max obs per site
means <- vector()
for (i in 1:length(temps$dts)){  # cycle through all the 14 day timesteps that we have model output for
  # pull abundances between each each timestep
  d <- HYOS.samp[which(HYOS.samp$Date >= temps$dts[i] & HYOS.samp$Date < temps$dts[i+1]), ]
  max_visits[i] <- length(d$Date) # number of observations
  means[i] <- mean(d$CountTotal) # means - this is just for checking NAs
}
# make data frame with the timestep, the mean
df <- cbind.data.frame(temps$dts, means)
# remove anythig where there are NA values of Density (means that volume or count is missing)
df <- df[!is.na(df$means), ]
# phenology may also play into dynamics, so include month column as well
#month <- month(HYOS.samp$Date)
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
time <- matrix(data = NA, nrow = R, ncol = J)
windspeed <- matrix(data = NA, nrow = R, ncol = J)
# habitat <- matrix(data = NA, nrow = R, ncol = J)
# weather <- matrix(data = NA, nrow = R, ncol = J)
for (i in 1:length(temps$dts)){
  d <- HYOS.samp[which(HYOS.samp$Date >= temps$dts[i] & HYOS.samp$Date < temps$dts[i+1]), ]
  flows[i] <- mean(d$`flow.magnitude$Discharge[HYOS.samp$vals]`)
  #date <- df[]
  temperature[i] <- mean(d$`temps$Temperature[HYOS.samp$vals]`)
  windspeed[i, ] <- c(d$WindSpeed, rep(NA, times = (J- length(d$CountTotal))))
  site_mat[i, ] <- c(d$CountTotal, rep(NA, times = (J- length(d$CountTotal))))
  # habitat[i, ] <- c(d$Habitat, rep(NA, times = (J- length(d$CountTotal))))
  time[i, ] <- c(d$TimeElapsed,rep(NA, times = (J- length(d$CountTotal))))
  # weather[i, ] <- c(d$Weather, rep(NA, times = (J- length(d$CountTotal))))
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
windspeed <- as.matrix((scale(windspeed[-nodata,])))
windspeed[is.na(windspeed)] <- mean(windspeed, na.rm = TRUE) # replace NAs with mean duration time since NAs not allowed in predictors or offsets
# habitat <- as.matrix(habitat[-nodata,])
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
#habitat[is.na(habitat)] <- Mode(na.omit(habitat)) # replace NAs with mean duration time since NAs not allowed in predictors or offsets
# weather <- as.matrix(weather[-nodata,])
# weather[is.na(weather)] <- Mode(na.omit(weather))

site_intercept <- rep(1, times = length(flows$V1)) 
site_covs<- as.matrix(cbind(site_intercept, flows, circdate)) #flows,temperature, circdate)
obs_covs <- array(data= NA, dim = c(length(flows$V1),J,2))
obs_covs[,,1] <- obs_intercept                                  
obs_covs[,,2] <- windspeed

#offset
offset <- as.matrix(scale(log(time[-nodata, ])))
offset[is.na(offset)] <- mean(offset, na.rm = TRUE) # replace NAs with mean duration time since NAs not allowed in predictors or offsets

sink("N-mixturePoisHYOS.jags")
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

Nmix_fit <- jags.model("N-mixturePoisHYOS.jags",data = jags_data, inits = jags_inits, n.chains = nc, n.adapt = 1000)

update(Nmix_fit, n.iter = 1000)

Nmix_fit_UI1 <- jagsUI::jags(data = jags_data, inits = jags_inits, parameters.to.save = parameters, model.file = "N-mixturePoisHYOS.jags",  n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, parallel = TRUE)
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
# mean(fit > fit.rep) # close to 1 so bad fit?
# plot(fit.rep ~ fit)
# abline(0, 1) # 1 to 1 line not even there
# 

# fit_df <- data.frame(y = c(c(unlist(site_mat)), c(unlist(y.rep))),
#                      data = rep(c("Observed", "Simulated"), each = length(site_mat)))
# library(ggplot2)
# ggplot(fit_df, aes(x = y, fill = data)) + geom_histogram() + facet_grid(.~data)# not getting all the 0s and missing the really high #s
# 
# cor.df <- left_join(N, means.list.HYOS, by=c('V1'="Date"), copy = T)
# cor.lm <- lm(cor.df$mean.abund ~ cor.df$V2)
# PoisN <- cor.test((cor.df$V2), (cor.df$mean.abund), method = "spearman")

# use data splitting method

cor.df <- left_join(lam, means.list.HYOS, by=c('V1'="Date"), copy = T)

cor.df1 <- cor.df  %>% slice(which(row_number() %% 5 == 0))
rho1 <- cor.test((cor.df1$V2), (cor.df1$mean.abund), method = "spearman")

cor.df2<- cor.df %>%  slice(which(row_number() %% 5 == 1))
rho2 <- cor.test((cor.df2$V2), (cor.df2$mean.abund), method = "spearman")

cor.df3 <- cor.df %>%  slice(which(row_number() %% 5 == 2))
rho3 <- cor.test((cor.df3$V2), (cor.df3$mean.abund), method = "spearman")

cor.df4 <- cor.df %>%  slice(which(row_number() %% 5 == 3))
rho4 <- cor.test((cor.df4$V2), (cor.df4$mean.abund), method = "spearman")

cor.df5 <- cor.df %>%  slice(which(row_number() %% 5 == 4))
rho5 <- cor.test((cor.df5$V2), (cor.df5$mean.abund), method = "spearman")

Poislam <- mean(c(rho1$estimate, rho2$estimate, rho3$estimate, rho4$estimate, rho5$estimate))

rm(zm)
rm(Nmix_fit)
# rm(Nmix_fit_UI1)
sink("N-mixtureZIPHYOS.jags")
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

Nmix_fit <- jags.model("N-mixtureZIPHYOS.jags",data = jags_data, inits = jags_inits, n.chains = nc, n.adapt = 1000)

update(Nmix_fit, n.iter = 1000)

Nmix_fit_UI2 <- jagsUI::jags(data = jags_data, inits = jags_inits, parameters.to.save = parameters, model.file = "N-mixtureZIPHYOS.jags",  n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, parallel = TRUE)
#
print(Nmix_fit_UI2)

zm = coda.samples(Nmix_fit, variable.names = c("lambda"), n.iter = ni, n.thin = nt)


lam <- MCMCpstr(zm, "lambda")
# N <- MCMCpstr(zm, "N")
# N <- as.data.frame(cbind(as.Date(temps$dts[-nodata]), unlist(N)))
# N$V1 <- as.Date(N$V1, origin = "1970-01-01")


lam <- as.data.frame(cbind(as.Date(temps$dts[-nodata]), unlist(lam)))
lam$V1 <- as.Date(lam$V1, origin = "1970-01-01")

# 
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
# cor.df <- left_join(N, means.list.HYOS, by=c('V1'="Date"), copy = T)
cor.df <- left_join(lam, means.list.HYOS, by=c('V1'="Date"), copy = T)

cor.df1 <- cor.df  %>% slice(which(row_number() %% 5 == 0))
rho1 <- cor.test((cor.df1$V2), (cor.df1$mean.abund), method = "spearman")

cor.df2<- cor.df %>%  slice(which(row_number() %% 5 == 1))
rho2 <- cor.test((cor.df2$V2), (cor.df2$mean.abund), method = "spearman")

cor.df3 <- cor.df %>%  slice(which(row_number() %% 5 == 2))
rho3 <- cor.test((cor.df3$V2), (cor.df3$mean.abund), method = "spearman")

cor.df4 <- cor.df %>%  slice(which(row_number() %% 5 == 3))
rho4 <- cor.test((cor.df4$V2), (cor.df4$mean.abund), method = "spearman")

cor.df5 <- cor.df %>%  slice(which(row_number() %% 5 == 4))
rho5 <- cor.test((cor.df5$V2), (cor.df5$mean.abund), method = "spearman")

Ziplam <- mean(c(rho1$estimate, rho2$estimate, rho3$estimate, rho4$estimate, rho5$estimate))

rm(zm)
rm(Nmix_fit)
# rm(Nmix_fit_UI2)

sink("N-mixtureZIPoverdispHYOS.jags")
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

Nmix_fit <- jags.model("N-mixtureZIPoverdispHYOS.jags",data = jags_data, inits = jags_inits, n.chains = nc, n.adapt = 1000)

update(Nmix_fit, n.iter = 1000)

Nmix_fit_UI3 <- jagsUI::jags(data = jags_data, inits = jags_inits, parameters.to.save = parameters, model.file = "N-mixtureZIPoverdispHYOS.jags",  n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, parallel = TRUE)

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
#
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
# cor.df <- left_join(N, means.list.HYOS, by=c('V1'="Date"), copy = T)
# cor.lm <- lm(cor.df$mean.abund ~ cor.df$V2)
# Zip_ovdN <- cor.test((cor.df$V2), (cor.df$mean.abund), method = "spearman")
#
cor.df <- left_join(lam, means.list.HYOS, by=c('V1'="Date"), copy = T)
# cor.lm <- lm(cor.df$mean.abund ~ cor.df$V2)
cor.df1 <- cor.df  %>% slice(which(row_number() %% 5 == 0))
rho1 <- cor.test((cor.df1$V2), (cor.df1$mean.abund), method = "spearman")

cor.df2<- cor.df %>%  slice(which(row_number() %% 5 == 1))
rho2 <- cor.test((cor.df2$V2), (cor.df2$mean.abund), method = "spearman")

cor.df3 <- cor.df %>%  slice(which(row_number() %% 5 == 2))
rho3 <- cor.test((cor.df3$V2), (cor.df3$mean.abund), method = "spearman")

cor.df4 <- cor.df %>%  slice(which(row_number() %% 5 == 3))
rho4 <- cor.test((cor.df4$V2), (cor.df4$mean.abund), method = "spearman")

cor.df5 <- cor.df %>%  slice(which(row_number() %% 5 == 4))
rho5 <- cor.test((cor.df5$V2), (cor.df5$mean.abund), method = "spearman")
Zip_ovdlam <- mean(c(rho1$estimate, rho2$estimate, rho3$estimate, rho4$estimate, rho5$estimate))

rm(zm)
rm(Nmix_fit)
# rm(Nmix_fit_UI3)

sink("N-mixturePoisoverdispHYOS.jags")
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

Nmix_fit <- jags.model("N-mixturePoisoverdispHYOS.jags",data = jags_data, inits = jags_inits, n.chains = nc, n.adapt = 1000)

update(Nmix_fit, n.iter = 1000)

Nmix_fit_UI4 <- jagsUI::jags(data = jags_data, inits = jags_inits, parameters.to.save = parameters, model.file = "N-mixturePoisoverdispHYOS.jags",  n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, parallel = TRUE)

print(Nmix_fit_UI4)

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

# cor.df <- left_join(N, means.list.HYOS, by=c('V1'="Date"), copy = T)
# cor.lm <- lm(cor.df$mean.abund ~ cor.df$V2)
# Pois_ovdN <- cor.test((cor.df$V2), (cor.df$mean.abund), method = "spearman")
# 
# 
cor.df <- left_join(lam, means.list.HYOS, by=c('V1'="Date"), copy = T)
cor.df1 <- cor.df  %>% slice(which(row_number() %% 5 == 0))
rho1 <- cor.test((cor.df1$V2), (cor.df1$mean.abund), method = "spearman")

cor.df2<- cor.df %>%  slice(which(row_number() %% 5 == 1))
rho2 <- cor.test((cor.df2$V2), (cor.df2$mean.abund), method = "spearman")

cor.df3 <- cor.df %>%  slice(which(row_number() %% 5 == 2))
rho3 <- cor.test((cor.df3$V2), (cor.df3$mean.abund), method = "spearman")

cor.df4 <- cor.df %>%  slice(which(row_number() %% 5 == 3))
rho4 <- cor.test((cor.df4$V2), (cor.df4$mean.abund), method = "spearman")

cor.df5 <- cor.df %>%  slice(which(row_number() %% 5 == 4))
rho5 <- cor.test((cor.df5$V2), (cor.df5$mean.abund), method = "spearman")
Pois_ovdlam <- mean(c(rho1$estimate, rho2$estimate, rho3$estimate, rho4$estimate, rho5$estimate))

# fit_df <- data.frame(y = c(c(unlist(site_mat)), c(unlist(y.rep))),
#                      data = rep(c("Observed", "Simulated"), each = length(site_mat)))
# library(ggplot2)
# ggplot(fit_df, aes(x = y, fill = data)) + geom_histogram() + facet_grid(.~data) #still not getting all the 0s and missing the really high #s
rm(zm)
rm(Nmix_fit)
# rm(Nmix_fit_U4)
sink("N-mixtureNBHYOS.jags")
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

Nmix_fit <- jags.model("N-mixtureNBHYOS.jags",data = jags_data, inits = jags_inits, n.chains = nc, n.adapt = 1000)

update(Nmix_fit, n.iter = 1000)



zm = coda.samples(Nmix_fit, variable.names = c("lambda"), n.iter = ni, n.thin = nt)


lam <- MCMCpstr(zm, "lambda")
# N <- MCMCpstr(zm, "N")
# N <- as.data.frame(cbind(as.Date(temps$dts[-nodata]), unlist(N)))
# N$V1 <- as.Date(N$V1, origin = "1970-01-01")


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
# 

Nmix_fit_UI5 <- jagsUI::jags(data = jags_data, inits = jags_inits, parameters.to.save = parameters, model.file = "N-mixtureNBHYOS.jags",  n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, parallel = TRUE)

print(Nmix_fit_UI5)

cor.df <- left_join(lam, means.list.HYOS, by=c('V1'="Date"), copy = T)
cor.df1 <- cor.df  %>% slice(which(row_number() %% 5 == 0))
rho1 <- cor.test((cor.df1$V2), (cor.df1$mean.abund), method = "spearman")

cor.df2<- cor.df %>%  slice(which(row_number() %% 5 == 1))
rho2 <- cor.test((cor.df2$V2), (cor.df2$mean.abund), method = "spearman")

cor.df3 <- cor.df %>%  slice(which(row_number() %% 5 == 2))
rho3 <- cor.test((cor.df3$V2), (cor.df3$mean.abund), method = "spearman")

cor.df4 <- cor.df %>%  slice(which(row_number() %% 5 == 3))
rho4 <- cor.test((cor.df4$V2), (cor.df4$mean.abund), method = "spearman")

cor.df5 <- cor.df %>%  slice(which(row_number() %% 5 == 4))
rho5 <- cor.test((cor.df5$V2), (cor.df5$mean.abund), method = "spearman")
nblam <- mean(c(rho1$estimate, rho2$estimate, rho3$estimate, rho4$estimate, rho5$estimate))


