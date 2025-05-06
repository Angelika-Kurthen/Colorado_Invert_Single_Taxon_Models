#####################################
# Ramped hydropeaking index, temperature, and flows
#####################################


library(purrr)
library(tidyverse)
library(lubridate)
library(plyr)
library(dplyr)
library(ggplot2)
# data retrieval tool from USGS
library(dataRetrieval)
source("1spFunctions.R")

# Code for HPC - tidyverse has some issues on our HPC because one of the packages is deprecated
# We have to manually load all tidyverse packages
# library(purrr, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
# library(tibble, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
# library(tidyr, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
# library(readr, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
# library(stringr, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
# library(forcats, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
# library(lubridate, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
# library(plyr, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
# library(dplyr, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
# library(ggplot2, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
# library(dataRetrieval, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")



#read in flow data from USGS gauge at Lees Ferry, AZ between 1985 to the end of the last water year
#flow <- readNWISdv("09380000", "00060", "1985-10-01", "2021-09-30")
#flow.magnitude <- TimestepDischarge(flow, 85000) #discharge data from USGS is in cfs, bankfull dischage (pers comm TK) 85000 cfs
# read in temp data
#temp <- temp <- readNWISdv("09380000", "00010", "2007-10-01", "2021-09-30")
# temp <- read.delim("gcmrc20230123125915.tsv", header=T)
# colnames(temp) <- c("Date", "Temperature")
# # peaklist <- c(0.01, 0.1, 0.2, 0.5)
# # tempslist <- c(0, 0.5, 1, 2.5, 5, 7.5)
# 
# n <- 77
# 
# # qr is the temp ramps I want to increase the average temp by 
# qr <- 0
# # how many years I want each temp ramp to last
# r <- 77
# 
# temps <- average.yearly.temp(temp, "Temperature", "Date")
# 
# temps <- rep.avg.year(temps, n, change.in.temp = qr, years.at.temp = r)
# Time <- c(1:1825)
# Date <- rep(c(1:365), times = 5)
# Day <- seq(as.Date("2022-01-01"), as.Date("2026-12-31"), by="days")
# Day <- Day[-which(Day == "2024-02-29")]
# 
# Temperature <-  -7.374528  * (cos(((2*pi)/365)*Date))  +  (-1.649263* sin(2*pi/(365)*Date))  + 10.956243
# 
# temp <- as.data.frame(cbind(Time, Day, Temperature))
# temp$Day <- as.Date(temp$Day, origin= "1970-01-01")
# colnames(temp) <- c("Time", "Date", "Temperature")
# temp <- TimestepTemperature(temp)
# temp <- temp[c(1,3)]
# peaklist <- 0 
# peakeach <- length(temp$Temperature)
# iteration <- 10
# baselineK <- 10000
# disturbanceK <- 40000
# extinct = 50
# flow.data <- discharge
# temp.data <- temp
# discharge <- rep(0.1, times = length(temp$Temperature))
# discharge[floor(runif(1, 90, 131))] <- runif(1, 0.25, 1)


# discharge <- rep(0.1, times = length(temps$Temperature))
HYOSmodel <- function(flow.data, temp.data, baselineK, disturbanceK, Qmin, extinct, iteration, peaklist = NULL, peakeach = NULL, stage_output = "all"){
#---------------------------------------------------------------
# set up model
source("HYOSSurvivorship.R")

Q <- as.numeric(flow.data)
temps <- temp.data

degreedays <- as.data.frame(cbind(temps$dts, temps$Temperature * 14))
colnames(degreedays) <- c("dts", "DegreeDay")
degreedays$DegreeDay[degreedays$DegreeDay<0] <- 0
degreedays$dts <- as.POSIXct(degreedays$dts, origin = "1970-01-01")

# need to make ramped increasing hydropeaking index 
hp <- c(rep(peaklist, each = peakeach))

# specify iterations
iterations <- iteration

# baseline K in the absence of disturbance
Kb <- as.numeric(baselineK)
# max K after a big disturbance
Kd <- as.numeric(disturbanceK)

# specify baseline transition probabilities for each species
# 3 stages - we have egg - larval instar V, pupae, and adult
# using eqs 3 to 7 from Birt et al 2009

# using Willis et al 1992, we consolidate life table infor for larval stages I and II (Stage 1),
# larval stages 
# and adults/eggs (Stage 3 with prebreeding census)

G1 = 0.016  # according to Willis et al 1992, average 3 timesteps (0.05/3)
G2 = 0.11/3  # move onto stage 3 (0.11/2) 2 instead of 20
P1 = 0.6666667  # remain in stage (1 1-(1/3))
P2 = 1 - (1/3) # remain in stage 2 (1 - (1/3)) 2 instead of 20


# want to run this for one year, in 14 day timesteps 
timestep <- seq(2, (length(temps$Temperature) + 1), by = 1)

# create an array to put our output into
output.N.array <- array(0, dim = c(length(timestep) + 1))

output.N.list <- list(output.N.array)

# create array to put the total N of all species into
Total.N <- array(0,
                 dim  <-c((length(timestep) +1 ), iterations),
                 dimnames <- list(1:(length(timestep) + 1), 1:iterations))

# create list of arrays w/ abundance data for each spp
reparray <- array(0,
                  
                  dim = c(length(timestep) + 1, 3, iterations),
                  dimnames = list(1:(length(timestep)+1), c("S1", "S2", "S3"), 1:iterations)
)

output.N.list <- reparray

Qmin <- Qmin
a <- 0.001
g <- 1
h <- surv.fit.HYOS$m$getPars()[2]  
k <- surv.fit.HYOS$m$getPars()[1] 

extinction <- extinct
#-------------------------
# Outer Loop of Iterations
#--------------------------


for (iter in c(1:iterations)) {
  K = Kb # need to reset K for each iteration
  
  # pull random values from a uniform distribution 
  output.N.list[1,1:3, iter]<- runif(3, min = 1, max = (0.3*K))
  
  # we often want to look at different parameter values after we run code, so we create some lists
  
  # list to input Ks
  Klist <- vector()
  Klist[1] <- K
  
  # list to imput flow morts
  flowmortlist <- vector()
  
  Flist <- vector()
  
  emergetime <- vector()
  
  sizelist <- vector()
  delta <- vector()
  TempSurvival <- vector()
  for(c in temps$Temperature){
    
    b <- TempSurv_HYOS(c)
    
    TempSurvival <- append(TempSurvival, b)
  }
  #-------------------------
  # Inner Loop of Timesteps
  #-------------------------
  
  for (t in timestep) {
    
    #----------------------------------------------------------
    # Calculate how many timesteps emerging adults have matured
    t <- t
    emergetime <- append(emergetime, back.count.degreedays(t, 1680, degreedays)) #mean from hauer and stanford 1728.889
    #delta <- append(delta, round(devtime(temps$Temperature[t-1])/14))
    
    #---------------------------------------------------------
    # Calculate fecundity per adult
    F3 = 235.6 * hydropeaking.mortality(lower = 0.4, upper = 0.6, h = hp[t-1])
    #F3 = rnorm(1, mean = 235.6, sd = 11.05102 ) * 0.5 * hydropeaking.mortality(lower = 0.4, upper = 0.6, h = hp[t-1])
    #from Willis Jr & Hendricks, sd calculated from 95% CI = 21.66 = 1.96*sd
    # * 0.5 assuming 50% female.
    
    # # we can scale fecundity based on the 95% CI of 21.66 (min = 213.94, max = 257.26) 
    if (t > 15) {
      size <- emergetime[t-1]
      sizelist <- append(sizelist, size)
      F3 <- ((7.219 * size) + 180.4) * hydropeaking.mortality(lower = 0.4, upper = 0.6, h = hp[t-1])
    }
    # size <- delta[t-1]
    # sizelist <- append(sizelist, size)
    # F3 <- ((41.86*size)+200) *0.5* hydropeaking.mortality(0.4, 0.6, h = hp[t-1])
    # 
    #---------------------------------------------------
    # Calculate the disturbance magnitude-K relationship 
    # Sets to 0 if below the Qmin
    Qf <- as.numeric(Qf.Function(Q[t-1], Qmin, a))
    
    #-------------------------------------------------------------------
    # Calculate K arrying capacity immediately following the disturbance
    # Calculate K arrying capacity immediately following the disturbance
    K0 <- as.numeric(K + ((Kd-K)*Qf))
    
    # Calculate final K for timestep, including relationship between K and time since disturbance
    K <- post.dist.K(K0, Kb, g, t, Q, Qmin)
    Klist <- append(Klist, K)
    
    
    #---------------------------------------------
    # Calculate effect of density dependence on fecundity
    
    # Logistic via Rogosch et al. Fish Model
    F3 <- Logistic.Dens.Dependence(F3, K, Total.N[t-1, iter])
    Flist <- append(Flist, F3)
    # 
    # Calculate new transition probabilities based on temperature
    # This is the growth v development tradeoff 
    # don't know if this exists for HYOS - they can emerge under a wide temp gradient (<5 - 25+ C) but relationship between growth and temp 
    # at cold temps, takes 20 timesteps to complete Stage 1
    # at warm temps, take 1 timestep to complete Stage 1
    # create linear eq to describe this, between temps of 5 and 25 C
    # Timesteps = -0.95(TEMP) + 24.75 
    
    #development measures# at cold temps
    # if (5 > temps$Temperature[t-1])  {
    #   G1 <- 0.36/20 *TempSurvival[t-1]
    #   P1 <- 1-(1/20)* TempSurvival[t-1]
    # }
    if (30 > temps$Temperature[t-1]){ # Spend about 20 timesteps in Stage 2
      G2 <- 0.04/20 * TempSurvival[t-1]
      P2 <- 1-(1/20) *TempSurvival[t-1]
    } else {
      G2 = 0.04/3 *TempSurvival[t-1]  # make them grow fast when super warm
      P2 = 1 - (1/3) *TempSurvival[t-1]
    }

    # if (temps$Temperature[t-1] > 30){
    #   G1 <- 0.36 *TempSurvival[t-1]
    #   P1 <- 0
    # }
    # Stage 1 can vary from a few timesteps to many (saw Instar 2 last from 5 to 250 days)
    if (5 <= temps$Temperature[t-1] & temps$Temperature[t-1] <=30 & is.na(emergetime[t] == F)){
      G1 <- 0.36/(emergetime[t-1]) *TempSurvival[t-1]
      P1 <- 1-(1/(emergetime[t-1])) *TempSurvival[t-1]
    }
    if (5 <= temps$Temperature[t-1] & temps$Temperature[t-1] <= 30 & is.na(emergetime[t] == T)) {
      G1 <- 0.36/((-0.95 * temps$Temperature[t-1]) + 24.75) *TempSurvival[t-1]
      P1 <- 1-(1/((-0.95 * temps$Temperature[t-1]) + 24.75))*TempSurvival[t-1]
    }
    
      # G1 <- 0.05/((delta[t-1]/2))
      # G2 <- 0.11/((delta[t-1]/2))
      # P1 <- 0.05*(1-(1/((delta[t-1]/2))))
      # P2 <- 0.11*(1-(1/((delta[t-1]/2))))

    #-----------------------------------------------
    # Create Lefkovitch Matrix
    
    H1 <- c(P1, 0, F3)
    H2 <- c(G1, P2, 0)
    H3 <- c(0, G2, 0) 
    
    A <- rbind( H1, H2, H3)
    
    #-----------------------------------------------
 
    # # growth (if below 10C, no growth can occur - everything basically freezes, if between 10 and 11, prob of remaining in same stage = 0.6395, if above 13, prob of transition to next stage is 0 
    
    #Glist <-append(Glist, AHYOS[3,2])
    #Plist <- append(Plist, AHYOS[2,2])
    
    #--------------------------------------
    # Calculate abundances for each stage
    
    output.N.list[t, 1:3, iter] <- A %*% output.N.list[t-1, 1:3, iter] 
    
    #------------------------------------------
    # Calculate immediate mortality due to temperature regime (outside of thermal optima)
    # output.N.list[t, 1, iter] <- output.N.list[t, 1, iter] * TempSurvival[t-1]
    # output.N.list[t, 2, iter] <- output.N.list[t, 2, iter] * TempSurvival[t-1]
    
    #Calculate immediate mortality due to flows
    # mortality due to flooding follows N0 = Nz*e^-hQ
    # but what about using a sigmoidal logistic function so we can control the threshold point and rate of growth
    # following m = 1/1+e^-h*(x-xf)
    # where h is is shape value
    # x is Q, and xf is threshold point (100% of pop dies)
    #plot(x = Q, y = 1/(1+exp(-0.02*(Q-100000))))
    
    #s1
    output.N.list[t, 1, iter] <- flood.mortality(output.N.list[t, 1, iter], k, h, Q[t-1], Qmin)
    #s2
    output.N.list[t,2,iter] <- flood.mortality(output.N.list[t,2,iter], k, h, Q[t-1], Qmin)
    #3
    #output.N.list[t,3,iter] <- flood.mortality(output.N.list[t,3,iter], k, h, Q[t-1], Qmin)
    
    flowmortlist <- append(flowmortlist, flood.mortality(1, k, h, Q[t-1], Qmin))
    #replist[[1]][,,1] <- output.N.list[[1]]
    Total.N[t,iter] <- sum(output.N.list[t,,iter])
    #check extinction threshold
    if (Total.N[t, iter] < extinction){
      output.N.list[t,,iter] <- 0
      Total.N[t, iter] <- 0  } #-------------------------
  }  # End Inner Loop  
  #------------------------- 
} #----------------------
# End Outer Loop
#----------------------
if (stage_output == "larvae"){
  return(output.N.list[ ,1:2, ])
}

if (stage_output == "all"){
  return(output.N.list[ , 1:3, ])
}
if (stage_output == "3"){
  return(output.N.list[ , 3, ])
}

if (stage_output == "size"){
  return(sizelist)
  
}

}


# out <- HYOSmodel(flow.data = discharge, temp.data = temp, disturbanceK = 40000, baselineK = 10000, Qmin = 0.25, extinct = 50, iteration = 1, peaklist = 0, peakeach = length(temp$Temperature))
# #------------------
# # Analyzing Results
# #-------------------
# # summarizing iterations
# means.list.HYOS <- mean.data.frame(out, burnin = 1, iteration = 1)
# means.list.HYOS <- cbind(means.list.HYOS[1:length(means.list.HYOS$mean.abund),], temps$dts[1:length(means.list.HYOS$mean.abund)])
# means.list.HYOS$`temps$dts` <- as.Date(means.list.HYOS$`temps$dts`)
# 
# 
# # means.list.HYOS <- as.data.frame(apply(out[, 3, ],MARGIN = 1, FUN = mean))
# # means.list.HYOS <- as.data.frame(cbind(means.list.HYOS[2:339,], temps$dts))
# # means.list.HYOS$`temps$dts` <- as.Date(means.list.HYOS$V2, origin = "1970-01-01")
# 
# # nt v nt+1
# plot(means.list.HYOS$mean.abund[300:600], means.list.HYOS$mean.abund[301:601], type = "b", xlab = "Nt", ylab= "Nt+1")
# # plot abundance over time
# 
# falls <- tibble(
#   x1 = c("2001-11-10", "2002-11-10", "2003-11-10", "2004-11-10", "2005-11-10", "2006-10-10"),
#   x2 = c("2001-11-10", "2002-11-10", "2003-11-10", "2004-11-10", "2005-11-10", "2006-10-10"),
#   y1 = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1),
#   y2 = c(0.23, 0.23, 0.23, 0.23, 0.23, 0.23)
#  )
# springs <- tibble(
#   x1 = c("2007-03-31", "2008-03-31", "2009-03-31", "2010-03-31", "2011-03-31", "2012-03-31"),
#   x2 = c("2007-03-31", "2008-03-31", "2009-03-31", "2010-03-31", "2011-03-31", "2012-03-31"),
#   y1 = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1),
#   y2 = c(0.23, 0.23, 0.23, 0.23, 0.23, 0.23)
# )
# # arrows <- tibble(
# #   x1 = c("2005-01-07", "2007-01-07", "2009-01-07", "2011-01-07"),
# #   x2 = c("2005-01-07", "2007-01-07", "2009-01-07", "2011-01-07"),
# #   y1 = c(1.45, 1.45, 1.45, 1.45), 
# #   y2 = c(1, 1, 1, 1)
# # )
# # 
# # arrows$x1 <- as.Date(arrows$x1)
# # arrows$x2 <- as.Date(arrows$x2)
# 
# falls$x1 <- as.Date(falls$x1)
# falls$x2 <- as.Date(falls$x2)
# springs$x1 <- as.Date(springs$x1)
# springs$x2 <- as.Date(springs$x2)
# 
# abund.trends.HYOS <- ggplot(data = means.list.HYOS[300:401,], aes(x =  `temps$dts`,
#                                                         y = mean.abund/10000, group = 1)) +
#   geom_point()+
#   # geom_ribbon(aes(ymin = mean.abund - 1.96 * se.abund,
#   #                ymax = mean.abund + 1.96 * se.abund),
#   #            colour = 'transparent',
#   #            alpha = .5,
#   #            show.legend = FALSE) +
#   geom_line(show.legend = FALSE) +
#   coord_cartesian(ylim = c(0,1.5)) +
#   ylab('Hydrospyche spp. Abundance/Reproductive Limit') +
#   xlab(" ")+
#   theme(text = element_text(size = 13), axis.text.x = element_text(angle=45, hjust = 1, size = 12.5), 
#         axis.text.y = element_text(size = 13))+
#   scale_x_date(date_labels="%B", date_breaks  ="6 months", limits = as.Date(c("2001-01-27", "2012-12-04"
#   )))+
#   
#   # annotate("segment", x = arrows$x1, y = arrows$y1, xend = arrows$x2, yend = arrows$y2,
#   #          arrow = arrow(type = "closed", length = unit(0.02, "npc")), color = "red")+
#   # annotate("text", x = arrows$x1[1], y = 1.5, label = "HI = 0.01", size = 4)+
#   # annotate("text", x = arrows$x1[2], y = 1.5, label = "HI = 0.1", size = 4)+
#   # annotate("text", x = arrows$x1[3], y = 1.5, label = "HI = 0.2", size = 4)+
#   # annotate("text", x = arrows$x1[4], y = 1.5, label = "HI = 0.5", size = 4 )
# 
#   annotate("segment", x = falls$x1, y = falls$y1, xend = falls$x2, yend = falls$y2,
#         arrow = arrow(type = "closed", length = unit(0.02, "npc")), color = "#a6611a")+
#   annotate("text", x = falls$x1[1], y = 0, label = "Fall HFE (0.45 Bankflow)" ,hjust = 0, size = 5, color = "#a6611a")+
#   annotate("segment", x = springs$x1, y = springs$y1, xend = springs$x2, yend = springs$y2,
#            arrow = arrow(type = "closed", length = unit(0.02, "npc")), color = "#018571")+
#   annotate("text", x = springs$x1[1], y = 0, label = "Spring HFE (0.45 Bankflow)" ,hjust = 0, size = 5, color = "#018571")
#   # 
#   ggsave(abund.trends.HYOS, filename = paste0("HYOSTempFlowHI", tempslist[te], peaklist[pe],".png"))
#   #ggsave(abund.trends.HYOS, filename = paste0("HYOSTempFlowHI", tempslist[te],".png"))
#   plotlist <- append(plotlist, paste0("HYOSTempFlowHI", tempslist[te], peaklist[pe],".png"))
# 
#   
# 
# plots <- lapply(ll <- plotlist ,function(x){
#   img <- as.raster(readPNG(x))
#   rasterGrob(img, interpolate = FALSE)
# })
#   ggsave(filename = paste0("HYOSTempFlowHI",peaklist[pe],".pdf"),width=8.5, height=11, 
#        marrangeGrob(grobs = plots, nrow = 3, ncol=2))
#   plotlist <- NULL
# 
# 
# 
# ggplot(data = NULL, mapping = aes(x = temps$dts, y = Total.N[2:2003]/10000))+
# geom_line(show.legend = FALSE) +
#   ylab('Hydrospyche spp. Abundance/Reproductive Limit') +
#   xlab(" ")
# #   for (te in 1:length(tempslist)){
# #   assign(paste0("p",te), readPNG(plotlist[te]))}
# #   grid.arrange(rasterGrob(p1), rasterGrob(p2), rasterGrob(p3), rasterGrob(p4), rasterGrob(p5), rasterGrob(p6), ncol = 3 )
# # 
# # # 
# # # 
# # pdf(paste0("HYOSTemp_", tempslist[te], "_Flood_HI_", peaklist[pe]), width = 8.27, height = 11.69)
# # gr <- grid.arrange(rasterGrob(p1), rasterGrob(p2), rasterGrob(p3), rasterGrob(p4), rasterGrob(p5), rasterGrob(p6), ncol = 3 )
# # 
# # library(patchwork) 
# # 
# # 
# # plotlist <- append(plotlist, paste0("HYOSTempFlowHI", tempslist[te], peaklist[pe],".png"))
# # 
# # 
# # ggsave(abund.trends.HYOS, filename = paste0("HYOSTempFlowHI", tempslist[te], peaklist[pe],".png"))
# # plotlist <- append(plotlist, paste0("HYOSTempFlowHI", tempslist[te], peaklist[pe],".png"))
# if (7 > temps$Temperature[t-1])  {
#   G1 <- 0.05/17
#   P1 <- 1-(1/17)
#   G2 <- 0 #no emergence
#   P2 <- 0.95# all remain in larval form
# }
# 
# if (temps$Temperature[t-1] > 17.5){
#   G1 <- 0.05
#   P1 <- 0
# }
# if (7 <= temps$Temperature[t-1] & temps$Temperature[t-1] <=17.5 & is.na(emergetime[t] == F)){
#   G1 <- 0.05/emergetime[t]
#   P1 <- 1-(1/emergetime[t])
# }
# if (7 <= temps$Temperature[t-1] & temps$Temperature[t-1] <= 17.5 & is.na(emergetime[t] == T)) {
#   G1 <- 0.05/((-0.95 * temps$Temperature[t-1]) + 24.75)
#   P1 <- 1-(1/((-0.95 * temps$Temperature[t-1]) + 24.75))}
