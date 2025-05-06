##########################
# NZMS 1 sp model with rho of 0.55
###########################
library(purrr)
library(tidyverse)
library(lubridate)
library(plyr)
library(dplyr)
library(ggplot2)
library(readxl)
# data retrieval tool from USGS
library(dataRetrieval)

# Code for HPC - tidyverse has some issues on our HPC because one of the packages is deprecated
# We have to manually load all tidyverse packages
# library(purrr, lib.loc = "/home/ib/kurthena/R_libs/4Phosphorus-mediated changes in life history traits of the invasive New Zealand mudsnail (Potamopyrgus antipodarum) .2.1")
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
source("1spFunctions.R")
source("NZMS_shell_length_fecundity.R")
# #read in flow data from USGS gauge at Lees Ferry, AZ between 1985 to the end of the last water year
# temp <- readNWISdv("09380000", "00010", "2007-10-01", "2021-09-30")
# temps <- average.yearly.temp(temp, "X_00010_00003", "Date")
# 
# 
# n <- 77
# # qr is the temp ramps I want to increase the average Lees Ferry temp by 
# # how many years I want each temp ramp to last
# qr <- 0
# r <- 77
# 
# temps <- rep.avg.year(temps, n, change.in.temp = qr, years.at.temp = r)
# discharge <- rep(0.1, times = length(temps$Temperature))
# #temps$Temperature <- rep(12, times = length(temps$dts))

NZMSmodel <- function(flow.data, temp.data, baselineK, disturbanceK, Qmin, extinct, iteration, peaklist = NULL, peakeach = NULL){
  # source functions
  source("NZMSSurvivorship.R")
  Q <- as.numeric(flow.data)
  temps <- temp.data
  
  degreedays <- as.data.frame(cbind(temps$dts, temps$Temperature * 14))
  colnames(degreedays) <- c("dts", "DegreeDay")
  degreedays$dts <- as.Date(degreedays$dts, origin = "1970-01-01")
  
  # need to make ramped increasing hydropeaking index 
  hp <- c(rep(peaklist, each = peakeach))
  
  # specify iterations
  iterations <- iteration
  
  # baseline K in the absence of disturbance
  Kb <- as.numeric(baselineK)
  # max K after a big disturbance
  Kd <- as.numeric(disturbanceK)
  
  # specify baseline transition probabilities
  # its speculated (Cross et al 2010) that survivorship is between 80 - 100% for NZMS in Grand Canyon - will say 90% survive, 10% baseline mortality
  # from timestep to timestep, we expect 90% to survive so for stage 1 (which lasts aproximately 14 timesteps, survival should be approx 09^14 = 0.2287679
  # from that, only 1/14th will transition out, 13/14 remain in stage (Birt et al 2009)
  
  # stage 1 G1 (prob transition to stage 2)
  # stage 1 P1 (prob remaining in stage 2) 
  # stage 2 P2 (prob remaining in stage 2)
  # stage 3 P3 (prob remaining in stage 3) 
  
  G1 = 0.95/14
  G2 = 0.95/14
  P1 = 13/14
  P2 = 13/14
  P3 = 6/7 
  
  
  
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
  h <- surv.fit.NZMS$m$getPars()[2]  
  k <- surv.fit.NZMS$m$getPars()[1] 
  
  extinction <- extinct
  #-------------------------
  # Outer Loop of Iterations
  #--------------------------
  
  
  for (iter in c(1:iterations)) {
    
    K = Kb # need to reset K for each iteration
    # we can pull random values  
    #output.N.list[1,1:3, iter]<- rnbinom(3, size = 6020.175, mu = 35.11555 )/3
    output.N.list[1,1:3, iter] <- runif(3, min = 1, max = (0.3*K))
    # we often want to look at different parameter values after we run code, so we create some lists
    # list to input Ks
    Klist <- vector()
    Klist[1] <- 10000
    Flist <- vector()
    # list to imput flow morts
    flowmortlist <- vector()
    Flist <- vector()
    emergetime <- vector()
    sizelist <- vector()
    TempSurvival <- vector()
    for(c in temps$Temperature){
      
      b <- TempSurv_NZMS(c)
      
      TempSurvival <- append(TempSurvival, b)
    }
    
    
    #-------------------------
    # Inner Loop of Timesteps
    #-------------------------
    
    for (t in timestep) {
      
      #---------------------------------------------------------
      # Calculate starting fecundity per adult
      #temps
      # fecundities estimated from McKenzie et al. 2013 - reduced fecundity above 24 C and below 9 C. 
      # optimal temp between 16 and 19 C, but we don't really have parameterization for that
      # F2 <- -0.209*(temps$Temperature[t-1])^2 + 7.32*temps$Temperature[t-1] - 48.936
      # F3 <- -0.563*(temps$Temperature[t-1])^2 + 19.708*temps$Temperature[t-1] -131.764
      
      F2 <- 8.87473 * (-0.0001427 *(temps$Temperature[t-1] - 17.5)^4 + 1) *  hydropeaking.mortality(lower = 0, upper = 1, h = hp[t-1])
      
      F3 <-  27.89665 *(-0.0001427 * (temps$Temperature[t-1] - 17.5)^4 + 1) * hydropeaking.mortality(lower = 0, upper = 1, h = hp[t-1])
      # # #   
      if (F2 < 0){
        F2 <- 0
      }
      if (F3 < 0){
        F3 <- 0
      }
      #       F2 <- *
      #       F3 <- 
      # other fecundity options (0.35, 3.78, 14.6, 22.1)
      #
      #---------------------------------------------------
      # Calculate the disturbance magnitude-K relationship
      # Sets to 0 if below the Qmin
      Qf <- Qf.Function(Q[t-1], Qmin, a)
      
      #-------------------------------------------------------------------
      # Calculate K arrying capacity immediately following the disturbance
      K0 <- K + ((Kd-K)*Qf)
      
      # Calculate final K for timestep, including relationship between K and time since disturbance
      K <- post.dist.K(K0, Kb, g, t, Q, Qmin)
      Klist <- append(Klist, K)
      #---------------------------------------------
      # Calculate effect of density dependnce on fecundity 
      
      # Logistic Density Dependence on Fecundity via Rogosch et al. Fish Model
      F2 <- Logistic.Dens.Dependence(F2, K, Total.N[t-1, iter])  
      F3 <- Logistic.Dens.Dependence(F3, K, Total.N[t-1, iter]) 
      
      if (F2 < 0){
        F2 <- 0
      }
      if (F3 < 0){
        F3 <- 0
      }
      
      Flist <- append(Flist, (F2 + F3))
      
      # specify baseline transition probabilities
      # its speculated (Cross et al 2010) that survivorship is between 80 - 100% for NZMS in Grand Canyon - will say 80% survival to be on the conservative side
      # from timestep to timestep, we expect 90% to survive so for stage 1 (which lasts aproximately 14 timesteps, survival should be approx 09^14 = 0.2287679
      # from that, only 1/14th will transition out, 13/14 remain in stage (Birt et al 2009)
      
      # stage 1 G1 (prob transition to stage 2)
      # stage 1 P1 (prob remaining in stage 2) 
      # stage 2 P2 (prob remaining in stage 2)
      # stage 3 P3 (prob remaining in stage 3) 
      
      stageduration1 <- timestep_to_mat(temps$Temperature[t-1])[[1]]
      stageduration2 <- timestep_to_mat(temps$Temperature[t-1])[[2]]
      
      G1 = (0.91/stageduration1) *TempSurvival[t-1]
      G2 = (0.91/stageduration2)*TempSurvival[t-1]
      P1 = (1-(1/stageduration1)) *TempSurvival[t-1]
      P2 = (1-(1/stageduration2)) *TempSurvival[t-1]
      P3 = (1-(1/7))*TempSurvival[t-1]
      
      
      
      #-----------------------------------------------
      # Create Lefkovitch Matrix
      
      S1 <- c(P1, F2, F3)
      S2 <- c(G1, P2, 0)
      S3 <- c(0, G2, P3) 
      
      A <- rbind( S1, S2, S3)
      
      #-----------------------------------------------
      # Calculate new transition probabilities based on temperature
      # This is the growth v development tradeoff
      #--------------------------------------
      # Calculate abundances for each stage
      
      output.N.list[t, 1:3, iter] <- A %*% output.N.list[t-1, 1:3, iter] 
      
      #------------------------------------------
      #Calculate immediate mortality due to flows
      
      #s1
      output.N.list[t, 1, iter] <- flood.mortality(output.N.list[t,1,iter], k, h, Q[t-1], Qmin)
      #s2
      output.N.list[t,2,iter] <- flood.mortality(output.N.list[t,2,iter], k, h, Q[t-1], Qmin)
      #3
      output.N.list[t,3,iter] <- flood.mortality(output.N.list[t,3,iter], k, h, Q[t-1], Qmin)
      
      # in case we want to look back on what the flow mortality rates were
      flowmortlist <- append(flowmortlist, flood.mortality(1, k, h, Q[t-1], Qmin))
      
      #-------------------------------------------------
      # Calculate sum of all stages (total population)
      Total.N[,iter] <- apply(output.N.list[,,iter],1,sum)
      # check extinction threshold
      if (Total.N[t, iter] < extinction){
        output.N.list[t,,iter] <- 0
        Total.N[t, iter] <- 0
      } #-------------------------
      # End Inner Loop  
      #------------------------- 
    } #----------------------
    # End Outer Loop
    #----------------------
  }
  return(output.N.list)
  return(sizelist)
  #return(Flist)
  #return(Klist)
}
