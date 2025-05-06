##########################
# CHIRis 1 sp model
###########################


library(purrr)
library(tidyverse)
library(lubridate)
library(plyr)
library(dplyr)
library(ggplot2)
# data retrieval tool from USGS
library(dataRetrieval)
CHIRmodel <- function(flow.data, temp.data, baselineK, disturbanceK, Qmin, extinct, iteration, peaklist = NULL, peakeach = NULL, stage_output = "all"){
  
  # set up model
  source("1spFunctions.R")
  source("CHIRSurvivorship.R")
  Q <- as.numeric(flow.data)
  temps <- temp.data
  
  degreedays <- as.data.frame(cbind(temps$dts, temps$Temperature * 14))
  colnames(degreedays) <- c("dts", "DegreeDay")
  degreedays$DegreeDay[degreedays$DegreeDay<0] <- 0
  degreedays$dts <- as.Date(degreedays$dts, origin = "1970-01-01")
  
  
  # need to make ramped increasing hydropeaking index 
  hp <- c(rep(peaklist, each = peakeach))
  
  # specify iterations
  iterations <- iteration
  
  # baseline K in the absence of disturbance
  Kb <- as.numeric(baselineK)
  # max K after a big disturbance
  Kd <- as.numeric(disturbanceK)
  
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
  h <- surv.fit.CHIR$m$getPars()[2]  
  k <- surv.fit.CHIR$m$getPars()[1] 
  
  extinction <- extinct
  
  
  #-------------------------
  # Outer Loop of Iterations
  #--------------------------
  
  
  for (iter in c(1:iterations)) {
    K = Kb # need to reset K for each iteration
    
    # pull random values from a uniform distribution 
    
    #output.N.list[1,1:3, iter] <- c(1166.201, 1041.528, 1417.288)
    
    output.N.list[1,1:3, iter] <- runif(3, min = 1, max = (0.3*K))
    
    # we often want to look at different parameter values after we run code, so we create some lists
    
    # list to input Ks
    Klist <- vector()
    Klist[1] <- 10000
    
    # list to imput flow morts
    flowmortlist <- vector()
    
    Flist <- vector()
    
    emergetime <- vector()
    
    sizelist <- vector()
    
    TempSurvival <- vector()
    for(c in temps$Temperature){
      b <- TempSurv_CHIR(c)
      TempSurvival <- append(TempSurvival, b)
    }
    #-------------------------
    # Inner Loop of Timesteps
    #-------------------------
    
    for (t in timestep) {
      
      #----------------------------------------------------------
      # Calculate how many timesteps emerging adults have matured
      t <- t
      emergetime <- append(emergetime, back.count.degreedays(t, 600, degreedays)) # value from Ali et al 1985
      #   if (t < 7){
      #     development <- append(development, delta[t-1]*(MaturationRate(devtime(temps$Temperature[t-1])/14)))
      #   } else {
      #   development <- append(development, development[t-delta[t]]*((MaturationRate(devtime(temps$Temperature[t-1])/14))/(MaturationRate(devtime(temps$Temperature[t-delta[t]])/14))))
      #   }
      # }
      #---------------------------------------------------------
      # Calculate fecundity per adult
      
      # we start by pulling fecundities from normal distribution
      # assuming 50 50 sex ration and 0.836 from Charles et al 2004
      F3 = 300 * 0.5 * hydropeaking.mortality(0.0, 0.6, h = hp[t-1])
      #CHIR egg # and % mortality from Charles et al 20040.836*
      # we can also relate fecundities to body size which is between 6 and 15 mm (also from Charles et al 2004)
      # we can "convert" emergetime to size by multiplying to get size between 6 and 15 mm and then convert to fecunity
      
      if (t > 19) {
        size <- 3*emergetime[t-1]-6
        sizelist <- append(sizelist, size)
        F3 <- (13.33*size)+180 * 0.5* hydropeaking.mortality(0.0, 0.6, h = hp[t-1])
      }
      # #--------------------------------------------------
      # Calculate the disturbance magnitude-K relationship
      # Sets to 0 if below the Qmin
      Qf <- Qf.Function(Q[t-1], Qmin, a)
      
      #-------------------------------------------------------------------
      # Calculate K carrying capacity immediately following the disturbance
      K0 <- K + ((Kd-K)*Qf)
      
      # Calculate final K for timestep, including relationship between K and time since disturbance
      K <- post.dist.K(K0, Kb, g, t, Q, Qmin)
      
      Klist <- append(Klist, K)
      #---------------------------------------------
      # Calculate effect of density dependence on fecundity 
      
      # Logistic via Rogosch et al. Fish Model
      # no immediate egg mortality incorporated
      F3 <- Logistic.Dens.Dependence(F3, K, Total.N[t-1, iter])
      # 
      # add F_CHIR to list
      Flist <- append(Flist, F3)
      #-----------------------------------------------
      # Calculate new transition probabilities based on temperature
      # This is the growth v development tradeoff
      # using Birt et al 2009 calcs
      # using survivals inspired by Thornton 1975 (lots of tests, survivals in the 70% range is reasonable since we don't know exact ph and salinity)
      
      # development measures
      # in this function, we assume that if below the min temp threshold (9) slow maturation
      # if above the max temp threshold (30), no one remains more than 1 timestep in each stage (fast maturation, small growth)
      
      if (5 > temps$Temperature[t-1]) {
        P1 <- (1-(1/20)) * TempSurvival[t-1]
        P2 <- P1
        G1 <- 0.74/20 * TempSurvival[t-1]
        G2 <- 0.66/20 * TempSurvival[t-1]
      }
      if (temps$Temperature[t-1] > 30){
        P1 <- 0
        P2 <- 0
        G1 <- 0.74 * TempSurvival[t-1]
        G2 <- 0.66 * TempSurvival[t-1]
      }
      
      if (5 <= temps$Temperature[t-1] & temps$Temperature[t-1] <= 30 & (is.na(emergetime[t-1]) == F)){
        G1 <- (0.74/((emergetime[t-1])/2)) * TempSurvival[t-1]
        G2 <- (0.66/((emergetime[t-1])/2)) * TempSurvival[t-1]
        P1 <- (1-(1/((emergetime[t-1])/2))) * TempSurvival[t-1]
        P2 <- P1
      }
      if (5 <= temps$Temperature[t-1] & temps$Temperature[t-1] <= 30 & (is.na(emergetime[t-1]) == T)) {
        G1 <- (0.74/((-0.136 * temps$Temperature[t-1]) + 5.088)) * TempSurvival[t-1]
        P1 <- (1-(1/((-0.136 * temps$Temperature[t-1]) + 5.088))) * TempSurvival[t-1]
        G2 <- (0.66/((-0.136 * temps$Temperature[t-1]) + 5.088)) * TempSurvival[t-1]
        P2 <- P1
      }
      #-----------------------------------------------
      # Create Lefkovitch Matrix
      
      T1 <- c(P1, 0, F3)
      T2 <- c(G1, P2, 0)
      T3 <- c(0, G2, 0) 
      
      A <- rbind( T1, T2, T3)
      
      #--------------------------------------
      # Calculate abundances for each stage
      
      output.N.list[t, 1:3, iter] <- A %*% output.N.list[t-1, 1:3, iter] 
      
      #------------------------------------------
      #Calculate immediate mortality due to flows
      # mortality due to flooding follows N0 = Nz*e^-hQ
      #s1
      output.N.list[t, 1, iter] <- flood.mortality(output.N.list[t, 1, iter], k, h, Q[t-1], Qmin)
      #s2Qt
      output.N.list[t,2,iter] <- flood.mortality(output.N.list[t,2,iter], k, h, Q[t-1], Qmin)
      
      #output.N.list[t,3,iter] <- flood.mortality(output.N.list[t,3,iter], k, h, Q[t-1], Qmin)
      
      flowmortlist <- append(flowmortlist, flood.mortality(1, k, h, Q[t-1], Qmin))
      
      #------------------------------------------------------
      # check extinction threshold and if below set to 0
      Total.N[t,iter] <- sum(output.N.list[t,,iter])
      if (Total.N[t,iter] < extinction){
        output.N.list[t,,iter] <- 0
        Total.N[t, iter] <- 0}
      
    } #-------------------------
    # End Inner Loop  
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
