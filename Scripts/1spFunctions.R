#############################
## 1 species Model Functions
#############################
library(minpack.lm)

## checkpos makes sure that the K-occupied term is positive, assigns 0 if not
checkpos <- function(x) {
  ifelse(x < 0, 0, x)
}

## function to calculate daily shell size, to help determine stage transitions and fecundity at different stages
shell.growth <- function(m, b, start.size){
  l <- seq(1, 1000, by = 1)
  lengths <- vector(length = 1000)
  lengths[1] <- start.size
  for (i in l){
    lengths[i + 1] <- m*lengths[i]+b
  }
  return(lengths)
}

# function to index flow data, summarize as mean discharge per timestep, and relativize to flow magnitude (aka disturbance magnitude)
TimestepDischarge <- function(flow, bankfull_discharge){
  # Make an index to be used for aggregating
  ID <- as.numeric(as.factor(flow$Date))-1
  # want it to be every 14 days, hence the 14
  ID <- ID %/% 14
  # aggregate over ID and TYPE for all numeric data.
  out <- aggregate(flow[sapply(flow,is.numeric)],
                   by=list(ID), #,flow$X_00060_00003
                   FUN=mean)
  # format output
  names(out)[1:2] <-c("dts","Discharge")
  # add the correct dates as the beginning of every period
  out$dts <- format(as.POSIXct(flow$Date[((out$dts*14) + 2)]), "%Y-%m-%d")
  # get mean Discharge data for every 14 days
  #out <- aggregate(out, by = list(out$dts), FUN = mean)
  out$Discharge <- out$Discharge/bankfull_discharge # standardize to disturbance magnitude by taking discharge/bankfull_discharge
  # order by date in chronological order
  out <- out[order(out$dts),]
  return(out)
}

# function to index and summarize temperature data over timesteps length
TimestepTemperature <- function(temp){
  # Make an index to be used for aggregating
  ID <- as.numeric(as.factor(temp$Date))-1
  # want it to be every 14 days, hence the 14
  ID <- ID %/% 14
  # aggregate over ID and TYPE for all numeric data.
  outs <- aggregate(temp[sapply(temp,is.numeric)],
                    by=list(ID),
                    FUN=mean, na.rm=TRUE)
  # format output
  names(outs)[1:2] <-c("dts","Temperature")
  # add the correct dates as the beginning of every period
  outs$dts <- as.POSIXct(temp$Date[((outs$dts*14)+1)], origin = "1970-01-01")
  # order by date in chronological order
  temps <- outs[order(outs$dts),]
  
  return(temps)
}

#function to calculate degree days accumulated every timestep
TimestepDegreeDay <- function(temp, river){
  # Make an index to be used for aggregating
  ID <- as.numeric(as.factor(temp$Date))-1
  # want it to be every 14 days, hence the 14
  ID <- ID %/% 14
  # aggregate over ID and T
  degreeday <- aggregate(temp[sapply(temp,is.numeric)],
                         by=list(ID),
                         FUN=sum)
  names(degreeday)[1:2] <-c("dts","DegreeDay")
  # add the correct dates as the beginning of every period
  degreeday$dts <- as.POSIXct(temp$Date[(degreeday$dts*14)+1])
  # order by date in chronological order
  degreeday <- degreeday[order(degreeday$dts),]
  # can't have negative numbers so turn those into 0s
  degreeday$DegreeDay[which(degreeday$DegreeDay < 0)] <-  0
  
  if (river == "Colorado River"){
    # there are less temperatures than discharge readings, so we will just repeat the same temp sequence for this exercise
    degreeday <- degreeday[1:363,] # last row isn't a full timestep
    DDs <- rbind(degreeday, degreeday, degreeday) 
    DDs <- DDs[1:length(flow.magnitude$Discharge), ]
  } 
  return(DDs)
}


# forward lookign degree day calculation 
forward.count.degreedays <- function(criticaldegreedays){
  emergetime <- vector()
  for(degree in 1:length(degreedays$DegreeDay)){
    degseq <- seq(degree, length(degreedays$DegreeDay), by = 1)
    # create an empty vector to put the total number of degree days accumulated
    vec <- 0
    # for each value in that sequence, we will add the degree day values of 
    #the timestep prior and check if it adds up to our threshold to emergence
    for (s in degseq) {
      if(vec <= criticaldegreedays) {vec <- degreedays$DegreeDay[s] + vec
      emerg <- NA}
      else {emerg <- s - degree
      break}
      # once we hit that threshold, we count the number of timesteps it took to reach that and add that to our emergetime vector
    }
    emergetime[degree] <- emerg
  }
  
  return(emergetime)
}




# backwards looking calculation of timesteps required to reach stage 3
back.count.degreedays <- function(time.now, criticaldegreedays, degreedays){
  # for each timestep, we want to back calculate the number of degree days
  if(time.now == 2) {
    emerg <- NA
  } else {
    # create a sequence of time from last t to current t
    degseq <- as.vector(seq(time.now - 1, 1, by = -1))
    # create an empty vector to put the total number of degree days accumulated
    vec <- 0
    # for each value in that sequence, we will add the degree day values of 
    #the timestep prior and check if it adds up to our threshold to emergence
    for (s in degseq) {
      if(vec <= criticaldegreedays) {
        vec <- as.numeric(degreedays$DegreeDay[s]) + vec
        emerg <- NA
      }
      else {emerg <- time.now - s
      break
      }
      # once we hit that threshold, we count the number of timesteps it took to reach that and add that to our emergetime vector
    }
  }
  return(emerg)
}



# function to calculate Qf from McMullen et al 2017. Sets to 0 if below the Qmin
Qf.Function <- function(Q, Qmin, a){
  if (Q < Qmin) {
    Qf <- 0
  } else {
    Qf <- (Q - Qmin)/(a + Q- Qmin)
  }
  return(Qf)
}


# Function to calc. K as a function of time post-disturbance at a particular disturbance intensity
post.dist.K <- function(K0, Kb, g, t, Q, Qmin){
  #calculate tau (times since most recent distubance)
  taut = (t-1) - max(which(Q[1:(t-1)] > Qmin))
  if (is_empty(taut)) { taut <-  0
  K <- K0} 
  if (is.na(taut)==T | taut == 0) { taut <-  0
  K <- K0} 
  if (is.infinite(taut)==T){ taut <- 0
  K <- K0}
  if (taut > 0) {
    K <- Kb + ((K0 - Kb)*exp(-g*taut))} # function from McMullen et al 2017, g is shape function
  return(K)
}


# Function to calculate logistic density dependence on fecundity, after Rogosch et al 2019
Logistic.Dens.Dependence <- function(Fecundity, K, N){
  f.rate <- Fecundity * checkpos((K - N)/K) 
  return(f.rate)
}


#Ricker model (after Recruitment = axe^-bx, see Bolker Ch 3 Deterministic Functions for
#Ecological Modeling)
Ricker.Dens.Dependence <- function(b, N, fecundity){
  f.rate <- fecundity * exp(-b * N)
  return(f.rate)}
# b = 0.005
#F_NZMS <- Ricker.Dens.Dependence(b, Total.N[t-1, iter], F_NZMS) 

# beverton holt is Nt+1 = rNt/1-Nt(r-1)/K
# it is supposed to be depensatory, so as t -> inf, Nt+1 -> K, BUT 
# the discrete nature of this causes it overshoot by a lot, 
# meaning it isn't any better or worse than traditional logistric growth models
Bev.Holt.Dens.Dependence <- function(r, N, K, fecundity){
  if(N < K){
    f.rate <- fecundity * (K - N/K)
  } else {
    f.rate <- fecundity * (1/K)
  }
}

# F_NZMS <- Bev.Holt.Dens.Dependence(r, Total.N[t-1, iter], K, F_NZMS)

# growth.development tradeoff function 
growth.development.tradeoff <- function(temp, thresholdtemp.min, thresholdtemp.max, min.rate, max.rate){  # m and b from y = mx+b 
  if (thresholdtemp.min > temp) rate <- min.rate
  if (temp > thresholdtemp.max) rate <- max.rate
  xs <- c(thresholdtemp.min, thresholdtemp.max)
  ys <- c(min.rate, max.rate)
  lm <- lm(ys~xs)
  if (thresholdtemp.min <= temp & temp <= thresholdtemp.max) rate <- (summary(lm)$coefficients[2,1] * temp) + summary(lm)$coefficients[1,1]
  return(rate)
}


# mortality due to flooding follows N0 = Nz*e^-h
flood.mortality <- function(N, k, h, Q, Qmin){
  if (Q <= Qmin){
    newN <- N
  } else {
    newN <- N * k * exp(-h * Q)
  }
  if (newN > N){ #only mortality or no effect 
    newN <- N
  }
  return(newN)
}

#function to  code into mean population abundance over iterations
mean.data.frame <- function(data, burnin, iteration){
  repdf <- plyr::adply(data, c(1,2,3))
  #repdf <- plyr::adply(data, c(1,2))
  names(repdf) <- c('timesteps', 'stage', 'rep', 'abund')
  repdf$timesteps <- as.numeric(as.character(repdf$timesteps))
  
  # totn <- plyr::adply(Total.N, c(1,2))
  # names(totmesteps <- as.numeric(as.character(totn$timesteps))n) <- c('timesteps', 'rep', 'tot.abund')
  # totn$ti
  
  # joining totn and repdf together
  # repdf <- dplyr::left_join(totn, repdf)
  
  ## calculating relative abundance
  # repdf <- mutate(repdf, rel.abund = abund/tot.abund)
  repdf$timesteps <- as.factor(repdf$timesteps)
  ## Taking mean results to cf w/ observed data'
  
  means.list<- repdf %>%
    # select(-tot.abund) %>%
    dplyr::group_by(timesteps, rep) %>% # combining stages
    dplyr::summarise(abund = sum(abund)) %>%
    ungroup() %>%
    dplyr::group_by(timesteps) %>%
    dplyr::summarise(mean.abund = mean(abund),
                     sd.abund = sd(abund),
                     se.abund = sd(abund)/sqrt(iteration)) %>%
    ungroup()
  
  if (is.null(burnin)== F){
    means.list <- means.list[burnin:length(means.list$timesteps), ]
  }
  return(means.list)
}


# want to make a function to find the yearly average time series of a set of flows
average.yearly.flows <- function(flowdata, flow.column_name, date.column_name){
  require(tidyverse)
  flowdata$Date <- as_datetime(flowdata[[date.column_name]])
  flowdata$Date <- yday(flowdata$Date)
  flowdata$Discharge <- flowdata[[flow.column_name]]
  flowdata <- flowdata %>% group_by(Date) %>% dplyr::summarise(Discharge = mean(Discharge), 
                                                               sd = sd(Discharge), 
                                                               count = n(), 
                                                               se = sd(Discharge/count))
  # Make an index to be used for aggregating
  ID <- as.numeric(as.factor(flowdata$Date)) 
  # want it to be every 14 days, hence the 14
  ID <- ID %/% 14
  
  ID[which(ID ==26)] <- 25
  # aggregate over ID and TYPEall numeric data.
  outs <- aggregate(flowdata[sapply(flowdata,is.numeric)],
                    by=list(ID),
                    FUN=mean)
  names(outs)[2:3] <-c("dts","Discharge")
  # add the correct dates as the beginning of every period
  outs$dts <- strptime(round(outs$dts), "%j") ###Note need to subtract 365 if larger than 365
  
  # order by date in chronological order#
  #outs <- outs[order(outs$dts),]
  outs$dts <- as_date(outs$dts)
  return(outs)
}

# want to make a function to find the yearly average time series of a set of temps
average.yearly.temp <- function(tempdata, temp.column_name, date.column_name){
  require(tidyverse)
  tempdata$Date <- as_datetime(tempdata[[date.column_name]])
  tempdata$Date <- yday(tempdata$Date)
  tempdata$Temperature <- tempdata[[temp.column_name]]
  tempdata <- tempdata %>% group_by(Date) %>% dplyr::summarise(Temperature = mean(Temperature), 
                                                               sd = sd(Temperature), 
                                                               count = n(), 
                                                               se = sd(Temperature/count))
  # Make an index to be used for aggregating
  ID <- as.numeric(as.factor(tempdata$Date)) 
  # want it to be every 14 days, hence the 14
  ID <- ID %/% 14
  
  ID[which(ID ==26)] <- 25
  # aggregate over ID and TYPEall numeric data.
  outs <- aggregate(tempdata[sapply(tempdata,is.numeric)],
                    by=list(ID),
                    FUN=mean)
  names(outs)[2:3] <-c("dts","Temperature")
  # add the correct dates as the beginning of every period
  outs$dts <- strptime(round(outs$dts), "%j") ###Note need to subtract 365 if larger than 365
  
  # order by date in chronological order#
  #outs <- outs[order(outs$dts),]
  outs$dts <- as_date(outs$dts)
  return(outs)
}


# function that creates a repeating yearly cyclical timeseries for n years, also allows temperatures to be manupulated
rep.avg.year <- function(data, n, change.in.temp = 0, years.at.temp = 0){
  yr <- seq(from = 2000, to = 2000+(n-1), by = 1)
  temp_seq <- do.call("rbind", replicate(n, data, simplify = FALSE))
  temp_seq$dts <- as.Date(temp_seq$dts)
  year(temp_seq$dts) <- rep(yr, each = 26) # want to make different years for each rep of the timestep
  #temp_seq$Temperature <- temp_seq$Temperature + rep(change.in.temp, years.at.temp*26)
  return(temp_seq)
}

hydrofunction <- function(x){exp(-x*2)}

hydropeaking.mortality <- function(lower, upper, hp){
  int <- integrate(hydrofunction, lower = lower, upper = upper) 
  int$value <- int$value + 0.6
  if (lower < hp & hp < upper) {
    hydro <- integrate(hydrofunction, lower = lower, upper = hp) 
    int$value <- int$value - hydro$value
  }
  if (hp >= upper){
    int$value <- 0
  }
  #H <- (1-int$value)*(1-hp)
  return(int$value)
}