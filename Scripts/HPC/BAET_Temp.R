############################
# Temperature % increase Sensitivity BAET
#############################
# Load required libraries
library(parallel)  # For parallel processing
library(dataRetrieval)
# Load required functions and model
source("1spFunctions.R")
source("BAET_1sp_Model.R")


# read in Lees Ferry temp and discharge data from 2007 to 2023
temp <- readNWISdv("09380000", "00010", "2007-10-01", "2023-05-01")
discharge <- readNWISdv("09380000", "00060", "2007-10-01", "2023-05-01")

# calculate average yearly flows
flow <- average.yearly.flows(flowdata = discharge, flow.column_name = "X_00060_00003", date.column_name = "Date" )
# calculate average yearly temperatures
temps <- average.yearly.temp(tempdata = temp, temp.column_name = "X_00010_00003", date.column_name = "Date")
# create a timeseries of average temperatures 100 years long
temps <- rep.avg.year(temps, n = 100, change.in.temp = 0, years.at.temp = 0)
# create a timeseries of average flows 100 years long
flows <- do.call("rbind", replicate(100, flow, simplify = FALSE))
# match dates
flows$dts <- as.Date(temps$dts)
# get discharge magnitude by dividing by bankfull discharge 
flows$Discharge <- flows$Discharge/85000

# create sequence of temperature modifiers
temp_seq <- c(1, 1.1, 1.2, 1.5)

# makes some vectors for data to go into
BAET_temp_abund <- data.frame(abundance=numeric(), temperature=factor(), taxa=factor())
BAET_temp_biomass <- data.frame(biomass=numeric(), temperature = factor(), taxa = factor())
BAET_temp_percapita <- data.frame(percapita=numeric(), temperature = factor(), taxa = factor())

results <- mclapply(temp_seq, function(te) {
  set.seed(123) # make reproducible
  # model sizes
  temps$Temperature <- temps$Temperature * te
  
  sizes <- BAETmodel(flow.data = flows$Discharge, temp.data = temps, baselineK = 5000, disturbanceK = 40000 , Qmin = 0.1, extinct = 50, iteration = 1000, peaklist = 0, peakeach = length(temps$Temperature), stage_output = "size")
  # model abundances
  out <- BAETmodel(flow.data = flows$Discharge, temp.data = temps, baselineK = 5000, disturbanceK = 40000 , Qmin = 0.1, extinct = 50, iteration = 1000, peaklist = 0, peakeach = length(temps$Temperature))
  temps$Temperature <- temps$Temperature / te
  
  # for each stage, calculate mean biomass
  s1s <- colMeans(out[-c(1:260), 1, ]) * (0.0053 * (mean(sizes[-c(1:260)])/2)^2.875)
  s2s <- colMeans(out[-c(1:260), 2,]) * (0.0053 * (mean(sizes[-c(1:260)]))^2.875)
  s3s <- colMeans(out[-c(1:260), 3,]) * (0.0053 * (mean(sizes[-c(1:260)]))^2.875)
  # sum the mean biomass of each stage to get mean timestep biomass
  sizes_list <- as.vector(s1s + s2s + s3s)
  # Store biomass data in a dataframe
  average_size <- cbind(sizes_list, rep(te, times = length(sizes_list)), rep("BAET",times = length(sizes_list)))
  colnames(average_size) <- colnames(BAET_temp_biomass)
  # put the percapita biomass of each stage in a dataframe
  #s1pc <- (0.0053 * (mean(sizes[-c(1:260)])/2)^2.875)
  #s2pc <-(0.0053 * (mean(sizes[-c(1:260)]))^2.875)
  s3pc <-  (0.0053 * (mean(sizes[-c(1:260)]))^2.875)
  #percapita <- c(s1pc, s2pc, s3pc)
  percapita <- c(s3pc)
  
  # store percapita biomass data in a dataframe
  percapita_biomass <- cbind(percapita, rep(te, times = length(percapita)), rep("BAET",times = length(percapita)))
  colnames(percapita_biomass) <- colnames(BAET_temp_percapita)
  
  
  # calculate mean abundances at each timestep
  means.list.BAET <- mean.data.frame(out, burnin = 260, iteration = 1000)
  means <- means.list.BAET$mean.abund
  # Store abundance data in a dataframe
  average_means <- cbind(means, rep(te, times = length(means)), rep("BAET",times = length(means)))
  colnames(average_means) <- colnames(BAET_temp_abund)
  # Append results to the main data storage
  BAET_temp_abund <- rbind(BAET_temp_abund, average_means)
  BAET_temp_biomass <- rbind(BAET_temp_biomass, average_size)
  BAET_temp_percapita <- rbind(BAET_temp_percapita, percapita_biomass)
  
  # Return results as a list
  return(list(BAET_temp_abund = average_means, BAET_temp_biomass = average_size, BAET_temp_percapita = percapita_biomass))
}, mc.cores = detectCores() - 1)

# Combine results from all temperature scenarios into final dataframes
BAET_temp_abund <- do.call(rbind, lapply(results, `[[`, "BAET_temp_abund"))
BAET_temp_biomass <- do.call(rbind, lapply(results, `[[`, "BAET_temp_biomass"))
BAET_temp_percapita <- do.call(rbind, lapply(results, `[[`, "BAET_temp_percapita"))

# Write results to CSV files
write.csv(BAET_temp_abund, "BAET_temp_abund.csv", row.names = FALSE)
write.csv(BAET_temp_biomass, "BAET_temp_biomass.csv", row.names = FALSE)
write.csv(BAET_temp_percapita, "BAET_temp_percapita.csv", row.names = FALSE)


# tabula rasa
rm(temp)
rm(temps)
rm(BAET_temp_abund)
rm(BAET_temp_biomass)

## Now we add the summer spike to temperatures 
# read in LF temp and discharge data from 2007 to 2023
temp <- readNWISdv("09380000", "00010", "2007-10-01", "2023-05-01")
# calculate average yearly temperatures
temps <- average.yearly.temp(tempdata = temp, temp.column_name = "X_00010_00003", date.column_name = "Date")
# create summertime spike (up to 21 C, then scale from there)
temps$Temperature[16:23] <- c(14, 16, 18, 21, 21, 18, 16, 14)

# create a timeseries of average temperatures 100 years long
temps <- rep.avg.year(temps, n = 100, change.in.temp = 0, years.at.temp = 0)

# add temperature spike
# makes some vectors for data to go into
BAET_temp_abund_spike <- data.frame(abundance=numeric(), temperature=factor(), taxa=factor())
BAET_temp_biomass_spike <- data.frame(biomass=numeric(), temperature = factor(), taxa = factor())
BAET_temp_percapita_spike <- data.frame(percapita=numeric(), temperature = factor(), taxa = factor())

results <- mclapply(temp_seq, function(te) {
  set.seed(123) # make reproducible
  # model sizes
  temps$Temperature <- temps$Temperature * te
  sizes <- BAETmodel(flow.data = flows$Discharge, temp.data = temps, baselineK = 5000, disturbanceK = 40000 , Qmin = 0.1, extinct = 50, iteration = 1000, peaklist = 0, peakeach = length(temps$Temperature), stage_output = "size")
  # model abundances
  out <- BAETmodel(flow.data = flows$Discharge, temp.data = temps, baselineK = 5000, disturbanceK = 40000 , Qmin = 0.1, extinct = 50, iteration = 1000, peaklist = 0, peakeach = length(temps$Temperature))
  temps$Temperature <- temps$Temperature / te
  
  # for each stage, calculate mean biomass
  s1s <- colMeans(out[-c(1:260), 1, ]) * (0.0053 * (mean(sizes[-c(1:260)])/2)^2.875)
  s2s <- colMeans(out[-c(1:260), 2,]) * (0.0053 * (mean(sizes[-c(1:260)]))^2.875)
  s3s <- colMeans(out[-c(1:260), 3,]) * (0.0053 * (mean(sizes[-c(1:260)]))^2.875)
  # sum the mean biomass of each stage to get mean timestep biomass
  sizes_list <- as.vector(s1s + s2s + s3s)
  # Store biomass data in a dataframe
  average_size <- cbind(sizes_list, rep(te, times = length(sizes_list)), rep("BAET",times = length(sizes_list)))
  colnames(average_size) <- colnames(BAET_temp_biomass_spike)
  # put the percapita biomass of each stage in a dataframe
  #s1pc <-  (0.0053 * (mean(sizes[-c(1:260)])/2)^2.875)
  #s2pc <- (0.0053 * (mean(sizes[-c(1:260)]))^2.875)
  s3pc <-(0.0053 * (mean(sizes[-c(1:260)]))^2.875)
  #percapita <- c(s1pc, s2pc, s3pc)
  percapita <- c(s3pc)
  
  
  # store percapita biomass data in a dataframe
  percapita_biomass <- cbind(percapita, rep(te, times = length(percapita)), rep("BAET",times = length(percapita)))
  colnames(percapita_biomass) <- colnames(BAET_temp_percapita_spike)
  
  
  # calculate mean abundances at each timestep
  means.list.BAET <- mean.data.frame(out, burnin = 260, iteration = 1000)
  means <- means.list.BAET$mean.abund
  # Store abundance data in a dataframe
  average_means <- cbind(means, rep(te, times = length(means)), rep("BAET",times = length(means)))
  colnames(average_means) <- colnames(BAET_temp_abund_spike)
  # Append results to the main data storage
  BAET_temp_abund_spike <- rbind(BAET_temp_abund_spike, average_means)
  BAET_temp_biomass_spike <- rbind(BAET_temp_biomass_spike, average_size)
  BAET_temp_percapita_spike <- rbind(BAET_temp_percapita_spike, percapita_biomass)
  
  # Return results as a list
  return(list(BAET_temp_abund_spike = average_means, BAET_temp_biomass_spike = average_size, BAET_temp_percapita_spike = percapita_biomass))
}, mc.cores = detectCores() - 1)

# Combine results from all temperature scenarios into final dataframes
BAET_temp_abund_spike <- do.call(rbind, lapply(results, `[[`, "BAET_temp_abund_spike"))
BAET_temp_biomass_spike <- do.call(rbind, lapply(results, `[[`, "BAET_temp_biomass_spike"))
BAET_temp_percapita_spike <- do.call(rbind, lapply(results, `[[`, "BAET_temp_percapita_spike"))

# Write results to CSV files
write.csv(BAET_temp_abund_spike, "BAET_temp_abund_spike.csv", row.names = FALSE)
write.csv(BAET_temp_biomass_spike, "BAET_temp_biomass_spike.csv", row.names = FALSE)
write.csv(BAET_temp_percapita_spike, "BAET_temp_percapita_spike.csv", row.names = FALSE)

