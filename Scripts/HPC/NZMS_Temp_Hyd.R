############################
# Temperature % increase  NZMS
#############################

# Load required libraries
library(parallel)  # For parallel processing
library(dataRetrieval)
# Load required functions and model
source("1spFunctions.R")
source("NZMS_1sp_Model.R")


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
NZMS_temp_hyd_abund <- data.frame(abundance=numeric(), temperature=factor(), taxa=factor())
NZMS_temp_hyd_biomass <- data.frame(biomass=numeric(), temperature = factor(), taxa = factor())
results <- mclapply(temp_seq, function(te) {
  set.seed(123) # make reproducible
  # model sizes
  temps$Temperature <- temps$Temperature * te
  out <- NZMSmodel(flow.data = flows$Discharge, temp.data = temps, baselineK = 5000, disturbanceK = 9000 , Qmin = 0.3, extinct = 50, iteration = 1000, peaklist = 0.17, peakeach = length(temps$Temperature))
  temps$Temperature <- temps$Temperature / te
  # calculate mean abundances at each timestep
  means.list.NZMS <- mean.data.frame(out, burnin = 260, iteration = 1000)
  means <- means.list.NZMS$mean.abund
  # Store abundance data in a dataframe
  average_means <- cbind(means, rep(te, times = length(means)), rep("NZMS",times = length(means)))
  colnames(average_means) <- colnames(NZMS_temp_hyd_abund)
  # for each stage, calculate mean biomass
  s1s <- colMeans(out[-c(1:260), 1, ]) * (0.02 * mean(c(0.5, 3.2))^2.4315)
  s2s <- colMeans(out[-c(1:260), 2, ]) * (0.02 * mean(c(3.2, 4))^2.4315)
  s3s <- colMeans(out[-c(1:260), 3, ]) * (0.02 * mean(c(4, 5.5))^2.4315)
  # sum the mean biomass of each stage to get mean timestep biomass
  sizes_list <- as.vector(s1s + s2s + s3s)
  # Store biomass data in a dataframe
  average_size <- cbind(sizes_list, rep(te, times = length(sizes_list)), rep("NZMS",times = length(sizes_list)))
  colnames(average_size) <- colnames(NZMS_temp_hyd_biomass)

  # Append results to the main data storage
  NZMS_temp_hyd_abund <- rbind(NZMS_temp_hyd_abund, average_means)
  NZMS_temp_hyd_biomass <- rbind(NZMS_temp_hyd_biomass, average_size)
  # Return results as a list
  return(list(NZMS_temp_hyd_abund = average_means, NZMS_temp_hyd_biomass = average_size))
}, mc.cores = detectCores() - 1)

# Combine results from all temperature scenarios into final dataframes
NZMS_temp_hyd_abund <- do.call(rbind, lapply(results, `[[`, "NZMS_temp_hyd_abund"))
NZMS_temp_hyd_biomass <- do.call(rbind, lapply(results, `[[`, "NZMS_temp_hyd_biomass"))

# Write results to CSV files
write.csv(NZMS_temp_hyd_abund, "NZMS_temp_hyd_abund.csv", row.names = FALSE)
write.csv(NZMS_temp_hyd_biomass, "NZMS_temp_hyd_biomass.csv", row.names = FALSE)

# tabula rasa
rm(temp)
rm(temps)
rm(NZMS_temp_hyd_abund)
rm(NZMS_temp_hyd_biomass)

## Now we add the summer spike to temperatures 
# read in LF temp and discharge data from 2007 to 2023
temp <- readNWISdv("09380000", "00010", "2007-10-01", "2023-05-01")
# calculate average yearly temperatures
temps <- average.yearly.temp(tempdata = temp, temp.column_name = "X_00010_00003", date.column_name = "Date")
# create summertime spike (up to 21 C, then scale from there)
temps$Temperature[16:23] <- c(14, 16, 18, 21, 21, 18, 16, 14)

# create a timeseries of average temperatures 100 years long
temps <- rep.avg.year(temps, n = 100, change.in.temp = 0, years.at.temp = 0)


# makes some vectors for data to go into
NZMS_temp_hyd_abund_spike <- data.frame(abundance=numeric(), temperature=factor(), taxa=factor())
NZMS_temp_hyd_biomass_spike <- data.frame(biomass=numeric(), temperature = factor(), taxa = factor())
results <- mclapply(temp_seq, function(te) {
  set.seed(123) # make reproducible
  # modify temperature regim
  temps$Temperature <- temps$Temperature * te
  # add temp spike in September
  out <- NZMSmodel(flow.data = flows$Discharge, temp.data = temps, baselineK = 5000, disturbanceK = 9000 , Qmin = 0.3, extinct = 50, iteration = 1000, peaklist = 0.17, peakeach = length(temps$Temperature))
  temps$Temperature <- temps$Temperature / te
  # calculate mean abundances at each timestep
  means.list.NZMS <- mean.data.frame(out, burnin = 260, iteration = 1000)
  means <- means.list.NZMS$mean.abund
  # Store abundance data in a dataframe
  average_means <- cbind(means, rep(te, times = length(means)), rep("NZMS",times = length(means)))
  colnames(average_means) <- colnames(NZMS_temp_hyd_abund_spike)
  # for each stage, calculate mean biomass
  s1s <- colMeans(out[-c(1:260), 1, ]) * (0.02 * mean(c(0.5, 3.2))^2.4315)
  s2s <- colMeans(out[-c(1:260), 2, ]) * (0.02 * mean(c(3.2, 4))^2.4315)
  s3s <- colMeans(out[-c(1:260), 3, ]) * (0.02 * mean(c(4, 5.5))^2.4315)
  # sum the mean biomass of each stage to get mean timestep biomass
  sizes_list <- as.vector(s1s + s2s + s3s)
  # Store biomass data in a dataframe
  average_size <- cbind(sizes_list, rep(te, times = length(sizes_list)), rep("NZMS",times = length(sizes_list)))
  colnames(average_size) <- colnames(NZMS_temp_hyd_biomass_spike)
  
  # Append results to the main data storage
  NZMS_temp_hyd_abund_spike <- rbind(NZMS_temp_hyd_abund_spike, average_means)
  NZMS_temp_hyd_biomass_spike <- rbind(NZMS_temp_hyd_biomass_spike, average_size)
  # Return results as a list
  return(list(NZMS_temp_hyd_abund_spike = average_means, NZMS_temp_hyd_biomass_spike = average_size))
}, mc.cores =  detectCores() -1)

# Combine results from all temperature scenarios into final dataframes
NZMS_temp_hyd_abund_spike <- do.call(rbind, lapply(results, `[[`, "NZMS_temp_hyd_abund_spike"))
NZMS_temp_hyd_biomass_spike <- do.call(rbind, lapply(results, `[[`, "NZMS_temp_hyd_biomass_spike"))

# Write results to CSV files
write.csv(NZMS_temp_hyd_abund_spike, "NZMS_temp_hyd_abund_spike.csv", row.names = FALSE)
write.csv(NZMS_temp_hyd_biomass_spike, "NZMS_temp_hyd_biomass_spike.csv", row.names = FALSE)

