############################
# Temperature % increase Sensitivity HYOS
#############################
# Load required libraries
library(parallel)  # For parallel processing
library(dataRetrieval)
# Load required functions and model
source("1spFunctions.R")
source("HYOS_1sp.R")


# read in LF temp and discharge data from 2007 to 2023
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
temp_seq <- c(1, 1.1,1.2, 1.5)

# makes some vectors for data to go into
HYOS_temp_abund <- data.frame(abundance=numeric(), temperature=factor(), taxa=factor())
HYOS_temp_biomass <- data.frame(biomass=numeric(), temperature = factor(), taxa = factor())
HYOS_temp_percapita <- data.frame(percapita=numeric(), temperature = factor(), taxa = factor())
HYOS_temp_propadult <- data.frame(propadult =numeric(), temperature = factor(), taxa = factor())

results <- mclapply(temp_seq, function(te) {
  set.seed(123) # make reproducible
  # model sizes
  temps$Temperature <- temps$Temperature * te
  
  # model sizes
  sizes <- HYOSmodel(flow.data = flows$Discharge, temp.data = temps, baselineK = 10000, disturbanceK = 40000 , Qmin = 0.3, extinct = 50, iteration = 1000, peaklist = 0, peakeach = length(temps$Temperature), stage_output = "size")
  # model abundances
  out <- HYOSmodel(flow.data = flows$Discharge, temp.data = temps, baselineK = 10000, disturbanceK = 40000 , Qmin = 0.3, extinct = 50, iteration = 1000, peaklist = 0, peakeach = length(temps$Temperature))
  
  temps$Temperature <- temps$Temperature / te
  
  # for each stage, calculate mean biomass
  s1s <- colMeans(out[-c(1:260), 1, ]) * (0.0046 * (mean(sizes[-c(1:260)]))^2.926)
  s2s <- colMeans(out[-c(1:260), 2, ]) * (0.0046 * (mean(sizes[-c(1:260)]+3))^2.926)
  s3s <- colMeans(out[-c(1:260), 3, ]) * (0.0046 * (mean(sizes[-c(1:260)]+3))^2.926)
  # sum the mean biomass of each stage to get mean timestep biomass
  sizes_list <- as.vector(s1s + s2s + s3s)
  # put the percapita biomass of each stage in a dataframe
  #s1pc <-(0.0046 * (mean(sizes[-c(1:260)]))^2.926)
  #s2pc <-(0.0046 * (mean(sizes[-c(1:260)]+3))^2.926)
  s3pc <- (0.0046 * (mean(sizes[-c(1:260)]+3))^2.926)
  #percapita <- c(s1pc, s2pc, s3pc)
  percapita <- c(s3pc)
  
  # calculate proportion of adults
  propadult <- (rowMeans(out[-c(1:260), 3, ])) / (rowMeans(out[-c(1:260), 1, ]) + rowMeans(out[-c(1:260), 2, ])+ rowMeans(out[-c(1:260), 3, ]))
  # store proportion adult biomass data in a dataframe
  propadult <- cbind(propadult, rep(te, times = length(propadult)), rep("HYOS",times = length(propadult)))
  colnames(propadult) <- colnames(HYOS_temp_propadult)
  
  
  # store percapita biomass data in a dataframe
  percapita_biomass <- cbind(percapita, rep(te, times = length(percapita)), rep("HYOS",times = length(percapita)))
  colnames(percapita_biomass) <- colnames(HYOS_temp_percapita)
  
  # Store biomass data in a dataframe
  average_size <- cbind(sizes_list, rep(te, times = length(sizes_list)), rep("HYOS",times = length(sizes_list)))
  colnames(average_size) <- colnames(HYOS_temp_biomass)
  # calculate mean abundances at each timestep
  means.list.HYOS <- mean.data.frame(out, burnin = 260, iteration = 1000)
  means <- (means.list.HYOS$mean.abund)
  # Store abundance data in a dataframe
  average_means <- cbind(means, rep(te, times = length(means)), rep("HYOS",times = length(means)))
  colnames(average_means) <- colnames(HYOS_temp_abund)
  # Append results to the main data storage
  HYOS_temp_abund <- rbind(HYOS_temp_abund, average_means)
  HYOS_temp_biomass <- rbind(HYOS_temp_biomass, average_size)
  HYOS_temp_percapita <- rbind(HYOS_temp_percapita, percapita_biomass)
  HYOS_temp_propadult <- rbind(HYOS_temp_propadult, propadult)
  
  # Return results as a list
  return(list(HYOS_temp_abund = average_means, HYOS_temp_biomass = average_size, HYOS_temp_percapita = percapita_biomass, HYOS_temp_propadult = propadult))
}, mc.cores = detectCores() - 1)

# Combine results from all temperature scenarios into final dataframes
HYOS_temp_abund <- do.call(rbind, lapply(results, `[[`, "HYOS_temp_abund"))
HYOS_temp_biomass <- do.call(rbind, lapply(results, `[[`, "HYOS_temp_biomass"))
HYOS_temp_percapita <- do.call(rbind, lapply(results, `[[`, "HYOS_temp_percapita"))
HYOS_temp_propadult <-  do.call(rbind, lapply(results, `[[`, "HYOS_temp_propadult"))

# Write results to CSV files
write.csv(HYOS_temp_abund, "HYOS_temp_abund.csv", row.names = FALSE)
write.csv(HYOS_temp_biomass, "HYOS_temp_biomass.csv", row.names = FALSE)
write.csv(HYOS_temp_percapita, "HYOS_temp_percapita.csv", row.names = FALSE)
write.csv(HYOS_temp_propadult, "HYOS_temp_propadult.csv", row.names = FALSE)


# tabula rasa
rm(temp)
rm(temps)
rm(HYOS_temp_abund)
rm(HYOS_temp_biomass)

## Now we add the summer spike to temperatures 
# read in LF temp and discharge data from 2007 to 2023
temp <- readNWISdv("09380000", "00010", "2007-10-01", "2023-05-01")
# calculate average yearly temperatures
temps <- average.yearly.temp(tempdata = temp, temp.column_name = "X_00010_00003", date.column_name = "Date")
# create summertime spike (up to 21 C, then scale from there)
temps$Temperature[16:23] <- c(14, 16, 18, 21, 21, 18, 16, 14)

# create a timeseries of average temperatures 100 years long
temps <- rep.avg.year(temps, n = 100, change.in.temp = 0, years.at.temp = 0)

# with temperature spike
# makes some vectors for data to go into
HYOS_temp_abund_spike <- data.frame(abundance=numeric(), temperature=factor(), taxa=factor())
HYOS_temp_biomass_spike <- data.frame(biomass=numeric(), temperature = factor(), taxa = factor())
HYOS_temp_percapita_spike <- data.frame(percapita=numeric(), temperature = factor(), taxa = factor())
HYOS_temp_propadult_spike <-  data.frame(propadult=numeric(), temperature = factor(), taxa = factor())

results <- mclapply(temp_seq, function(te) {
  set.seed(123) # make reproducible
  # model sizes
  temps$Temperature <- temps$Temperature * te
  
  # model sizes
  sizes <- HYOSmodel(flow.data = flows$Discharge, temp.data = temps, baselineK = 10000, disturbanceK = 40000 , Qmin = 0.3, extinct = 50, iteration = 1000, peaklist = 0, peakeach = length(temps$Temperature), stage_output = "size")
  # model abundances
  out <- HYOSmodel(flow.data = flows$Discharge, temp.data = temps, baselineK = 10000, disturbanceK = 40000 , Qmin = 0.3, extinct = 50, iteration = 1000, peaklist = 0, peakeach = length(temps$Temperature))
  
  temps$Temperature <- temps$Temperature / te
  # for each stage, calculate mean biomass
  s1s <- colMeans(out[-c(1:260), 1, ]) * (0.0046 * (mean(sizes[-c(1:260)]))^2.926)
  s2s <- colMeans(out[-c(1:260), 2, ]) * (0.0046 * (mean(sizes[-c(1:260)]+3))^2.926)
  s3s <- colMeans(out[-c(1:260), 3, ]) * (0.0046 * (mean(sizes[-c(1:260)]+3))^2.926)
  # sum the mean biomass of each stage to get mean timestep biomass
  sizes_list <- as.vector(s1s + s2s + s3s)
  # Store biomass data in a dataframe
  average_size <- cbind(sizes_list, rep(te, times = length(sizes_list)), rep("HYOS",times = length(sizes_list)))
  colnames(average_size) <- colnames(HYOS_temp_biomass_spike)
  
  # put the percapita biomass of each stage in a dataframe
  #s1pc <- (0.0046 * (mean(sizes[-c(1:260)]))^2.926)
  #s2pc <- (0.0046 * (mean(sizes[-c(1:260)]+3))^2.926)
  s3pc <- (0.0046 * (mean(sizes[-c(1:260)]+3))^2.926)
  #percapita <- c(s1pc, s2pc, s3pc)
  percapita <- c(s3pc)
  
  # store percapita biomass data in a dataframe
  percapita_biomass <- cbind(percapita, rep(te, times = length(percapita)), rep("HYOS",times = length(percapita)))
  colnames(percapita_biomass) <- colnames(HYOS_temp_percapita_spike)
  
  # calculate proportion of adults
  propadult <- (rowMeans(out[-c(1:260), 3, ])) / (rowMeans(out[-c(1:260), 1, ]) + rowMeans(out[-c(1:260), 2, ])+ rowMeans(out[-c(1:260), 3, ]))
  # store proportion adult biomass data in a dataframe
  propadult <- cbind(propadult, rep(te, times = length(propadult)), rep("HYOS",times = length(propadult)))
  colnames(propadult) <- colnames(HYOS_temp_propadult_spike)
  
  
  # calculate mean abundances at each timestep
  means.list.HYOS <- mean.data.frame(out, burnin = 260, iteration = 1000)
  means <- (means.list.HYOS$mean.abund)
  # Store abundance data in a dataframe
  average_means <- cbind(means, rep(te, times = length(means)), rep("HYOS",times = length(means)))
  colnames(average_means) <- colnames(HYOS_temp_abund_spike)
  # Append results to the main data storage
  HYOS_temp_abund_spike <- rbind(HYOS_temp_abund_spike, average_means)
  HYOS_temp_biomass_spike <- rbind(HYOS_temp_biomass_spike, average_size)
  HYOS_temp_percapita_spike <- rbind(HYOS_temp_percapita_spike, percapita_biomass)
  HYOS_temp_propadult_spike <- rbind(HYOS_temp_propadult_spike, propadult)
  
  # Return results as a list
  return(list(HYOS_temp_abund_spike = average_means, HYOS_temp_biomass_spike = average_size, HYOS_temp_percapita_spike = percapita_biomass, HYOS_temp_propadult_spike = propadult))
}, mc.cores = detectCores() - 1)

# Combine results from all temperature scenarios into final dataframes
HYOS_temp_abund_spike <- do.call(rbind, lapply(results, `[[`, "HYOS_temp_abund_spike"))
HYOS_temp_biomass_spike <- do.call(rbind, lapply(results, `[[`, "HYOS_temp_biomass_spike"))
HYOS_temp_percapita_spike <- do.call(rbind, lapply(results, `[[`, "HYOS_temp_percapita_spike"))
HYOS_temp_propadult_spike <-  do.call(rbind, lapply(results, `[[`, "HYOS_temp_propadult_spike"))

# Write results to CSV files
write.csv(HYOS_temp_abund_spike, "HYOS_temp_abund_spike.csv", row.names = FALSE)
write.csv(HYOS_temp_biomass_spike, "HYOS_temp_biomass_spike.csv", row.names = FALSE)
write.csv(HYOS_temp_percapita_spike, "HYOS_temp_percapita_spike.csv", row.names = FALSE)
write.csv(HYOS_temp_propadult_spike, "HYOS_temp_propadult_spike.csv", row.names = FALSE)
