############################
# Hydropeaking Sensitivity BAET
#############################
# Load necessary libraries
library(doParallel)
library(foreach)

# Load custom functions and data
source("1spFunctions.R")
source("BAET_1sp_Model.R")

# Read in Lees Ferry temp and discharge data from 2007 to 2023
temp <- readNWISdv("09380000", "00010", "2007-10-01", "2023-05-01")
discharge <- readNWISdv("09380000", "00060", "2007-10-01", "2023-05-01")

# Calculate average yearly flows
flow <- average.yearly.flows(flowdata = discharge, flow.column_name = "X_00060_00003", date.column_name = "Date")
# Calculate average yearly temperatures
temps <- average.yearly.temp(tempdata = temp, temp.column_name = "X_00010_00003", date.column_name = "Date")
# Create a timeseries of average temperatures 100 years long
temps <- rep.avg.year(temps, n = 100, change.in.temp = 0, years.at.temp = 0)
# Create a timeseries of average flows 100 years long
flows <- do.call("rbind", replicate(100, flow, simplify = FALSE))
# Match dates
flows$dts <- as.Date(temps$dts)
# Get discharge magnitude by dividing by bankfull discharge
flows$Discharge <- flows$Discharge / 85000

# Create sequence of hydropeaking intensities
hydropeak <- seq(0.00, 0.7, by = 0.05)

# Detect available cores for parallel execution
cores <- as.numeric(Sys.getenv("SLURM_CPUS_ON_NODE", 1)) # Default to 1 if not running under SLURM
cl <- makeCluster(cores)
registerDoParallel(cl)

# Parallel computation for hydropeaking scenarios
results <- foreach(hyd = 1:length(hydropeak), .combine = rbind) %dopar% {
  set.seed(123 + hyd) # Ensure reproducibility
  
  # Model sizes
  sizes <- BAETmodel(flow.data = flows$Discharge, temp.data = temps, baselineK = 5000, disturbanceK = 40000, Qmin = 0.1, extinct = 50, iteration = 1000, peaklist = hydropeak[hyd], peakeach = length(temps$Temperature), stage_output = "size")
  
  # Model abundances
  out <- BAETmodel(flow.data = flows$Discharge, temp.data = temps, baselineK = 5000, disturbanceK = 40000, Qmin = 0.1, extinct = 50, iteration = 1000, peaklist = hydropeak[hyd], peakeach = length(temps$Temperature))
  
  # For each stage, calculate mean biomass
  s1s <- mean(out[-c(1:260), 1, ]) * (0.0053 * (mean(sizes[-c(1:260)]) / 2)^2.875)
  s2s <- mean(out[-c(1:260), 2, ]) * (0.0053 * (mean(sizes[-c(1:260)]))^2.875)
  s3s <- mean(out[-c(1:260), 3, ]) * (0.0053 * (mean(sizes[-c(1:260)]))^2.875)
  
  # Sum the mean biomass of each stage to get mean timestep biomass
  sizemean <- sum(c(s1s, s2s, s3s))
  sizesd <- sd(c(s1s, s2s, s3s))
  
  # Calculate average annual biomass of stage 3 (emergent adults)
  s3ss <- as.data.frame(cbind(rowMeans(out[-c(1:260), 3, ]) * (0.0053 * ((sizes[-c(1:260)]))^2.875), year(temps$dts[-c(1:259)])))
  s3sYr <- aggregate(V1 ~ V2, data = s3ss, FUN = sum, na.rm = TRUE)
  S3Yrprod <- mean(s3sYr$V1, na.rm = TRUE)
  S3Yrprodsd <- sd(s3sYr$V1, na.rm = TRUE)
  # Calculate mean abundances at each timestep
  means.list.BAET <- mean.data.frame(out, burnin = 260, iteration = 1000)
  mean_abund <- mean(means.list.BAET$mean.abund)
  sd_abund <- sd(means.list.BAET$mean.abund, na.rm = TRUE)
  
  # Return the results for this hydropeaking scenario
  c(hydropeak[hyd], mean_abund, sd_abund, sizemean, sizesd, S3Yrprod,S3Yrprodsd)
}

# Stop the cluster
stopCluster(cl)

# Compile results into data frames
results_df <- as.data.frame(results)
colnames(results_df) <- c("Hydropeak", "MeanAbund", "SdAbund", "SizeMean", "SizeSd", "S3Yrprod", "S3Yrprodsd")
# Save the results as a CSV file
write.csv(results_df, file = "baet_hydropeaking_results.csv", row.names = FALSE)
