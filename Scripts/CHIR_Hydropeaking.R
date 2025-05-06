############################
# Hydropeaking Sensitivity CHIR
#############################
# Load necessary libraries
library(doParallel)
library(foreach)

# Load custom functions and data
source("1spFunctions.R")
source("CHIR_1sp_Model.R")


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

# create sequence of hydropeaking intensities
hydropeak <- seq(0.0, 0.7, by = 0.05)
# Detect available cores for parallel execution
cores <- as.numeric(Sys.getenv("SLURM_CPUS_ON_NODE", 1)) # Default to 1 if not running under SLURM
cl <- makeCluster(cores)
registerDoParallel(cl)

# Parallel computation for hydropeaking scenarios
results <- foreach(hyd = 1:length(hydropeak), .combine = rbind) %dopar% {
  set.seed(123 + hyd) # Ensure reproducibility
  # model sizes
  sizes <- CHIRmodel(flow.data = flows$Discharge, temp.data = temps, baselineK = 10000, disturbanceK = 40000 , Qmin = 0.2, extinct = 50, iteration = 2, peaklist = hydropeak[hyd], peakeach = length(temps$Temperature), stage_output = "size")
  # model abundances
  out <- CHIRmodel(flow.data = flows$Discharge, temp.data = temps, baselineK = 10000, disturbanceK = 40000 , Qmin = 0.2, extinct = 50, iteration = 2, peaklist = hydropeak[hyd], peakeach = length(temps$Temperature))
  s1s <- mean(out[-c(1:260), 1, ]) * (0.0018 * (mean((sizes[-c(1:260)])/2))^2.617)
  s2s <- mean(out[-c(1:260), 2,]) * (0.0018 * (mean(sizes[-c(1:260)]))^2.617)
  s3s <- mean(out[-c(1:260), 3,]) * (0.0018 * (mean(sizes[-c(1:260)]))^2.617)
  # sum the mean biomass of each stage to get mean timestep biomass
  sizemean<- sum(c(s1s, s2s, s3s))
  # produce standard devation
  sizesd<- sd(c(s1s, s2s, s3s))
  
  
  # now instead of getting the mean, calculate biomass at every timestep
  s1ss <- as.data.frame(cbind(rowMeans(out[-c(1:260), 1, ]) * (0.0018 * ((sizes[-c(1:260)]/2))^2.617), year(temps$dts[-c(1:259)])))
  s2ss <- as.data.frame(cbind(rowMeans(out[-c(1:260), 2, ]) * (0.0018 * ((sizes[-c(1:260)]))^2.617), year(temps$dts[-c(1:259)])))
  s3ss <- as.data.frame(cbind(rowMeans(out[-c(1:260), 3, ]) * (0.0018 * ((sizes[-c(1:260)]))^2.617), year(temps$dts[-c(1:259)])))
  
  
  # find annual stage specific biomass based on year
  s1sYr <- aggregate(V1 ~ V2, data = s1ss, FUN = sum, na.rm = TRUE)
  s2sYr <- aggregate(V1 ~ V2, data = s2ss, FUN = sum, na.rm = TRUE)
  s3sYr <- aggregate(V1 ~ V2, data = s3ss, FUN = sum, na.rm = TRUE)
  
  # add all stages together to get average annual biomass (aka secondary production)
  #Yrprod[hyd] <- sum(mean(c(s1sYr$V1, s2sYr$V1, s3sYr$V1), na.rm = T))
  
  # average annual biomass of only stage 3 (emergent adults)
  S3Yrprod <- mean(s3sYr$V1, na.rm = T)
  S3Yrprodsd <- sd(s3sYr$V1, na.rm = T)
  # calculate mean abundances at each timestep
  means.list.CHIR <- mean.data.frame(out, burnin = 260, iteration = 2)
  # calculate the average of mean abundances at each hydropeaking intensity
  mean_abund <- mean(means.list.CHIR$mean.abund)
  # calculate the standard deviation of mean abundances at each hydropeaking intensity
  sd_abund  <- sd(means.list.CHIR$mean.abund, na.rm = T)
  
  c(hydropeak[hyd], mean_abund, sd_abund, sizemean, sizesd, S3Yrprod, S3Yrprodsd)
  
}
#Stop the cluster
stopCluster(cl)

# Compile results into data frames
results_df <- as.data.frame(results)
colnames(results_df) <- c("Hydropeak", "MeanAbund", "SdAbund", "SizeMean", "SizeSd", "S3Yrprod", "S3Yrprodsd")
# Save the results as a CSV file
write.csv(results_df, file = "chir_hydropeaking_results.csv", row.names = FALSE)

# compile abundance data
# CHIR_hyd_means <- as.data.frame(cbind(hydropeak, means, sd, rep("CHIR", length(means))))

# compile timeseries biomass data
# CHIR_hyd_size <- as.data.frame(cbind(hydropeak, sizemeans, sizesd, rep("CHIR", length(sizemeans))))

# compile annual production biomass data
# CHIR_hyd_yrprod <- as.data.frame(cbind(hydropeak, S3Yrprod, rep("CHIR", length(sizemeans))))
# 
# ggplot(data = CHIR_hyd_means, aes(x = hydropeak,  y = means, group = 1, color = "CHIR")) +
#   geom_ribbon(aes(ymin = means - sd,
#                   ymax = means + sd),
#               colour = 'transparent',
#               alpha = .15,
#               show.legend = T) +
#   geom_point(show.legend = T, linewidth = 1, alpha = 0.8) +
#   #geom_line(data = flow.magnitude, aes(x = as.Date(dts), y = X_00060_00003), color = "blue") +
#   #geom_line(data = temps, aes(x = as.Date(dts), y = Temperature*1000), color = "green")+
#   #coord_cartesian(ylim = c(0,6000)) +S1 and S2 (inds/m2)
#   ylab('CHIR spp. Abund.') +
#   xlab("")+
#   labs(colour=" ")+
#   theme_bw()+
#   scale_color_manual(values = colors)+
#   #scale_y_continuous(
#   # sec.axis = sec_axis(~., name="CHIRidae Larvae (inds/m2)"
#   # ))+
#   theme(text = element_text(size = 13), axis.text.x = element_text(angle=45, hjust = 1, size = 12.5), 
#         axis.text.y = element_text(size = 13), )


# ggplot(data = CHIR_hyd_size, aes(x = hydropeak,  y = sizemeans, group = 1, color = "CHIR")) +
#   geom_ribbon(aes(ymin = sizemeans - sizesd,
#                   ymax = sizemeans + sizesd),
#               colour = 'transparent',
#               alpha = .15,
#               show.legend = T) +
#   geom_line(show.legend = T, linewidth = 1, alpha = 0.8) +
#   #geom_line(data = flow.magnitude, aes(x = as.Date(dts), y = X_00060_00003), color = "blue") +
#   #geom_line(data = temps, aes(x = as.Date(dts), y = Temperature*1000), color = "green")+
#   #coord_cartesian(ylim = c(0,6000)) +S1 and S2 (inds/m2)
#   ylab('CHIR spp. Biomass (mg)') +
#   xlab("")+
#   labs(colour=" ")+
#   theme_bw()+
#   scale_color_manual(values = colors)+
#   #scale_y_continuous(
#   # sec.axis = sec_axis(~., name="CHIRidae Larvae (inds/m2)"
#   # ))+
#   theme(text = element_text(size = 13), axis.text.x = element_text(angle=45, hjust = 1, size = 12.5), 
#         axis.text.y = element_text(size = 13), )
# 
