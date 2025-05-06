##############################
#Figure Generation for Manuscript 2
#############################
library(ggpubr)
library(tidyverse)
library(patchwork)
library(svglite)
library(dataRetrieval)
library(cowplot)
library(png)
library(grid)
library(lm.beta)
library(effectsize)
library(flextable)
library(officer)
# Set global options for digits
options(digits = 3)

#############################
# Temperature plot
source("1spFunctions.R")
temp <- readNWISdv("09380000", "00010", "2015-10-01", "2022-05-01")

temperature <-ggplot(temp, aes(x = as.Date(Date), y = X_00010_00003)) +
  geom_line(aes(color = "Colorado \n River at \n Lee's \n Ferry"), linewidth = 1) +  # Correct legend mapping
  labs(x = " ", y = "Temperature", title = "Water Temperature", color = " ") +  # Legend title
  theme_bw() +
  scale_color_manual(values = c("Colorado \n River at \n Lee's \n Ferry" = "black")) +  # Assign legend color
  scale_x_date(date_labels = "%Y") +
  theme(
    text = element_text(size = 15),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 15), 
    axis.text.y = element_text(size = 15)
  )

  
#################################
# Comparision of model output and empirical data
source("BAETValidation.R")
source("CHIR_Validation.R")
source("GAMM_Validation.R")
source("HYOSValidation.R")
source("NZMS_Validation.R")
Fig1 <- ggarrange(NZMSts, BAETts, GAMMts,temperature, 
          labels = c("a", "b", "c", "d"),
          ncol = 2, nrow = 2, common.legend =F)
ggsave(filename = "fig3.1.png", Fig1, device = "png", dpi = "retina", height = 9, width = 11)


# # hydropeaking intensity 
source("BAET_Hydropeaking.R")
source("HYOS_Hydropeaking.R")
source("NZMS_Hydropeaking.R")
source("GAMM_Hydropeaking.R")
source("CHIR_Hydropeaking.R")

hyos_hydropeaking <- read_csv("hyos_hydropeaking_results.csv")
hyos_hydropeaking$Taxa <- rep("HYOS", times = length(hyos_hydropeaking$Hydropeak))
chir_hydropeaking <- read_csv("chir_hydropeaking_results.csv")
chir_hydropeaking$Taxa <- rep("CHIR", times = length(chir_hydropeaking$Hydropeak))
baet_hydropeaking_results <- read_csv("baet_hydropeaking_results.csv")
baet_hydropeaking_results$Taxa <- rep("BAET", times = length(baet_hydropeaking_results$Hydropeak))
nzms_hydropeaking <- read_csv("nzms_hydropeaking_results.csv")
nzms_hydropeaking$Taxa <- rep("NZMS", times = length(nzms_hydropeaking$Hydropeak))
gamm_hydropeaking <- read_csv("gamm_hydropeaking_results.csv")
gamm_hydropeaking$Taxa <- rep("GAMM", times = length(gamm_hydropeaking$Hydropeak))
#Hydro_Abund <- as.data.frame(rbind(BAET_hyd_means, HYOS_hyd_means, NZMS_hyd_means, CHIR_hyd_means, GAMM_hyd_means))
Hydro_Abund <- as.data.frame(rbind(baet_hydropeaking_results, hyos_hydropeaking, nzms_hydropeaking, chir_hydropeaking, gamm_hydropeaking))
# combine all abundance data
Hydro_Abund$Hydropeak <- as.numeric(Hydro_Abund$Hydropeak)

# make sure in the correct format

taxa_offsets <- c("Baetidae" = 0,
                  "Chironomidae" = -0.01,
                  "G. lacustris" = 0,
                  "Hydropsyche" = 0.001,
                  "P. antipodarum" = 0)

# Add offset column to both datasets
Hydro_Abund <- Hydro_Abund %>%
  mutate(Hydropeak_offset = Hydropeak + taxa_offsets[(Taxa)])

Hydro_Abund$MeanAbund <- as.numeric(Hydro_Abund$MeanAbund)
Hydro_Abund$SdAbund <- as.numeric(Hydro_Abund$SdAbund)
Hydro_Abund$Taxa <- as.factor(Hydro_Abund$Taxa)

first_zero <- first_zero %>%
  mutate(Hydropeak_offset = Hydropeak + taxa_offsets[(Taxa)])

hyd_abund <- ggplot(data = Hydro_Abund, aes(x = Hydropeak_offset, y = log(MeanAbund), group = Taxa, color = Taxa)) +
  geom_line(linewidth = 1, alpha = 0.8) +
  geom_ribbon(aes(ymin = log(MeanAbund - SdAbund),ymax = log(MeanAbund + SdAbund),
                  fill = Taxa), alpha = .15, color = "transparent", show.legend = FALSE) +
  geom_point(data = first_zero,
             aes(x = Hydropeak_offset, y = 0.25, color = Taxa), shape = 8, size = 3,
             stroke = 2, show.legend = FALSE) +
  geom_vline(aes(xintercept = 0.17), linetype = "dashed", linewidth = 1) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4), limits = c(0, 11), expand = c(0, 0)) +
  scale_color_manual(
    name = "Taxon",
    labels = c(expression(paste(italic("Baetidae")," spp.")),
               expression(paste(italic("Chironomidae"), " spp.")),
               expression(italic("G. lacustris")),
               expression(paste(italic("Hydropsyche"), " spp.")),
               expression(italic("P. antipodarum"))),
    values = c("#66CCEE", "#228833", "#CCBB44", "#AA3377", "#4477AA")) +
  scale_fill_manual(
    name = "Taxon",
    labels = c("Baetidae spp.", "Chironomidae spp.", "G. lacustris", "Hydropsyche spp.", "P. antipodarum"),
    values = c("#66CCEE", "#228833", "#CCBB44", "#AA3377", "#4477AA")) +
  theme_bw() +
  xlab("Hydropeaking Intensity") +
  ylab("Log Abundance per Reach") +
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5),
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))


# read in mean timestep biomass data
#Hydro_TS_Biomass <- rbind(BAET_hyd_size, HYOS_hyd_size, NZMS_hyd_size, CHIR_hyd_size, GAMM_hyd_size)
# make sure in correct format
# Hydro_TS_Biomass$hydropeak <- as.numeric(Hydro_TS_Biomass$hydropeak)
# Hydro_TS_Biomass$sizemeans <- as.numeric(Hydro_TS_Biomass$sizemeans)
# Hydro_TS_Biomass$sizesd <- as.numeric(Hydro_TS_Biomass$sizesd)
# Hydro_TS_Biomass$V4 <- as.factor(Hydro_TS_Biomass$V4)
#plot
hyd_ts <- ggplot(data = Hydro_Abund, aes(x = Hydropeak_offset, y = log(SizeMean), group = Taxa, color = Taxa)) +
  geom_line(linewidth = 1, alpha = 0.8) +
  geom_ribbon(aes(ymin = log(SizeMean - SizeSd), ymax = log(SizeMean + SizeSd),
                  fill = Taxa), alpha = .15, color = "transparent", show.legend = FALSE) +
  geom_point(data = first_zero,
             aes(x = Hydropeak_offset, y = 0.25, color = Taxa), shape = 8, size = 3,
             stroke = 2, show.legend = FALSE) +
  geom_vline(aes(xintercept = 0.17), linetype = "dashed", linewidth = 1) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4), limits = c(0, 12), expand = c(0, 0)) +
  scale_color_manual(
    name = "Taxon",
    labels = c(expression(paste(italic("Baetidae")," spp.")),
               expression(paste(italic("Chironomidae"), " spp.")),
               expression(italic("G. lacustris")),
               expression(paste(italic("Hydropsyche"), " spp.")),
               expression(italic("P. antipodarum"))),
    values = c("#66CCEE", "#228833", "#CCBB44", "#AA3377", "#4477AA")) +
  scale_fill_manual(
    name = "Taxon",
    labels = c("Baetidae spp.", "Chironomidae spp.", "G. lacustris", "Hydropsyche spp.", "P. antipodarum"),
    values = c("#66CCEE", "#228833", "#CCBB44", "#AA3377", "#4477AA")) +
  theme_bw() +
  xlab("Hydropeaking Intensity") +
  ylab("Log Standing Biomass per Reach") +
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5),
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))


# read in mean annual productivity/biomass data
#Hydro_Yr_Prod <- as.data.frame(rbind(BAET_hyd_yrprod, HYOS_hyd_yrprod, CHIR_hyd_yrprod, NZMS_hyd_yrprod, GAMM_hyd_yrprod))
# make sure in correct format

# Hydro_Yr_Prod$hydropeak <- as.numeric(Hydro_Yr_Prod$hydropeak)
# Hydro_Yr_Prod$S3Yrprod <- as.numeric(Hydro_Yr_Prod$S3Yrprod)
# Hydro_Yr_Prod$V3 <- as.factor(Hydro_Yr_Prod$V3)

Hyd_yr <- ggplot(data = Hydro_Abund, aes(x = Hydropeak_offset, y = log(S3Yrprod ), group = Taxa, color = Taxa)) +
  geom_line(linewidth = 1, alpha = 0.8) +
  geom_ribbon(aes(ymin = log(S3Yrprod  - S3Yrprodsd), ymax = log(S3Yrprod  + S3Yrprodsd),
                  fill = Taxa), alpha = .15, color = "transparent", show.legend = FALSE) +
  geom_point(data = first_zero,
             aes(x = Hydropeak_offset, y = 0.25, color = Taxa), shape = 8, size = 3,
             stroke = 2, show.legend = FALSE) +
  geom_vline(aes(xintercept = 0.17), linetype = "dashed", linewidth = 1) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4), limits = c(0, 11), expand = c(0, 0)) +
  scale_color_manual(
    name = "Taxon",
    labels = c(expression(paste(italic("Baetidae")," spp.")),
               expression(paste(italic("Chironomidae"), " spp.")),
               expression(italic("G. lacustris")),
               expression(paste(italic("Hydropsyche"), " spp.")),
               expression(italic("P. antipodarum"))),
    values = c("#66CCEE", "#228833", "#CCBB44", "#AA3377", "#4477AA")) +
  scale_fill_manual(
    name = "Taxon",
    labels = c("Baetidae spp.", "Chironomidae spp.", "G. lacustris", "Hydropsyche spp.", "P. antipodarum"),
    values = c("#66CCEE", "#228833", "#CCBB44", "#AA3377", "#4477AA")) +
  theme_bw() +
  xlab("Hydropeaking Intensity") +
  ylab("Log Emergent Standing Biomass (mg)") +
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5),
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

Fig3.8 <- ggarrange(hyd_abund, hyd_ts, Hyd_yr,
          labels = c("a", "b", "c"),
          ncol = 3, nrow = 1, common.legend = T)
ggsave("Fig3.8.png", plot = Fig3.8, device = "png", dpi = "retina", width = 8, height = 6, units = c("in"))

##################################################
# Plots for different temperature scenarios

# read in Lees Ferry temp and discharge data from 2007 to 2023
temp <- readNWISdv("09380000", "00010", "2007-10-01", "2023-05-01")
temps <- average.yearly.temp(tempdata = temp, temp.column_name = "X_00010_00003", date.column_name = "Date")
# create a timeseries of average temperatures 100 years long
temps <- rep.avg.year(temps, n = 100, change.in.temp = 0, years.at.temp = 0)


# Define the temperature multipliers and labels
scenarios <- data.frame(
  multiplier = c(1, 1.1, 1.2, 1.5),
  label = c("base", "+10%", "+20%", "+50%")
)

# Apply modifications to create the different temperature regimes
temperature.regime <- map2_dfr(scenarios$multiplier, scenarios$label, ~ {
  temps[1:25, ] %>%
    mutate(Temperature = Temperature * .x, se = .y)
})

regime <- ggplot(data = temperature.regime, aes(x = as.Date(dts), y = Temperature, color = se))+
  geom_line(size = 1)+
  theme_bw()+
  scale_color_manual(name = "Warming Regime", labels=c("+10%", "+20%", "+50%", "Baseline"), values=c("#abd9e9", "#fdae61","#d7191c","#2c7bb6"))+
  scale_x_date(date_labels="%b", date_breaks  ="4 month")+
  xlab("")+
  theme(text = element_text(size = 11), axis.text.x = element_text(hjust = 1, size = 10, angle = 0), 
        axis.text.y = element_text(size = 10), legend.key = element_rect(fill = "transparent"))

hyos_abund <- read_csv("HYOS_temp_abund.csv")
baet_abund <- read_csv("BAET_temp_abund.csv")
chir_abund <- read_csv("CHIR_temp_abund.csv")
gamm_abund <- read_csv("GAMM_temp_abund.csv")
nzms_abund <- read_csv("NZMS_temp_abund.csv")


# with spike

# read in LF temp and discharge data from 2007 to 2023
temp <- readNWISdv("09380000", "00010", "2007-10-01", "2023-05-01")
# calculate average yearly temperatures
temps <- average.yearly.temp(tempdata = temp, temp.column_name = "X_00010_00003", date.column_name = "Date")
# create summertime spike (up to 21 C, then scale from there)
temps$Temperature[16:22] <- c(14, 16, 18, 21, 21, 18, 16, 14)

# create a timeseries of average temperatures 100 years long
temps <- rep.avg.year(temps, n = 100, change.in.temp = 0, years.at.temp = 0)

# Define the temperature multipliers and labels
scenarios <- data.frame(
  multiplier = c(1, 1.1, 1.2, 1.5),
  label = c("base", "+10%", "+20%", "+50%")
)

# Apply modifications to create the different temperature regimes
temperature.regime <- map2_dfr(scenarios$multiplier, scenarios$label, ~ {
  temps[1:25, ] %>%
    mutate(Temperature = Temperature * .x, se = .y)
})

#add temperature spike in September
spikeregime <- ggplot(data = temperature.regime, aes(x = as.Date(dts), y = Temperature, color = se))+
  geom_line(size = 1)+
  theme_bw()+
  scale_color_manual(name = "Warming Regime", labels=c("+10%", "+20%", "+50%", "Baseline"), values=c("#abd9e9", "#fdae61","#d7191c","#2c7bb6"))+
  scale_x_date(date_labels="%b", date_breaks  ="4 month")+
  xlab("")+
  theme(text = element_text(size = 11), axis.text.x = element_text(hjust = 1, size = 10, angle = 0), 
        axis.text.y = element_text(size = 10), legend.key = element_rect(fill = "transparent"))

############################
# Temperature Survival Curves
source("GAMMSurvivorship.R")
source("HYOSSurvivorship.R")
source("BAETSurvivorship.R")
source("NZMSSurvivorship.R")
source("CHIRSurvivorship.R")

gamm_surv <- ggplot(data = GAMMSurvRates, aes(x = Temperature, y = Survival, color = "Taxa"))+
  geom_point(aes(color = "data"), size = 1.5, alpha = 1, show.legend = T)+
  geom_line(data = tempGamm, aes(x = tem, y = V2), linewidth = 1, show.legend = T, alpha = 0.8)+
  theme_bw()+
  ylim(c(0,1))+
  ylab("G. lacustris Survival")+
  scale_color_manual(name = " ", labels=c("Empirical Data", "Fitted"), values=c("#CCBB44","black"))+
  theme(text = element_text(size = 13), axis.text.x = element_text(hjust = 1, angle=45, size = 12), 
        axis.text.y = element_text(size = 13),legend.position = "bottom", legend.key = element_rect(fill = "transparent"))

hyos_surv <- ggplot(data = HYOSSurvRates, aes(x = Temperature, y = Survival, color = "Taxa"))+
  geom_point(aes(color = "data"), size = 1.5, alpha = 1, show.legend = T)+
  geom_line(data = tempHyos, aes(x = tem, y = V2), linewidth = 1, show.legend = T, alpha = 0.8)+
  theme_bw()+
  ylab("Hydropsyche spp. Survival")+
  scale_color_manual(name = " ", labels=c("Empirical Data", "Fitted"), values=c("#AA3377","black"))+
  theme(text = element_text(size = 13), axis.text.x = element_text(hjust = 1, angle=45, size = 12), 
        axis.text.y = element_text(size = 13),legend.position = "bottom", legend.key = element_rect(fill = "transparent"))

baet_surv <- ggplot(data = BAETSurvRate, aes(x = Temperature, y = Survival, color = "Taxa"))+
  geom_point(aes(color = "data"), size = 1.5, alpha = 1, show.legend = T)+
  geom_line(data = tempBaet, aes(x = tem, y = V2), linewidth = 1, show.legend = T, alpha = 0.8)+
  theme_bw()+
  ylab("Baetidae spp. Survival")+
  scale_color_manual(name = " ", labels=c("Empirical Data", "Fitted"), values=c("#66CCEE","black"))+
  theme(text = element_text(size = 13), axis.text.x = element_text(hjust = 1, angle=45, size = 12), 
        axis.text.y = element_text(size = 13),legend.position = "bottom", legend.key = element_rect(fill = "transparent"))

nzms_surv <- ggplot(data = NZMSSurvRates, aes(x = Temperature, y = Survival, color = "Taxa"))+
  geom_point(aes(color = "data"), size = 1.5, alpha = 1, show.legend = T)+
  geom_line(data = tempHyos, aes(x = tem, y = V2), linewidth = 1, show.legend = T, alpha = 0.8)+
  theme_bw()+
  ylab("P. antipordarum Survival")+
  scale_color_manual(name = " ", labels=c("Empirical Data", "Fitted"), values=c("#4477AA","black"))+
  theme(text = element_text(size = 13), axis.text.x = element_text(hjust = 1, angle=45, size = 12), 
        axis.text.y = element_text(size = 13),legend.position = "bottom", legend.key = element_rect(fill = "transparent"))

chir_surv <- ggplot(data = CHIRSurvRate, aes(x = Temp, y = Survival, color = "Taxa"))+
  geom_point(aes(color = "data"), size = 1.5, alpha = 1, show.legend = T)+
  geom_line(data = tempHyos, aes(x = tem, y = V2), linewidth = 1, show.legend = T, alpha = 0.8)+
  theme_bw()+
  ylim(c(0,1))+
  ylab("Chironomidae spp. Survival")+
  xlab("Temperature")+
  scale_color_manual(name = " ", labels=c("Empirical Data", "Fitted"), values=c("#228833","black"))+
  theme(text = element_text(size = 13), axis.text.x = element_text(hjust = 1, angle=45, size = 12), 
        axis.text.y = element_text(size = 13),legend.position = "bottom", legend.key = element_rect(fill = "transparent"))

tempsurv_plot <- ggarrange(nzms_surv, baet_surv, gamm_surv, hyos_surv, chir_surv,  
          labels = c("a", "b", "c", "d", "e"),
          ncol = 3, nrow = 2, common.legend =F)

ggsave("TempSurvPlot.png", plot = tempsurv_plot, device = "png", width = 14, height = 12, dpi = "retina")


##############################################
# Hydropysche spp Scenario
hyos_abund <- read_csv("HYOS_temp_abund.csv")
hyos_pc <- read_csv("HYOS_temp_percapita.csv")
hyos_biomass <- read_csv("HYOS_temp_biomass.csv")
hyos_abund_sp <- read_csv("HYOS_temp_abund_spike.csv")
hyos_pc_sp <- read_csv("HYOS_temp_percapita_spike.csv")
hyos_biomass_sp <- read_csv("HYOS_temp_biomass_spike.csv")
hyos_abund_t.h <- read_csv("HYOS_temp_hyd_abund.csv")
hyos_biomass_t.h <- read_csv("HYOS_temp_hyd_biomass.csv")
hyos_abund_t.h_spike <- read_csv("HYOS_temp_hyd_abund_spike.csv")
hyos_biomass_t.h_spike <- read_csv("HYOS_temp_hyd_biomass_spike.csv")
hyos_abund_HFE <- read_csv("HYOS_temp_abund_HFE.csv")
hyos_biomass_HFE <- read_csv("HYOS_temp_biomass_HFE.csv")
hyos_biomass_HFE_sp <- read_csv("HYOS_temp_biomass_HFE_spike.csv")
hyos_abund_HFE_sp <- read_csv("HYOS_temp_abund_HFE_spike.csv")
hyos_abund_hyd_HFE <- read_csv("HYOS_temp_hyd_abund_HFE.csv")
hyos_biomass_hyd_HFE <- read_csv("HYOS_temp_hyd_biomass_HFE.csv")
hyos_abund_hyd_HFE.sp <- read_csv("HYOS_temp_hyd_abund_HFE_spike.csv")
hyos_biomass_hyd_HFE.sp <- read_csv("HYOS_temp_hyd_biomass_HFE_spike.csv")
hyos_pc_hfe <- read_csv("HYOS_temp_percapita_HFE.csv")
hyos_pc_hfe_sp <- read_csv("HYOS_temp_percapita_HFE_spike.csv")

# abundance
hyos_abund_combo <- bind_rows(
  hyos_abund %>% mutate(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 0),
  hyos_abund_sp %>% mutate(source = "Temperature & spike", temp_spike = 1, hydropeaking = 0, HFE = 0),
  hyos_abund_t.h %>% mutate(source = "Temperature & hydropeaking", temp_spike = 0, hydropeaking = 1, HFE = 0),
  hyos_abund_t.h_spike %>% mutate(source = "Temperature & spike & hydropeaking", temp_spike = 1, hydropeaking = 1, HFE = 0),
  hyos_abund_HFE %>% mutate(source = "Temperature & HFE",  temp_spike = 0, hydropeaking = 0, HFE = 1),
  hyos_abund_HFE_sp %>% mutate(source = "Temperature & spike & HFE",  temp_spike = 1, hydropeaking = 0, HFE = 1),
  hyos_abund_hyd_HFE %>% mutate(source = "Temperature & hydropeaking & HFE",  temp_spike = 0, hydropeaking = 1, HFE = 1),
  hyos_abund_hyd_HFE.sp %>% mutate(source = "Temperature & spike & hydropeaking & HFE",  temp_spike = 1, hydropeaking = 1, HFE = 1)
)

# biomass
hyos_biomass_combo <- bind_rows(
  hyos_biomass %>% mutate(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 0),
  hyos_biomass_sp %>% mutate(source = "Temperature & spike", temp_spike = 1, hydropeaking = 0, HFE = 0),
  hyos_biomass_t.h %>% mutate(source = "Temperature & hydropeaking", temp_spike = 0, hydropeaking = 1, HFE = 0),
  hyos_biomass_t.h_spike %>% mutate(source = "Temperature & spike & hydropeaking",  temp_spike = 1, hydropeaking = 1, HFE = 0),
  hyos_biomass_HFE %>% mutate(source = "Temperature & HFE",  temp_spike = 0, hydropeaking = 0, HFE = 1),
  hyos_biomass_HFE_sp %>% mutate(source = "Temperature & spike & HFE",  temp_spike =1, hydropeaking = 0, HFE = 1),
  hyos_biomass_hyd_HFE %>% mutate(source = "Temperature & hydropeaking & HFE",  temp_spike = 0, hydropeaking = 1, HFE = 1),
  hyos_biomass_hyd_HFE.sp %>% mutate(source = "Temperature & spike & hydropeaking & HFE",  temp_spike = 1, hydropeaking = 1, HFE = 1)
)

hyos_pc_combo <- bind_rows(
  hyos_pc %>% mutate(source = "Temperature"),
  hyos_pc_sp %>% mutate(source = "Temperature & spike"),
)


paired.colors <- c("#56B4E9", "#365C8DFF",
                   "#FDE725FF","#E69F00", 
                   "#9FDA3AFF","#009E73", 
                   "#CC79A7","#46337EFF")

label.names <- c(
  "Temperature", 
  "Temperature +\nSpike", 
  "Temperature +\nHydropeaking", 
  "Temperature +\nSpike + Hydropeaking",
  "Temperature +\nHFE", 
  "Temperature +\nSpike + HFE", 
  "Temperature +\nHydropeaking + HFE", 
  "Temperature +\nSpike + Hydropeaking + HFE"
)

hyos_abund_combo$abundance <- scale(hyos_abund_combo$abundance)
hyos_abund_combo$source <- ordered(hyos_abund_combo$source, levels = c(
  "Temperature", 
  "Temperature & spike", 
  "Temperature & hydropeaking", 
  "Temperature & spike & hydropeaking",
  "Temperature & HFE", 
  "Temperature & spike & HFE", 
  "Temperature & hydropeaking & HFE", 
  "Temperature & spike & hydropeaking & HFE"))  
hyos_means_df_abund <- hyos_abund_combo %>%
  group_by(source, temperature) %>%
  summarise(mean_abund = mean(abundance), .groups = "drop") %>%
  mutate(is_min = mean_abund == min(mean_abund)) %>%
  ungroup()


h.abund <- ggplot(data = hyos_abund_combo, aes(y = abundance, x = as.factor(temperature), color = source, group = source)) +
  stat_summary(fun.data = mean_sdl, geom = "errorbar", linewidth = 0.75, width = 0.2, 
               position = position_dodge(0.5), show.legend = F) +
  geom_point(data = hyos_means_df_abund, aes(x = as.factor(temperature), 
                                               y = mean_abund, shape = is_min, 
                                               stroke = 1, color = source, group = source), 
             size = 3, position = position_dodge(0.5), show.legend = FALSE) + # make extinct values a star
  scale_shape_manual(values = c("TRUE" = 8, "FALSE" = 16)) +  
  scale_x_discrete(labels=c("1" = "Baseline", "1.1" = "+10%",
                            "1.2" = "+20%", "1.5" = "+50%"))+
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1))+
  theme_bw()+
  scale_color_manual(name = " ", labels=label.names, values=paired.colors)+
  xlab("Temperature")+
  labs(y=expression(paste(italic("Hydropsyche"), "  Scaled Abundance")))+
  theme(text = element_text(size = 13), axis.text.x = element_text(hjust = 1, angle=45, size = 12), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))
# hyos biomass
hyos_biomass_combo$biomass <- scale(hyos_biomass_combo$biomass)
hyos_biomass_combo$source <- ordered(hyos_biomass_combo$source, levels = c(
  "Temperature", 
  "Temperature & spike", 
  "Temperature & hydropeaking", 
  "Temperature & spike & hydropeaking",
  "Temperature & HFE", 
  "Temperature & spike & HFE", 
  "Temperature & hydropeaking & HFE", 
  "Temperature & spike & hydropeaking & HFE"))  
library(dplyr)
library(ggplot2)

# group means
hyos_means_df_biomass <- hyos_biomass_combo %>%
  group_by(source, temperature) %>%
  summarise(mean_biomass = mean(biomass), .groups = "drop") %>%
  mutate(is_min = mean_biomass == min(mean_biomass)) %>%
  ungroup()

h.biomass <- ggplot(data = hyos_biomass_combo, aes(y = biomass, x = as.factor(temperature), color = source, group = source)) +
  stat_summary(fun.data = mean_sdl, geom = "errorbar", linewidth = 0.75, width = 0.2,
               position = position_dodge(0.5)) +
  stat_summary(fun = mean, geom = "line", size = 0.75, position = position_dodge(0.5)) +
  geom_point(data = hyos_means_df_biomass, aes(x = as.factor(temperature), 
                                               y = mean_biomass, shape = is_min, 
                                               stroke = 1, color = source, group = source), 
             size = 3, position = position_dodge(0.5), show.legend = FALSE) + # make extinct values a star
  scale_shape_manual(values = c("TRUE" = 8, "FALSE" = 16)) +  
  scale_x_discrete(labels = c("1" = "Baseline", "1.1" = "+10%", "1.2" = "+20%", "1.5" = "+50%")) +
  scale_color_manual(name = " ", labels = label.names, values = paired.colors) +
  theme_bw() +
  xlab("Temperature") +
  labs(y = expression(paste(italic("Hydropsyche"), " spp. Scaled Biomass"))) +
  theme(text = element_text(size = 13),
        axis.text.x = element_text(hjust = 1, angle = 45, size = 11),
        axis.text.y = element_text(size = 13),
        legend.key = element_rect(fill = "transparent"))

h_plot <- wrap_plots(h.abund, h.biomass) +
  plot_annotation(tag_levels = 'a')+
  plot_layout(guides = 'collect') &
  theme(
    legend.text = element_text(size = 11, lineheight = 0.9),
    legend.spacing.y = unit(0.5, "cm"),
    legend.key.height = unit(1.25, "cm")
  )
# Read the PNG image
caddis <- readPNG("caddisfly.png")
caddis_grob <- rasterGrob(caddis, x = 0.78, y = 0.92, width = 0.135, height = 0.135)  # Adjust x, y, width, and height

# Overlay the image on the plot
h.fig <- ggdraw() +
  draw_plot(h_plot) +
  draw_grob(caddis_grob)

ggsave("hfig.png", plot = h.fig, device = "png", width = 8.44, height = 6, units = c("in"), dpi = "retina")

h.pc <- ggplot(data = hyos_pc_combo, aes(x = as.factor(temperature), y = percapita, group = source, color = source))+
  stat_summary(fun.data = mean_sdl, geom = "errorbar", linewidth = 0.75, width = 0.2) +
  stat_summary(fun = mean, geom = "point", size = 3)+
  stat_summary(fun = mean, geom = "line", size = 0.75) +
  scale_x_discrete(labels=c("1" = "Baseline", "1.1" = "+10%",
                            "1.2" = "+20%", "1.5" = "+50%"))+
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1))+
  theme_bw()+
  scale_color_manual(name = " ", labels= c("Temperature", 
                                           "Temperature +\nSpike"), values= c("#56B4E9", "#365C8DFF"))+
  xlab("Temperature")+
  guides(color=guide_legend(title="Scenario"))+
  labs(y="Adult Biomass (mg)", title = expression(paste(italic("Hydropsyche"), " spp.")))+
  theme(text = element_text(size = 13), axis.text.x = element_text(hjust = 1, angle=45, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

###########################################
## Baetidae spp. Scenarios
baet_abund <- read_csv("BAET_temp_abund.csv")
baet_pc <- read_csv("BAET_temp_percapita.csv")
baet_biomass <- read_csv("BAET_temp_biomass.csv")
baet_abund_sp <- read_csv("BAET_temp_abund_spike.csv")
baet_pc_sp <- read_csv("BAET_temp_percapita_spike.csv")
baet_abund_t.h <- read_csv("BAET_temp_hyd_abund.csv")
baet_biomass_sp <- read_csv("BAET_temp_biomass_spike.csv")
baet_biomass_t.h <- read_csv("BAET_temp_hyd_biomass.csv")
baet_abund_t.h_spike <- read_csv("BAET_temp_hyd_abund_spike.csv")
baet_biomass_t.h_spike <- read_csv("BAET_temp_hyd_biomass_spike.csv")
baet_abund_HFE <- read_csv("BAET_temp_abund_HFE.csv")
baet_biomass_HFE <- read_csv("BAET_temp_biomass_HFE.csv")
baet_abund_HFE_sp <- read_csv("BAET_temp_abund_HFE_spike.csv")
baet_biomass_HFE.sp <- read_csv("BAET_temp_biomass_HFE_spike.csv")
baet_abund_HFE_sp <- read_csv("BAET_temp_abund_HFE_spike.csv")
baet_abund_hyd_HFE <- read_csv("BAET_temp_hyd_abund_HFE.csv")
baet_biomass_hyd_HFE <- read_csv("BAET_temp_hyd_biomass_HFE.csv")
baet_abund_hyd_HFE.sp <- read_csv("BAET_temp_hyd_abund_HFE_spike.csv")
baet_biomass_hyd_HFE.sp <- read_csv("BAET_temp_hyd_biomass_HFE_spike.csv")


# abundance
baet_abund_combo <- bind_rows(
  baet_abund %>% mutate(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 0),
  baet_abund_sp %>% mutate(source = "Temperature & spike", temp_spike = 1, hydropeaking = 0, HFE = 0),
  baet_abund_t.h %>% mutate(source = "Temperature & hydropeaking", temp_spike = 0, hydropeaking = 1, HFE = 0),
  baet_abund_t.h_spike %>% mutate(source = "Temperature & spike & hydropeaking", temp_spike =1, hydropeaking = 1, HFE = 0),
  baet_abund_HFE %>% mutate(source = "Temperature & HFE", temp_spike = 0, hydropeaking = 0, HFE = 1),
  baet_abund_HFE_sp %>% mutate(source = "Temperature & spike & HFE", temp_spike = 1, hydropeaking = 0, HFE = 1),
  baet_abund_hyd_HFE %>% mutate(source = "Temperature & hydropeaking & HFE", temp_spike = 0, hydropeaking = 1, HFE = 1),
  baet_abund_hyd_HFE.sp %>% mutate(source = "Temperature & spike & hydropeaking & HFE", temp_spike = 1, hydropeaking = 1, HFE = 1),
)

# biomass
baet_biomass_combo <- bind_rows(
  baet_biomass %>% mutate(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 0),
  baet_biomass_sp %>% mutate(source = "Temperature & spike", temp_spike = 1, hydropeaking = 0, HFE = 0),
  baet_biomass_t.h %>% mutate(source = "Temperature & hydropeaking", temp_spike = 0, hydropeaking = 1, HFE = 0),
  baet_biomass_t.h_spike %>%mutate(source = "Temperature & spike & hydropeaking", temp_spike =1, hydropeaking = 1, HFE = 0),
  baet_biomass_HFE %>%mutate(source = "Temperature & HFE", temp_spike = 0, hydropeaking = 0, HFE = 1),
  baet_biomass_HFE.sp %>%mutate(source = "Temperature & spike & HFE", temp_spike = 1, hydropeaking = 0, HFE = 1),
  baet_biomass_hyd_HFE %>% mutate(source = "Temperature & hydropeaking & HFE", temp_spike = 0, hydropeaking = 1, HFE = 1),
  baet_biomass_hyd_HFE.sp %>% mutate(source = "Temperature & spike & hydropeaking & HFE", temp_spike = 1, hydropeaking = 1, HFE = 1),
)

baet_pc_combo <- bind_rows(
  baet_pc %>% mutate(source = "Temperature"),
  baet_pc_sp %>% mutate(source = "Temperature & spike"),
)

baet_abund_combo$abundance <- scale(baet_abund_combo$abundance) 
baet_abund_combo$source <- ordered(baet_abund_combo$source, levels = c(
  "Temperature", 
  "Temperature & spike", 
  "Temperature & hydropeaking", 
  "Temperature & spike & hydropeaking",
  "Temperature & HFE", 
  "Temperature & spike & HFE", 
  "Temperature & hydropeaking & HFE", 
  "Temperature & spike & hydropeaking & HFE"))  
baet_means_df_abund <- baet_abund_combo %>%
  group_by(source, temperature) %>%
  summarise(mean_abund = mean(abundance), .groups = "drop") %>%
  mutate(is_min = mean_abund == min(mean_abund)) %>%
  ungroup()


b.abund <- ggplot(data = baet_abund_combo, aes(y = abundance, x = as.factor(temperature), color = source, group = source)) +
  stat_summary(fun.data = mean_sdl, geom = "errorbar", linewidth = 0.75, width = 0.2, 
               position = position_dodge(0.5), show.legend = F) +
  geom_point(data = baet_means_df_abund, aes(x = as.factor(temperature), 
                                             y = mean_abund, shape = is_min, 
                                             stroke = 1, color = source, group = source), 
             size = 3, position = position_dodge(0.5), show.legend = FALSE) + # make extinct values a star
  scale_shape_manual(values = c("TRUE" = 8, "FALSE" = 16)) +    scale_x_discrete(labels=c("1" = "Baseline", "1.1" = "+10%",
                            "1.2" = "+20%", "1.5" = "+50%"))+
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1))+
  theme_bw()+
  scale_color_manual(name = " ", labels=label.names, values=paired.colors)+
  xlab("Temperature")+
  labs(y=expression(paste(italic("Baetidae"), "  Scaled Abundance")))+
  theme(text = element_text(size = 13), axis.text.x = element_text(hjust = 1, angle=45, size = 12), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))
# baet biomass
baet_biomass_combo$biomass <- scale(baet_biomass_combo$biomass)
baet_biomass_combo$source <- ordered(baet_biomass_combo$source, levels = c(
  "Temperature", 
  "Temperature & spike", 
  "Temperature & hydropeaking", 
  "Temperature & spike & hydropeaking",
  "Temperature & HFE", 
  "Temperature & spike & HFE", 
  "Temperature & hydropeaking & HFE", 
  "Temperature & spike & hydropeaking & HFE"))

baet_means_df_biomass <- baet_biomass_combo %>%
  group_by(source, temperature) %>%
  summarise(mean_biomass = mean(biomass), .groups = "drop") %>%
  mutate(is_min = mean_biomass == min(mean_biomass)) %>%
  ungroup()

b.biomass <- ggplot(data = baet_biomass_combo, aes(y = biomass, x = as.factor(temperature), color = source, group = source)) +
  stat_summary(fun.data = mean_sdl, geom = "errorbar", linewidth = 0.75, width = 0.2, 
               position = position_dodge(0.5)) +
  stat_summary(fun = mean, geom = "line", size = 0.75, position = position_dodge(0.5)) +
  geom_point(data = baet_means_df_biomass, aes(x = as.factor(temperature), 
                                               y = mean_biomass, shape = is_min, 
                                               stroke = 1, color = source, group = source), 
             size = 3, position = position_dodge(0.5), show.legend = FALSE) + # make extinct values a star
  scale_shape_manual(values = c("TRUE" = 8, "FALSE" = 16)) +  
  scale_x_discrete(labels=c("1" = "Baseline", "1.1" = "+10%",
                            "1.2" = "+20%", "1.5" = "+50%"))+
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1))+
  theme_bw()+
  scale_color_manual(name = " ", labels=label.names, values=paired.colors)+
  xlab("Temperature")+
  labs(y=expression(paste(italic("Baetidae"), " spp. Scaled Biomass")))+
  theme(text = element_text(size = 13), axis.text.x = element_text(hjust = 1, angle=45, size = 11), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

b_plot <- wrap_plots(b.abund, b.biomass) +
  plot_annotation(tag_levels = 'a')+
  plot_layout(guides = 'collect') &
  theme(
    legend.text = element_text(size = 11, lineheight = 0.9),
    legend.spacing.y = unit(0.5, "cm"),
    legend.key.height = unit(1.25, "cm")
  )
# Read the PNG image
mayfly <- readPNG("mayfly.png")
mayfly_grob <- rasterGrob(mayfly, x = 0.825, y = 0.92, width = 0.135, height = 0.135)  # Adjust x, y, width, and height

# Overlay the image on the plot
b.fig <- ggdraw() +
  draw_plot(b_plot) +
  draw_grob(mayfly_grob)

ggsave("bfig.png", plot = b.fig, device = "png", width = 8.44, height = 6, units = c("in"), dpi = "retina")


baet_pc_combo <- bind_rows(
  baet_pc %>% mutate(source = "Temperature"),
  baet_pc_sp %>% mutate(source = "Temperature & spike"),
)
b.pc <- ggplot(data = baet_pc_combo, aes(x = as.factor(temperature), y = percapita, group = source, color = source))+
  stat_summary(fun.data = mean_sdl, geom = "errorbar", linewidth = 0.75, width = 0.2, show.legend = F) +
  stat_summary(fun = mean, geom = "point", size = 3, show.legend = F)+
  stat_summary(fun = mean, geom = "line", size = 0.75, show.legend = F) +
  scale_x_discrete(labels=c("1" = "Baseline", "1.1" = "+10%",
                            "1.2" = "+20%", "1.5" = "+50%"))+
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1))+
  theme_bw()+
  scale_color_manual(name = " ", labels= c("Temperature", 
                                           "Temperature +\nSpike"), values= c("#56B4E9", "#365C8DFF"))+
  xlab("Temperature")+
  labs(y="Adult Biomass (mg)", title = expression(paste(italic("Baetidae"), " spp.")))+
  theme(text = element_text(size = 13), axis.text.x = element_text(hjust = 1, angle=45, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

################################################
# Chironomidae spp Scenarios
chir_abund <- read_csv("CHIR_temp_abund.csv")
chir_pc <- read_csv("CHIR_temp_percapita.csv")
chir_biomass <- read_csv("CHIR_temp_biomass.csv")
chir_abund_sp <- read_csv("CHIR_temp_abund_spike.csv")
chir_pc_sp <- read_csv("CHIR_temp_percapita_spike.csv")
chir_biomass_sp <- read_csv("CHIR_temp_biomass_spike.csv")
chir_abund_t.h <- read_csv("CHIR_temp_hyd_abund.csv")
chir_biomass_t.h <- read_csv("CHIR_temp_hyd_biomass.csv")
chir_abund_t.h_spike <- read_csv("CHIR_temp_hyd_abund_spike.csv")
chir_biomass_t.h_spike <- read_csv("CHIR_temp_hyd_biomass_spike.csv")
chir_abund_HFE <- read_csv("CHIR_temp_abund_HFE.csv")
chir_biomass_HFE <- read_csv("CHIR_temp_biomass_HFE.csv")
chir_abund_HFE_sp <- read_csv("CHIR_temp_abund_HFE_spike.csv")
chir_biomass_HFE.sp <- read_csv("CHIR_temp_biomass_HFE_spike.csv")
chir_abund_HFE_sp <- read_csv("CHIR_temp_abund_HFE_spike.csv")
chir_abund_hyd_HFE <- read_csv("CHIR_temp_hyd_abund_HFE.csv")
chir_biomass_hyd_HFE <- read_csv("CHIR_temp_hyd_biomass_HFE.csv")
chir_abund_hyd_HFE.sp <- read_csv("CHIR_temp_hyd_abund_HFE_spike.csv")
chir_biomass_hyd_HFE.sp <- read_csv("CHIR_temp_hyd_biomass_HFE_spike.csv")

# abundance
chir_abund_combo <- bind_rows(
  chir_abund %>% mutate(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 0),
  chir_abund_sp %>% mutate(source = "Temperature & spike", temp_spike = 1, hydropeaking = 0, HFE = 0),
  chir_abund_t.h %>% mutate(source = "Temperature & hydropeaking", temp_spike = 0, hydropeaking = 1, HFE = 0),
  chir_abund_t.h_spike %>% mutate(source = "Temperature & spike & hydropeaking", temp_spike = 1, hydropeaking = 1, HFE = 0),
  chir_abund_HFE %>% mutate(source = "Temperature & HFE", temp_spike = 0, hydropeaking = 0, HFE = 1),
  chir_abund_HFE_sp %>% mutate(source = "Temperature & spike & HFE", temp_spike = 1, hydropeaking = 0, HFE = 1),
  chir_abund_hyd_HFE %>% mutate(source = "Temperature & hydropeaking & HFE", temp_spike = 0, hydropeaking = 1, HFE = 1),
  chir_abund_hyd_HFE.sp %>% mutate(source = "Temperature & spike & hydropeaking & HFE", temp_spike = 1, hydropeaking = 1, HFE = 1),
)

# biomass
chir_biomass_combo <- bind_rows(
  chir_biomass %>% mutate(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 0),
  chir_biomass_sp %>% mutate(source = "Temperature & spike",temp_spike = 1, hydropeaking = 0, HFE = 0),
  chir_biomass_t.h %>% mutate(source = "Temperature & hydropeaking", temp_spike = 0, hydropeaking = 1, HFE = 0),
  chir_biomass_t.h_spike %>% mutate(source =  "Temperature & spike & hydropeaking", temp_spike = 1, hydropeaking = 1, HFE = 0),
  chir_biomass_HFE %>% mutate(source = "Temperature & HFE", temp_spike = 0, hydropeaking = 0, HFE = 1),
  chir_biomass_HFE.sp %>% mutate(source =  "Temperature & spike & HFE", temp_spike = 1, hydropeaking = 0, HFE = 1),
  chir_biomass_hyd_HFE %>% mutate(source = "Temperature & hydropeaking & HFE", temp_spike = 0, hydropeaking = 1, HFE = 1),
  chir_biomass_hyd_HFE.sp %>% mutate(source = "Temperature & spike & hydropeaking & HFE", temp_spike = 1, hydropeaking = 1, HFE = 1),
)

chir_pc_combo <- bind_rows(
  chir_pc %>% mutate(source = "Temperature"),
  chir_pc_sp %>% mutate(source = "Temperature & spike"),
  #chir_pc_hfe %>% mutate(source = "Temperature & HFE"),
  #chir_pc_hfe_sp %>% mutate(source = "Temperature & spike & HFE")
)

chir_abund_combo$abundance <- scale(chir_abund_combo$abundance) 
chir_abund_combo$source <- ordered(chir_abund_combo$source, levels = c(
  "Temperature", 
  "Temperature & spike", 
  "Temperature & hydropeaking", 
  "Temperature & spike & hydropeaking",
  "Temperature & HFE", 
  "Temperature & spike & HFE", 
  "Temperature & hydropeaking & HFE", 
  "Temperature & spike & hydropeaking & HFE"))
c.abund <- ggplot(data = chir_abund_combo, aes(y = abundance, x = as.factor(temperature), color = source, group = source)) +
  stat_summary(fun.data = mean_sdl, geom = "errorbar", linewidth = 0.75, width = 0.2, 
               position = position_dodge(0.5), show.legend = F) +
  stat_summary(fun = mean, geom = "point", size = 3, position = position_dodge(0.5), show.legend = F) +
  scale_x_discrete(labels=c("1" = "Baseline", "1.1" = "+10%",
                            "1.2" = "+20%", "1.5" = "+50%"))+
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1))+
  theme_bw()+
  scale_color_manual(name = " ", labels=label.names, values=paired.colors)+
  xlab("Temperature")+
  labs(y=expression(paste(italic("Chironomidae"), "  Scaled Abundance")))+
  theme(text = element_text(size = 13), axis.text.x = element_text(hjust = 1, angle=45, size = 12), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))
# baet biomass
chir_biomass_combo$biomass <- scale(chir_biomass_combo$biomass)
chir_biomass_combo$source <- ordered(chir_biomass_combo$source, levels = c(
  "Temperature", 
  "Temperature & spike", 
  "Temperature & hydropeaking", 
  "Temperature & spike & hydropeaking",
  "Temperature & HFE", 
  "Temperature & spike & HFE", 
  "Temperature & hydropeaking & HFE", 
  "Temperature & spike & hydropeaking & HFE"))
c.biomass <- ggplot(data = chir_biomass_combo, aes(y = biomass, x = as.factor(temperature), color = source, group = source)) +
  stat_summary(fun.data = mean_sdl, geom = "errorbar", linewidth = 0.75, width = 0.2, 
               position = position_dodge(0.5)) +
  stat_summary(fun = mean, geom = "point", size = 3, position = position_dodge(0.5)) +
  stat_summary(fun = mean, geom = "line", size = 0.75, position = position_dodge(0.5)) +
  scale_x_discrete(labels=c("1" = "Baseline", "1.1" = "+10%",
                            "1.2" = "+20%", "1.5" = "+50%"))+
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1))+
  theme_bw()+
  scale_color_manual(name = " ", labels=label.names, values=paired.colors)+
  xlab("Temperature")+
  labs(y=expression(paste(italic("Chironomidae"), " spp. Scaled Biomass")))+
  theme(text = element_text(size = 13), axis.text.x = element_text(hjust = 1, angle=45, size = 11), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

c_plot <- wrap_plots(c.abund, c.biomass) +
  plot_annotation(tag_levels = 'a')+
  plot_layout(guides = 'collect') &
  theme(
    legend.text = element_text(size = 11, lineheight = 0.9),
    legend.spacing.y = unit(0.5, "cm"),
    legend.key.height = unit(1.25, "cm")
  )
# Read the PNG image
chiro <- readPNG("chironomidae.png")
chiro_grob <- rasterGrob(chiro, x = 0.825, y = 0.92, width = 0.135, height = 0.135)  # Adjust x, y, width, and height

# Overlay the image on the plot
c.fig <- ggdraw() +
  draw_plot(c_plot) +
  draw_grob(chiro_grob)
ggsave("cfig.png", plot = c.fig, device = "png", width = 8.44, height = 6, units = c("in"), dpi = "retina")


chir_pc_combo <- bind_rows(
  baet_pc %>% mutate(source = "Temperature"),
  baet_pc_sp %>% mutate(source = "Temperature & spike"),
)

c.pc <- ggplot(data = chir_pc_combo, aes(x = as.factor(temperature), y = percapita, group = source, color = source))+
  stat_summary(fun.data = mean_sdl, geom = "errorbar", linewidth = 0.75, width = 0.2, show.legend = F) +
  stat_summary(fun = mean, geom = "point", size = 3, show.legend = F)+
  stat_summary(fun = mean, geom = "line", size = 0.75, show.legend = F) +
  scale_x_discrete(labels=c("1" = "Baseline", "1.1" = "+10%",
                            "1.2" = "+20%", "1.5" = "+50%"))+
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1))+
  theme_bw()+
  scale_color_manual(name = " ", labels= c("Temperature", 
                                           "Temperature +\nSpike"), values= c("#56B4E9", "#365C8DFF"))+
  xlab("Temperature")+
  labs(y="Adult Biomass (mg)", title = expression(paste(italic("Chironomidae"), " spp.")))+
  theme(text = element_text(size = 13), axis.text.x = element_text(hjust = 1, angle=45, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))


###########################################
# G. lacustris Scenaruis
gamm_abund <- read_csv("GAMM_temp_abund.csv")
gamm_pc <- read_csv("GAMM_temp_percapita.csv")
gamm_biomass <- read_csv("GAMM_temp_biomass.csv")
gamm_abund_sp <- read_csv("GAMM_temp_abund_spike.csv")
gamm_pc_sp <- read_csv("GAMM_temp_percapita_spike.csv")
gamm_biomass_sp <- read_csv("GAMM_temp_biomass_spike.csv")
gamm_abund_t.h <- read_csv("GAMM_temp_hyd_abund.csv")
gamm_biomass_t.h <- read_csv("GAMM_temp_hyd_biomass.csv")
gamm_abund_t.h_spike <- read_csv("GAMM_temp_hyd_abund_spike.csv")
gamm_biomass_t.h_spike <- read_csv("GAMM_temp_hyd_biomass_spike.csv")
gamm_abund_HFE <- read_csv("GAMM_temp_abund_HFE.csv")
gamm_biomass_HFE <- read_csv("GAMM_temp_biomass_HFE.csv")
gamm_abund_HFE_sp <- read_csv("GAMM_temp_abund_HFE_spike.csv")
gamm_biomass_HFE.sp <- read_csv("GAMM_temp_biomass_HFE_spike.csv")
gamm_abund_HFE_sp <- read_csv("GAMM_temp_abund_HFE_spike.csv")
gamm_abund_hyd_HFE <- read_csv("GAMM_temp_hyd_abund_HFE.csv")
gamm_biomass_hyd_HFE <- read_csv("GAMM_temp_hyd_biomass_HFE.csv")
gamm_abund_hyd_HFE.sp <- read_csv("GAMM_temp_hyd_abund_HFE_spike.csv")
gamm_biomass_hyd_HFE.sp <- read_csv("GAMM_temp_hyd_biomass_HFE_spike.csv")

# abundance
gamm_abund_combo <- bind_rows(
  gamm_abund %>% mutate(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 0),
  gamm_abund_sp %>%  mutate(source = "Temperature & spike", temp_spike = 1, hydropeaking = 0, HFE = 0),
  gamm_abund_t.h %>% mutate(source = "Temperature & hydropeaking", temp_spike = 0, hydropeaking = 1, HFE = 0),
  gamm_abund_t.h_spike %>% mutate(source = "Temperature & spike & hydropeaking", temp_spike = 1, hydropeaking = 1, HFE = 0),
  gamm_abund_HFE %>% mutate(source = "Temperature & HFE", temp_spike = 0, hydropeaking = 0, HFE = 1),
  gamm_abund_HFE_sp %>% mutate(source = "Temperature & spike & HFE",  temp_spike = 1, hydropeaking = 0, HFE = 1),
  gamm_abund_hyd_HFE %>% mutate(source = "Temperature & hydropeaking & HFE",  temp_spike = 0, hydropeaking = 1, HFE = 1),
  gamm_abund_hyd_HFE.sp %>% mutate(source = "Temperature & spike & hydropeaking & HFE",  temp_spike = 1, hydropeaking = 1, HFE = 1)
)

# biomass
gamm_biomass_combo <- bind_rows(
  gamm_biomass %>%  mutate(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 0),
  gamm_biomass_sp %>%  mutate(source = "Temperature & spike", temp_spike = 1, hydropeaking = 0, HFE = 0),
  gamm_biomass_t.h %>% mutate(source = "Temperature & hydropeaking", temp_spike = 0, hydropeaking = 1, HFE = 0),
  gamm_biomass_t.h_spike %>% mutate(source = "Temperature & spike & hydropeaking", temp_spike = 1, hydropeaking = 1, HFE = 0),
  gamm_biomass_HFE %>% mutate(source = "Temperature & HFE", temp_spike = 0, hydropeaking = 0, HFE = 1),
  gamm_biomass_HFE.sp %>% mutate(source = "Temperature & spike & HFE",  temp_spike = 1, hydropeaking = 0, HFE = 1),
  gamm_biomass_hyd_HFE %>% mutate(source = "Temperature & hydropeaking & HFE",  temp_spike = 0, hydropeaking = 1, HFE = 1),
  gamm_biomass_hyd_HFE.sp %>% mutate(source = "Temperature & spike & hydropeaking & HFE",  temp_spike = 1, hydropeaking = 1, HFE = 1)
)

gamm_pc_combo <- bind_rows(
  gamm_pc %>% mutate(source = "Temperature"),
  gamm_pc_sp %>% mutate(source = "Temperature & spike")
)

gamm_abund_combo$abundance <- scale(gamm_abund_combo$abundance)
gamm_abund_combo$source <- ordered(gamm_abund_combo$source, levels = c(
  "Temperature", 
  "Temperature & spike", 
  "Temperature & hydropeaking", 
  "Temperature & spike & hydropeaking",
  "Temperature & HFE", 
  "Temperature & spike & HFE", 
  "Temperature & hydropeaking & HFE", 
  "Temperature & spike & hydropeaking & HFE"))  

gamm_means_df_abund <- gamm_abund_combo %>%
  group_by(source, temperature) %>%
  summarise(mean_abund = mean(abundance), .groups = "drop") %>%
  mutate(is_min = mean_abund == min(mean_abund)) %>%
  ungroup()


g.abund <- ggplot(data = gamm_abund_combo, aes(y = abundance, x = as.factor(temperature), color = source, group = source)) +
  stat_summary(fun.data = mean_sdl, geom = "errorbar", linewidth = 0.75, width = 0.2, 
               position = position_dodge(0.5), show.legend = F) +
  geom_point(data = gamm_means_df_abund, aes(x = as.factor(temperature), 
                                             y = mean_abund, shape = is_min, 
                                             stroke = 1, color = source, group = source), 
             size = 3, position = position_dodge(0.5), show.legend = FALSE) + # make extinct values a star
  scale_shape_manual(values = c("TRUE" = 8, "FALSE" = 16)) +  
  scale_x_discrete(labels=c("1" = "Baseline", "1.1" = "+10%",
                            "1.2" = "+20%", "1.5" = "+50%"))+
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1))+
  theme_bw()+
  scale_color_manual(name = " ", labels=label.names, values=paired.colors)+
  xlab("Temperature")+
  labs(y=expression(paste(italic("G. lacustris"), " Scaled Abundance")))+
  theme(text = element_text(size = 13), axis.text.x = element_text(hjust = 1, angle=45, size = 12), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

# gamm biomass
gamm_biomass_combo$biomass <- scale(gamm_biomass_combo$biomass)
gamm_biomass_combo$source <- ordered(gamm_biomass_combo$source, levels = c(
  "Temperature", 
  "Temperature & spike", 
  "Temperature & hydropeaking", 
  "Temperature & spike & hydropeaking",
  "Temperature & HFE", 
  "Temperature & spike & HFE", 
  "Temperature & hydropeaking & HFE", 
  "Temperature & spike & hydropeaking & HFE"))  

gamm_means_df_biomass <- gamm_biomass_combo %>%
  group_by(source, temperature) %>%
  summarise(mean_biomass = mean(biomass), .groups = "drop") %>%
  mutate(is_min = mean_biomass == min(mean_biomass)) %>%
  ungroup()

g.biomass <- ggplot(data = gamm_biomass_combo, aes(y = biomass, x = as.factor(temperature), color = source, group = source)) +
  stat_summary(fun.data = mean_sdl, geom = "errorbar", linewidth = 0.75, width = 0.2, 
               position = position_dodge(0.5)) +
  stat_summary(fun = mean, geom = "line", size = 0.75, position = position_dodge(0.5)) +
  geom_point(data = gamm_means_df_biomass, aes(x = as.factor(temperature), 
                                               y = mean_biomass, shape = is_min, 
                                               stroke = 1, color = source, group = source), 
             size = 3, position = position_dodge(0.5), show.legend = FALSE) + # make extinct values a star
  scale_shape_manual(values = c("TRUE" = 8, "FALSE" = 16)) +  
  scale_x_discrete(labels=c("1" = "Baseline", "1.1" = "+10%",
                            "1.2" = "+20%", "1.5" = "+50%"))+
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1))+
  theme_bw()+
  scale_color_manual(name = " ", labels=label.names, values=paired.colors)+
  xlab("Temperature")+
  labs(y=expression(paste(italic("G. lacustris"), " Scaled Biomass")))+
  theme(text = element_text(size = 13), axis.text.x = element_text(hjust = 1, angle=45, size = 11), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

g_plot <- wrap_plots(g.abund, g.biomass) +
  plot_annotation(tag_levels = 'a')+
  plot_layout(guides = 'collect') &
  theme(
    legend.text = element_text(size = 11, lineheight = 0.9),
    legend.spacing.y = unit(0.5, "cm"),
    legend.key.height = unit(1.25, "cm")
  )
# Read the PNG image
gamm <- readPNG("Gammarus.png")
gamm_grob <- rasterGrob(gamm, x = 0.825, y = 0.92, width = 0.135, height = 0.135)  # Adjust x, y, width, and height

# Overlay the image on the plot
g.fig <- ggdraw() +
  draw_plot(g_plot) +
  draw_grob(gamm_grob)

ggsave("gfig.png", plot = g.fig, device = "png", width = 8.44, height = 6, units = c("in"), dpi = "retina")


g.pc <- ggplot(data = gamm_pc_combo, aes(x = as.factor(temperature), y = percapita, group = source, color = source))+
  stat_summary(fun.data = mean_sdl, geom = "errorbar", linewidth = 0.75, width = 0.2, show.legend = F) +
  stat_summary(fun = mean, geom = "point", size = 3, show.legend = F)+
  stat_summary(fun = mean, geom = "line", size = 0.75, show.legend = F) +
  scale_x_discrete(labels=c("1" = "Baseline", "1.1" = "+10%",
                            "1.2" = "+20%", "1.5" = "+50%"))+
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1))+
  theme_bw()+
  scale_color_manual(name = " ", labels= c("Temperature", 
                                           "Temperature +\nSpike"), values= c("#56B4E9", "#365C8DFF"))+
  xlab("Temperature")+
  labs(y="Proportion Adult", title = expression(paste(italic("G. lacustris"))))+
  theme(text = element_text(size = 13), axis.text.x = element_text(hjust = 1, angle=45, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

#############################################
# NZMS Scenarios

nzms_abund <- read_csv("NZMS_temp_abund.csv")
nzms_pc <- read_csv("NZMS_temp_percapita.csv")
nzms_biomass <- read_csv("NZMS_temp_biomass.csv")
nzms_abund_sp <- read_csv("NZMS_temp_abund_spike.csv")
nzms_pc_sp <- read_csv("NZMS_temp_percapita_spike.csv")
nzms_biomass_sp <- read_csv("NZMS_temp_biomass_spike.csv")
nzms_abund_t.h <- read_csv("NZMS_temp_hyd_abund.csv")
nzms_biomass_t.h <- read_csv("NZMS_temp_hyd_biomass.csv")
nzms_abund_t.h_spike <- read_csv("NZMS_temp_hyd_abund_spike.csv")
nzms_biomass_t.h_spike <- read_csv("NZMS_temp_hyd_biomass_spike.csv")
nzms_abund_HFE <- read_csv("NZMS_temp_abund_HFE.csv")
nzms_biomass_HFE <- read_csv("NZMS_temp_biomass_HFE.csv")
nzms_abund_HFE_sp<- read_csv("NZMS_temp_abund_HFE_spike.csv")
nzms_biomass_HFE.sp <- read_csv("NZMS_temp_biomass_HFE_spike.csv")
nzms_abund_HFE_sp<- read_csv("NZMS_temp_abund_HFE_spike.csv")
nzms_abund_hyd_HFE <- read_csv("NZMS_temp_hyd_abund_HFE.csv")
nzms_biomass_hyd_HFE <- read_csv("NZMS_temp_hyd_biomass_HFE.csv")
nzms_abund_hyd_HFE.sp <- read_csv("NZMS_temp_hyd_abund_HFE_spike.csv")
nzms_biomass_hyd_HFE.sp <- read_csv("NZMS_temp_hyd_biomass_HFE_spike.csv")
# 
# nzms_abund <- read_csv("NZMS_temp_abund_2000.csv")
# nzms_pc <- read_csv("NZMS_temp_percapita_2000.csv")
# nzms_biomass <- read_csv("NZMS_temp_biomass_2000.csv")
# nzms_abund_sp <- read_csv("NZMS_temp_abund_spike_2000.csv")
# nzms_pc_sp <- read_csv("NZMS_temp_percapita_spike_2000.csv")
# nzms_biomass_sp <- read_csv("NZMS_temp_biomass_spike_2000.csv")
# nzms_abund_t.h <- read_csv("NZMS_temp_hyd_abund_2000.csv")
# nzms_biomass_t.h <- read_csv("NZMS_temp_hyd_biomass_2000.csv")
# nzms_abund_t.h_spike <- read_csv("NZMS_temp_hyd_abund_spike_2000.csv")
# nzms_biomass_t.h_spike <- read_csv("NZMS_temp_hyd_biomass_spike_2000.csv")
# nzms_abund_HFE <- read_csv("NZMS_temp_abund_HFE_2000.csv")
# nzms_biomass_HFE <- read_csv("NZMS_temp_biomass_HFE_2000.csv")
# nzms_abund_HFE_sp<- read_csv("NZMS_temp_abund_HFE_spike_2000.csv")
# nzms_biomass_HFE.sp <- read_csv("NZMS_temp_biomass_HFE_spike_2000.csv")
# nzms_abund_HFE_sp<- read_csv("NZMS_temp_abund_HFE_spike_2000.csv")
# nzms_abund_hyd_HFE <- read_csv("NZMS_temp_hyd_abund_HFE_2000.csv")
# nzms_biomass_hyd_HFE <- read_csv("NZMS_temp_hyd_biomass_HFE_2000.csv")
# nzms_abund_hyd_HFE.sp <- read_csv("NZMS_temp_hyd_abund_HFE_spike_2000.csv")
# nzms_biomass_hyd_HFE.sp <- read_csv("NZMS_temp_hyd_biomass_HFE_spike_2000.csv")

#abundance
nzms_abund_combo <- bind_rows(
  nzms_abund %>% mutate(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 0),
  nzms_abund_sp %>% mutate(source = "Temperature & spike", temp_spike = 1, hydropeaking = 0, HFE = 0),
  nzms_abund_t.h %>% mutate(source = "Temperature & hydropeaking", temp_spike = 0, hydropeaking = 1, HFE = 0),
  nzms_abund_t.h_spike %>% mutate(source = "Temperature & spike & hydropeaking", temp_spike = 1, hydropeaking = 1, HFE = 0),
  nzms_abund_HFE %>% mutate(source = "Temperature & HFE", temp_spike = 0, hydropeaking = 0, HFE = 1),
  nzms_abund_HFE_sp %>% mutate(source = "Temperature & spike & HFE", temp_spike = 1, hydropeaking = 0, HFE = 1),
  nzms_abund_hyd_HFE %>% mutate(source = "Temperature & hydropeaking & HFE", temp_spike = 0, hydropeaking = 1, HFE = 1),
  nzms_abund_hyd_HFE.sp %>% mutate(source = "Temperature & spike & hydropeaking & HFE", temp_spike = 1, hydropeaking = 1, HFE = 1)
)

# biomass
nzms_biomass_combo <- bind_rows(
  nzms_biomass %>% mutate(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 0),
  nzms_biomass_sp %>% mutate(source = "Temperature & spike", temp_spike = 1, hydropeaking = 0, HFE = 0),
  nzms_biomass_t.h %>% mutate(source = "Temperature & hydropeaking", temp_spike = 0, hydropeaking = 1, HFE = 0),
  nzms_biomass_t.h_spike %>% mutate(source = "Temperature & spike & hydropeaking", temp_spike = 1, hydropeaking = 1, HFE = 0),
  nzms_biomass_HFE %>% mutate(source = "Temperature & HFE", temp_spike = 0, hydropeaking = 0, HFE = 1),
  nzms_biomass_HFE.sp %>% mutate(source = "Temperature & spike & HFE", temp_spike = 1, hydropeaking = 0, HFE = 1),
  nzms_biomass_hyd_HFE %>% mutate(source = "Temperature & hydropeaking & HFE", temp_spike = 0, hydropeaking = 1, HFE = 1),
  nzms_biomass_hyd_HFE.sp %>% mutate(source = "Temperature & spike & hydropeaking & HFE", temp_spike = 1, hydropeaking = 1, HFE = 1)
)

nzms_pc_combo <- bind_rows(
  nzms_pc %>% mutate(source = "Temperature"),
  nzms_pc_sp %>% mutate(source = "Temperature & spike")
)

nzms_abund_combo$abundance <- scale(nzms_abund_combo$abundance)
nzms_abund_combo$source <- ordered(nzms_abund_combo$source, levels = c(
  "Temperature", 
  "Temperature & spike", 
  "Temperature & hydropeaking", 
  "Temperature & spike & hydropeaking",
  "Temperature & HFE", 
  "Temperature & spike & HFE", 
  "Temperature & hydropeaking & HFE", 
  "Temperature & spike & hydropeaking & HFE"))  

n.abund <- ggplot(data = nzms_abund_combo, aes(y = abundance, x = as.factor(temperature), color = source, group = source)) +
  stat_summary(fun.data = mean_sdl, geom = "errorbar", linewidth = 0.75, width = 0.2, 
               position = position_dodge(0.5), show.legend = F) +
  stat_summary(fun = mean, geom = "point", size = 3, position = position_dodge(0.5), show.legend = F) +
  scale_x_discrete(labels=c("1" = "Baseline", "1.1" = "+10%",
                            "1.2" = "+20%", "1.5" = "+50%"))+
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1))+
  theme_bw()+
  scale_color_manual(name = " ", labels=label.names, values=paired.colors)+
  xlab("Temperature")+
  labs(y=expression(paste(italic("P. antipodarum"), " Scaled Abundance")))+
  theme(text = element_text(size = 13), axis.text.x = element_text(hjust = 1, angle=45, size = 12), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))
# nzms biomass
nzms_biomass_combo$biomass <- scale(nzms_biomass_combo$biomass)
nzms_biomass_combo$source <- ordered(nzms_biomass_combo$source, levels = c(
  "Temperature", 
  "Temperature & spike", 
  "Temperature & hydropeaking", 
  "Temperature & spike & hydropeaking",
  "Temperature & HFE", 
  "Temperature & spike & HFE", 
  "Temperature & hydropeaking & HFE", 
  "Temperature & spike & hydropeaking & HFE"))  
n.biomass <- ggplot(data = nzms_biomass_combo, aes(y = biomass, x = as.factor(temperature), color = source, group = source)) +
  stat_summary(fun.data = mean_sdl, geom = "errorbar", linewidth = 0.75, width = 0.2, 
               position = position_dodge(0.5)) +
  stat_summary(fun = mean, geom = "point", size = 3, position = position_dodge(0.5)) +
  stat_summary(fun = mean, geom = "line", size = 0.75, position = position_dodge(0.5)) +
  scale_x_discrete(labels=c("1" = "Baseline", "1.1" = "+10%",
                            "1.2" = "+20%", "1.5" = "+50%"))+
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1))+
  theme_bw()+
  scale_color_manual(name = " ", labels=label.names, values=paired.colors)+
  xlab("Temperature")+
  labs(y=expression(paste(italic("P. antipodarum"), " Scaled Biomass")))+
  theme(text = element_text(size = 11), axis.text.x = element_text(hjust = 1, angle=45, size = 11), 
        axis.text.y = element_text(size = 11), legend.key = element_rect(fill = "transparent"))

n_plot <- wrap_plots(n.abund, n.biomass) +
  plot_annotation(tag_levels = 'a')+
  plot_layout(guides = 'collect') &
  theme(
    legend.text = element_text(size = 11, lineheight = 0.9),
    legend.spacing.y = unit(0.5, "cm"),
    legend.key.height = unit(1.25, "cm")
  )
# Read the PNG image
nzms <- readPNG("nzms.png")
nzms_grob <- rasterGrob(nzms, x = 0.825, y = 0.92, width = 0.135, height = 0.135)  # Adjust x, y, width, and height

# Overlay the image on the plot
n.fig <- ggdraw() +
  draw_plot(n_plot) +
  draw_grob(nzms_grob)

ggsave("nfig.png", plot = n.fig, device = "png", width = 8.44, height = 6, units = c("in"), dpi = "retina")


n.pc <- ggplot(data = nzms_pc_combo, aes(x = as.factor(temperature), y = percapita, group = source, color = source))+
  stat_summary(fun = mean, geom = "point", size = 3, show.legend = F)+
  stat_summary(fun = mean, geom = "line", size = 0.75, show.legend = F) +
  scale_x_discrete(labels=c("1" = "Baseline", "1.1" = "+10%",
                            "1.2" = "+20%", "1.5" = "+50%"))+
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1))+
  theme_bw()+
  scale_color_manual(name = " ", labels= c("Temperature", 
                                           "Temperature +\nSpike"), values= c("#56B4E9", "#365C8DFF"))+
  xlab("Temperature")+
  labs(y="Proportion Adult", title = expression(paste(italic("P. antipodarum"))))+
  theme(text = element_text(size = 13), axis.text.x = element_text(hjust = 1, angle=45, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(color = "transparent"))


# read in Lees Ferry temp and discharge data from 2007 to 2023
temp <- readNWISdv("09380000", "00010", "2007-10-01", "2023-05-01")
temps <- average.yearly.temp(tempdata = temp, temp.column_name = "X_00010_00003", date.column_name = "Date")
# create a timeseries of average temperatures 100 years long
temps <- rep.avg.year(temps, n = 100, change.in.temp = 0, years.at.temp = 0)

# Define the temperature multipliers and labels
scenarios <- data.frame(
  multiplier = c(1, 1.1, 1.2, 1.5),
  label = c("base", "+10%", "+20%", "+50%")
)


# Apply modifications to create the different temperature regimes
temperature.regime <- map2_dfr(scenarios$multiplier, scenarios$label, ~ {
  temps[1:25, ] %>%
    mutate(Temperature = Temperature * .x, se = .y)
})

regime <- ggplot(data = temperature.regime, aes(x = as.Date(dts), y = Temperature, color = se))+
  geom_line(size = 1)+
  theme_bw()+
  scale_color_manual(name = "Warming Regime", labels=c("+10%", "+20%", "+50%", "Baseline"), values=c("#abd9e9", "#fdae61","#d7191c","#2c7bb6"))+
  scale_x_date(date_labels="%b", date_breaks  ="4 month")+
  xlab("")+
  theme(text = element_text(size = 11), axis.text.x = element_text(hjust = 1, size = 10, angle = 0), 
        axis.text.y = element_text(size = 10), legend.key = element_rect(fill = "transparent"))


temp <- readNWISdv("09380000", "00010", "2007-10-01", "2023-05-01")
# calculate average yearly temperatures
temps <- average.yearly.temp(tempdata = temp, temp.column_name = "X_00010_00003", date.column_name = "Date")
# create summertime spike (up to 21 C, then scale from there)
temps$Temperature[16:23] <- c(14, 16, 18, 21, 21, 18, 16, 14)

# create a timeseries of average temperatures 100 years long
temps <- rep.avg.year(temps, n = 100, change.in.temp = 0, years.at.temp = 0)

# Define the temperature multipliers and labels
scenarios <- data.frame(
  multiplier = c(1, 1.1, 1.2, 1.5),
  label = c("base", "+10%", "+20%", "+50%")
)

# Apply modifications to create the different temperature regimes
temperature.regime2 <- map2_dfr(scenarios$multiplier, scenarios$label, ~ {
  temps[1:25, ] %>%
    mutate(Temperature = Temperature * .x, se = .y)
})


#add temperature spike in September
spikeregime <- ggplot(data = temperature.regime2, aes(x = as.Date(dts), y = Temperature, color = se))+
  geom_line(size = 1)+
  theme_bw()+
  scale_color_manual(name = "Warming Regime", labels=c("+10%", "+20%", "+50%", "Baseline"), values=c("#abd9e9", "#fdae61","#d7191c","#2c7bb6"))+
  scale_x_date(date_labels="%b", date_breaks  ="4 month")+
  xlab("")+
  theme(text = element_text(size = 11), axis.text.x = element_text(hjust = 1, size = 10, angle = 0), 
        axis.text.y = element_text(size = 10), legend.key = element_rect(fill = "transparent"))

temp_combo <- bind_rows(
  temperature.regime %>% mutate(source = "Temperature"),
  temperature.regime2 %>% mutate(source = "Temperature & spike"),
)



#add temperature spike in September
regimes <- ggplot(data = temp_combo, aes(x = as.Date(dts), y = Temperature, color = se, linetype = source))+
  geom_line(size = 1)+
  theme_bw()+
  scale_color_manual(name = "Warming Regime", labels=c("+10%", "+20%", "+50%", "Baseline"), values=c("#abd9e9", "#fdae61","#d7191c","#2c7bb6"))+
  scale_x_date(date_labels="%b", date_breaks  ="4 month")+
  xlab("Date")+
  labs(title = "Temperature")+
  guides(linetype=guide_legend(title="Scenario"))+
  theme(text = element_text(size = 13), axis.text.x = element_text(hjust = 1, size = 12.5, angle = 0), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))


Fig2.2 <- wrap_plots(b.pc, h.pc, c.pc, g.pc, n.pc, regimes)+
   plot_layout(guides = 'collect', nrow = 3)+
  plot_annotation(tag_levels = 'a')

ggsave("fig2.2.png", Fig2.2, device = "png", height = 11, width = 15, unit = c("in"), dpi = "retina")

