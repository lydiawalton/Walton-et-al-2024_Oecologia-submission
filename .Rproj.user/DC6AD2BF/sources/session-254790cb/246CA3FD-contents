##############################################################################################
# Title: Investigating the impact of contrasting water and air temperatures on the feeding
#        activity of juvenile Pisaster ochraceus (Part 1)
#
# Author: Lydia Walton
# Last updated: 19-Mar-2024
############################################################################################

dev.off() #Clear open plots
rm(list = ls(all=T)) #Clear environment

#Set working directory
library(here)
#Packages for data manipulation
library(dplyr)
library(tidyr)
#Packages for statistics and modeling
library(AICcmodavg)
library(glmmTMB)
library(EnvStats)
library(broom)
library(DHARMa)
library(emmeans)
#Packages for data visualization
library(ggplot2)
library(ggpmisc)
library(tibble)
library(quantreg)
library(kableExtra)
library(sjPlot)
library(RColorBrewer)
library(extrafont)
library(ggpubr)
library(patchwork)

##Load data files
seastar.collection <- read.csv("MS Data/Pisaster_Collections_June2023.csv")
seastar.adj.heat <- read.csv("MS Data/Pisaster_HeatstressPeriod_June2023.csv")
seastar.recovery <- read.csv("MS Data/Pisaster_RecoveryPeriod_June2023.csv")


###Create theme for plotting
#Font change to Times New Roman as "A"
#windowsFonts(A = windowsFont("Times New Roman"))

#Lydia's theme
LW_theme <- theme_classic() +
  theme(#text = element_text(family = "A"),
        axis.title.x = element_text(vjust = -1, size = 12),
        axis.title.y = element_text(vjust = 2, size = 11),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        panel.background = element_rect(fill = NA, color = "black"))

###########################################
## Data cleaning and manipulation
###########################################

#Subset the data by the variables we want 
subset.collection <- seastar.collection %>% #Collection data
  subset(select = c(star.ID, arm.length, disc.diameter, collection.weight))


subset.adj.heat <- seastar.adj.heat %>% #Adj day (day9) + heat stress
  subset(select=c(star.ID, water.treatment, air.treatment, trial.period, heat.trial.day, total.mussels.consumed, cum.mussels, trial.type, group))


subset.recovery <- seastar.recovery %>% #Recovery
  subset(select=c(star.ID, water.treatment, air.treatment, trial.period, rec.trial.day, rec.mussels.consumed, cum.mussels, trial.type, group))

#Change the column names to match
colnames(subset.adj.heat) <- c("star.ID", "water.treatment", "air.treatment", "trial.period", "day.of.trial", "mussels.consumed","cum.mussels", "trial.type", "group")

colnames(subset.recovery) <- c("star.ID", "water.treatment", "air.treatment", "trial.period", "day.of.trial", "mussels.consumed","cum.mussels", "trial.type", "group")

#Merge dfs together
master.number <- rbind(subset.adj.heat, subset.recovery)

#Add the collection data
master.number <- merge(master.number, subset.collection, by = "star.ID")

#Make sure all of the variables are the right structure
master.number$water.treatment <- as.factor(master.number$water.treatment)
master.number$air.treatment <- as.factor(master.number$air.treatment)

#Remove #67 (4 armed star) and the sea stars that died later
master.number <- master.number %>% 
  filter(star.ID != "SB_2023_ind67", star.ID != "SB_2023_ind9", star.ID != "SB_2023_ind11", star.ID != "SB_2023_ind18", star.ID != "SB_2023_ind19", star.ID != "SB_2023_ind21", star.ID != "SB_2023_ind23", star.ID != "SB_2023_ind53", star.ID != "SB_2023_ind58")

#Remove NAs
master.number <- na.omit(master.number)

##RESCALE body size for modelling
master.number$disc.diameter_z <- arm::rescale(master.number$disc.diameter)


##----------------MODELLING THE RESULTS-----------------------------------------

##############################
# Daily mussels consumed
## Total experimental period
##############################

daily.heatrecovery.modelling <- master.number %>% 
  filter(trial.period != "Adjustment") %>% 
  group_by(star.ID, day.of.trial, trial.type, disc.diameter_z, water.treatment, air.treatment, group) %>% 
  summarise(total.number.mussels = sum(mussels.consumed))

daily.heatrecovery.modelling <- daily.heatrecovery.modelling %>% 
  mutate(Treatment = case_when(
    water.treatment =="15" & air.treatment=="20" ~ "15°C water/20°C air", 
    water.treatment =="15" & air.treatment=="25" ~ "15°C water/25°C air", 
    water.treatment =="15" & air.treatment=="30" ~ "15°C water/30°C air", 
    water.treatment =="20" & air.treatment=="20" ~ "20°C water/20°C air", 
    water.treatment =="20" & air.treatment=="25" ~ "20°C water/25°C air", 
    water.treatment =="20" & air.treatment=="30" ~ "20°C water/30°C air", 
    TRUE ~ "other"
  )) #Add a new Treatment column


#-------------Making the GLMM
daily.heatrecovery.glmm <- glmmTMB(total.number.mussels ~ trial.type*Treatment
                    + disc.diameter_z 
                     + (1|day.of.trial)
                     + (1|group),
                     data = daily.heatrecovery.modelling, 
                     family = poisson)

summary(daily.heatrecovery.glmm) 

daily.heatrecovery.dharma <- simulateResiduals(fittedModel = daily.heatrecovery.glmm, plot = F)

plot(daily.heatrecovery.dharma)

testDispersion(daily.heatrecovery.glmm, alternative = c("two.sided", "greater", "less"),
               plot = T, type = "DHARMa")

tab_model(daily.heatrecovery.glmm,
          show.est = TRUE) #visualizing model outputs

#############################################
# Total mussels within each trial period
## (Heat stress and Recovery)
#############################################

##----------HEAT STRESS-----------------------------------------------------------------
heatstress.modelling <- master.number %>% 
  filter(trial.type == "Heat stress") %>% 
  group_by(star.ID, disc.diameter_z, water.treatment, air.treatment, group) %>% 
  summarise(total.number.mussels = sum(mussels.consumed))

#Making the GLMM
heatstress.glmm <- glmmTMB(total.number.mussels ~ water.treatment*air.treatment +
                                   + disc.diameter_z
                                   + (1|group),
                                   data = heatstress.modelling, 
                                   family = poisson)

summary(heatstress.glmm) 

heatstress.dharma <- simulateResiduals(fittedModel = heatstress.glmm, plot = F)

plot(heatstress.dharma)

tab_model(heatstress.glmm,
          show.est = TRUE)

###------HEAT STRESS + 8 day blocks-----------------------------------------------------------

###############
# DAY 1 to 8
###############
heatstress.day1to8.modelling <- master.number %>% 
  filter(trial.type == "Heat stress") %>% 
  filter(day.of.trial == "1" | day.of.trial == "2" | day.of.trial == "3" | day.of.trial == "4" |
           day.of.trial == "5" | day.of.trial == "6" | day.of.trial == "7" | day.of.trial == "8") %>% 
  group_by(star.ID, disc.diameter_z, water.treatment, air.treatment, group) %>% 
  summarise(total.number.mussels = sum(mussels.consumed))

#Making the GLMM
heatstress.day1to8.glmm <- glmmTMB(total.number.mussels ~ water.treatment*air.treatment +
                             + disc.diameter_z
                           + (1|group),
                           data = heatstress.day1to8.modelling, 
                           family = poisson)

summary(heatstress.day1to8.glmm) 

heatstress.day1to8.dharma <- simulateResiduals(fittedModel = heatstress.day1to8.glmm, plot = F)

plot(heatstress.day1to8.dharma)

tab_model(heatstress.day1to8.glmm,
          show.est = TRUE)

#Create an object with the estimates and CI 
heatstress.day1to8.estimates <- plot_model(heatstress.day1to8.glmm, type = "int")


##############
# Day 9 -16
##############

heatstress.day9to16.modelling <- master.number %>% 
  filter(trial.type == "Heat stress") %>% 
  filter(day.of.trial == "9" | day.of.trial == "10" | day.of.trial == "11" | day.of.trial == "12" |
           day.of.trial == "13" | day.of.trial == "14" | day.of.trial == "15" | day.of.trial == "16") %>% 
  group_by(star.ID, disc.diameter_z, water.treatment, air.treatment, group) %>% 
  summarise(total.number.mussels = sum(mussels.consumed))

#Making the GLMM
heatstress.day9to16.glmm <- glmmTMB(total.number.mussels ~ water.treatment*air.treatment +
                                     + disc.diameter_z
                                   + (1|group),
                                   data = heatstress.day9to16.modelling, 
                                   family = poisson)

summary(heatstress.day9to16.glmm) 

heatstress.day9to16.dharma <- simulateResiduals(fittedModel = heatstress.day9to16.glmm, plot = F)

plot(heatstress.day9to16.dharma)

testDispersion(heatstress.day9to16.dharma, alternative = c("two.sided", "greater", "less"),
               plot = T, type = "DHARMa")

tab_model(heatstress.day9to16.glmm,
          show.est = TRUE)

#Create an object with the estimates and CI 
heatstress.day9to16.estimates <- plot_model(heatstress.day9to16.glmm, type = "int")


##----------------RECOVERY--------------------------------------------------------------

##############
# Day 17-24
##############

recovery.modelling <- master.number %>% 
  filter(trial.type == "Recovery") %>% 
  group_by(star.ID, disc.diameter_z, water.treatment, air.treatment, group) %>% 
  summarise(total.number.mussels = sum(mussels.consumed))

#Making the GLMM
recovery.glmm <- glmmTMB(total.number.mussels ~ water.treatment*air.treatment +
                           + disc.diameter_z
                           + (1|group),
                           data = recovery.modelling, 
                           family = poisson)

summary(recovery.glmm) 

recovery.dharma <- simulateResiduals(fittedModel = recovery.glmm, plot = F)

plot(recovery.dharma)

testDispersion(recovery.glmm, alternative = c("two.sided", "greater", "less"),
               plot = T, type = "DHARMa")

tab_model(recovery.glmm,
          show.est = TRUE)

#Create an object with the estimates and CI 
recovery.estimates <- plot_model(recovery.glmm, type = "int")



#######################################
#
## MAKING PLOTS OF THE DATA
#
#######################################

##Set up the data for plotting
heatrecovery_plotting <- master.number %>% 
  filter(trial.type != "Adjustment") %>% 
  mutate(trial.block = case_when(
    day.of.trial %in% 1:8 ~ "Day 1-8",
    day.of.trial %in% 9:16 ~ "Day 9-16",
    day.of.trial %in% 17:24 ~ "Day 17-24",
    TRUE ~ "other"
  )) %>% 
  group_by(star.ID, trial.block, trial.type, water.treatment, air.treatment) %>% 
  summarise(total.number.mussels = sum(mussels.consumed)) 

##Make trial block an ordered factor
heatrecovery_plotting$trial.block <- factor(heatrecovery_plotting$trial.block,
                                               ordered = TRUE,
                                               levels = c("Day 1-8", "Day 9-16", "Day 17-24"))

##Add a treatment column
heatrecovery_plotting <- heatrecovery_plotting %>% 
  mutate(Treatment = case_when(
    water.treatment =="15" & air.treatment=="20" ~ "15°C water/20°C air", 
    water.treatment =="15" & air.treatment=="25" ~ "15°C water/25°C air", 
    water.treatment =="15" & air.treatment=="30" ~ "15°C water/30°C air", 
    water.treatment =="20" & air.treatment=="20" ~ "20°C water/20°C air", 
    water.treatment =="20" & air.treatment=="25" ~ "20°C water/25°C air", 
    water.treatment =="20" & air.treatment=="30" ~ "20°C water/30°C air", 
    TRUE ~ "other"
  )) 


##################################
# Jitter Plot with Estimates
#
# Panel A
#################################

#Load the model estimate data
PO_modelestimates_byBlock <- read.csv("MS Data/Pisaster_heat-recovery_model-estimates_by-Trial.Block.csv")

#Make appropriate columns factors
PO_modelestimates_byBlock$trial.block <- factor(PO_modelestimates_byBlock$trial.block,
                                                ordered = TRUE,
                                                levels = c("Day 1-8", "Day 9-16", "Day 17-24"))

PO_modelestimates_byBlock$water.treatment <- as.factor(PO_modelestimates_byBlock$water.treatment)
PO_modelestimates_byBlock$air.treatment <- as.factor(PO_modelestimates_byBlock$air.treatment)

#Add treatment column
PO_modelestimates_byBlock <- PO_modelestimates_byBlock %>% 
  mutate(Treatment = case_when(
    water.treatment =="15" & air.treatment=="20" ~ "15°C water/20°C air", 
    water.treatment =="15" & air.treatment=="25" ~ "15°C water/25°C air", 
    water.treatment =="15" & air.treatment=="30" ~ "15°C water/30°C air", 
    water.treatment =="20" & air.treatment=="20" ~ "20°C water/20°C air", 
    water.treatment =="20" & air.treatment=="25" ~ "20°C water/25°C air", 
    water.treatment =="20" & air.treatment=="30" ~ "20°C water/30°C air", 
    TRUE ~ "other"
  ))

###########
# PLOT
###########

#Make the baseplot with the boxes
heatrecovery_JITTERPLOT <- heatrecovery_plotting %>% 
  ggplot() +
  geom_point(aes(x = Treatment, y = total.number.mussels, 
                 shape = trial.block, color = Treatment),
             position = position_jitterdodge(dodge.width = 1),
             alpha = 0.5) +
  scale_color_manual(name = "Treatment (°C)", 
                     values = c("#9DC7C8", "#197BBD", "#05299E" , "#FC6DAB", "#FC814A", "#8D0801")) +
  scale_shape_manual(name = "Trial block",
                     values = c(17, 25, 8)) +
  LW_theme +
  theme(legend.position = "top",
        legend.box = "vertical",
        legend.spacing.y = unit(-0.1, "cm")) +
  labs(x = "", y = "Total # of mussels consumed per individual") +
  scale_y_continuous(breaks = c(0,1,2,3,4,5,6,7,8), limits = c(-0.6,8))

##ADDING THE ESTIMATES

heatrecovery_JITTERPLOT.estimates <- heatrecovery_JITTERPLOT +
  geom_point(data = PO_modelestimates_byBlock,
             aes(x = Treatment, y = estimate, group = trial.block, 
                 color = Treatment),
             position = position_dodge(width = 1), size = 4) +
  scale_color_manual(name = "Treatment (°C)", 
                     values = c("#9DC7C8", "#197BBD", "#05299E" , "#FC6DAB", "#FC814A", "#8D0801")) +
  #Add confidence interval for the model coefficient
  geom_errorbar(data = PO_modelestimates_byBlock,
                aes(x = Treatment, y = estimate, group = trial.block,
                    ymin = conf.low, ymax = conf.high, color = Treatment),
                position = position_dodge(width = 1), width = 0.5, linewidth = 1) +
  #Add sample sizes
  geom_text(label = "n = 11", x = 1, y = -0.7, size = 4) +
  geom_text(label = "n = 12", x = 2, y = -0.7, size = 4) +
  geom_text(label = "n = 7", x = 3, y = -0.7, size = 4) +
  geom_text(label = "n = 12", x = 4, y = -0.7, size = 4) +
  geom_text(label = "n = 12", x = 5, y = -0.7, size = 4) +
  geom_text(label = "n = 9", x = 6, y = -0.7, size = 4) +
  
  ##Add recovery period shading
  annotate("rect", fill = "grey", alpha = 0.3, 
           xmin = 1.2, xmax = 1.5,
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "grey", alpha = 0.3, 
           xmin = 2.2, xmax = 2.5,
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "grey", alpha = 0.3, 
           xmin = 3.2, xmax = 3.5,
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "grey", alpha = 0.3, 
           xmin = 4.2, xmax = 4.5,
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "grey", alpha = 0.3, 
           xmin = 5.2, xmax = 5.5,
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "grey", alpha = 0.3, 
           xmin = 6.2, xmax = 6.5,
           ymin = -Inf, ymax = Inf)

heatrecovery_JITTERPLOT.estimates


##################################################
# Boxplots by trial.block with model estimates
# Panel B, C, D
#################################################


###########################
# BOXPLOT Day 1-8
# Panel B
##########################

###SET UP THE DATA FOR PLOTTING

#For the boxplot
day1to8_boxplot <- heatrecovery_plotting %>% 
  filter(trial.block =="Day 1-8")
  

#For the estimates
day1to8_estimates <- PO_modelestimates_byBlock %>% 
  filter(trial.block =="Day 1-8") 

#Look at the model estimates
summary(heatstress.day1to8.glmm)

#Making the boxplot

day1to8_PLOT_PanelB <- day1to8_boxplot %>% 
  ggplot() +
  geom_point(aes(x = Treatment, y = total.number.mussels, 
                 shape = trial.block, color = Treatment),
             position = position_jitterdodge(dodge.width = 1),
             alpha = 0.5) +
  scale_color_manual(name = "Treatment (°C)", 
                     values = c("#9DC7C8", "#197BBD", "#05299E" , "#FC6DAB", "#FC814A", "#8D0801")) +
  scale_shape_manual(name = "Trial block",
                     values = c(17)) +
  LW_theme +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(x = "Treatment (°C)",y = "Total # of mussels consumed per individual") +
  scale_y_continuous(breaks = c(0,1,2,3,4,5,6,7,8), limits = c(-0.6,8))

##ADDING THE ESTIMATES

day1to8_PLOT_PanelB <- day1to8_PLOT_PanelB +
  geom_point(data = day1to8_estimates,
               aes(x = Treatment, y = estimate, group = trial.block, 
                 color = Treatment),
             position = position_dodge(width = 0.75), size = 4) +
  scale_color_manual(name = "Treatment (°C)", 
                     values = c("#9DC7C8", "#197BBD", "#05299E" , "#FC6DAB", "#FC814A", "#8D0801")) +
  #Add confidence interval for the model coefficient
  geom_errorbar(data = day1to8_estimates,
                aes(x = Treatment, y = estimate, group = trial.block,
                    ymin = conf.low, ymax = conf.high, color = Treatment),
                position = position_dodge(width = 0.75), width = 0.5, linewidth = 1) +
  #Add sample sizes
 # geom_text(label = "N = 11", x = 1, y = -0.7, size = 4, family = "A") +
 # geom_text(label = "N = 12", x = 2, y = -0.7, size = 4, family = "A") +
#  geom_text(label = "N = 7", x = 3, y = -0.7, size = 4, family = "A") +
#  geom_text(label = "N = 12", x = 4, y = -0.7, size = 4, family = "A") +
#  geom_text(label = "N = 12", x = 5, y = -0.7, size = 4, family = "A") +
#  geom_text(label = "N = 9", x = 6, y = -0.7, size = 4, family = "A") +
  #Add plot title
  geom_label(label = "Day 1-8", x = 3.5, y = 8, size = 4)

 
 #Add significant levels 
  #geom_text(label= "**", x=2, y=4.5, size=5, color = "black") + #Air 25 sign overall
  #geom_text(label= "**", x=5, y=6.5, size=5, color = "black") + #Air 25 sign overall
  #geom_text(label= "***", x=3, y=2, size=5, color = "black") + #Air 30 sign overall
  #geom_text(label= "***", x=6, y=4.5, size=5, color = "black") + #Air 30 sign overall
  #geom_text(label= "*", x=5, y=7.5, size=5, color = "red") + #water 20 and air 25 sign interaction
  #geom_text(label= "*", x=6, y=5.5, size=5, color = "red") #water 20 and air 30 sign interaction
  
day1to8_PLOT_PanelB



###########################
# BOXPLOT Day 9-16
# Panel C
##########################

###SET UP THE DATA FOR PLOTTING

#For the boxplot
day9to16_boxplot <- heatrecovery_plotting %>% 
  filter(trial.block =="Day 9-16")

#For the estimates
day9to16_estimates <- PO_modelestimates_byBlock %>% 
  filter(trial.block =="Day 9-16") 

#Look at the model estimates
summary(heatstress.day9to16.glmm)

#Making the boxplot

day9to16_PLOT_PanelC <- day9to16_boxplot %>% 
  ggplot() +
  geom_point(aes(x = Treatment, y = total.number.mussels, 
                 shape = trial.block, color = Treatment),
             position = position_jitterdodge(dodge.width = 1),
             alpha = 0.5) +
  scale_color_manual(name = "Treatment (°C)", 
                     values = c("#9DC7C8", "#197BBD", "#05299E" , "#FC6DAB", "#FC814A", "#8D0801")) +
  scale_shape_manual(name = "Trial block",
                     values = c(25)) +
  LW_theme +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(x = "Treatment (°C)", y = "Total # of mussels consumed per individual") +
  scale_y_continuous(breaks = c(0,1,2,3,4,5,6,7,8), limits = c(-0.6,8))


##ADDING THE ESTIMATES

day9to16_PLOT_PanelC <- day9to16_PLOT_PanelC +
  geom_point(data = day9to16_estimates,
             aes(x = Treatment, y = estimate, group = trial.block, 
                 color = Treatment),
             position = position_dodge(width = 0.75), size = 4) +
  scale_color_manual(name = "Treatment (°C)", 
                     values = c("#9DC7C8", "#197BBD", "#05299E" , "#FC6DAB", "#FC814A", "#8D0801")) +
  #Add confidence interval for the model coefficient
  geom_errorbar(data = day9to16_estimates,
                aes(x = Treatment, y = estimate, group = trial.block,
                    ymin = conf.low, ymax = conf.high, color = Treatment),
                position = position_dodge(width = 0.75), width = 0.5, linewidth = 1) +
  #Add sample sizes
 # geom_text(label = "N = 11", x = 1, y = -0.7, size = 4, family = "A") +
#  geom_text(label = "N = 12", x = 2, y = -0.7, size = 4, family = "A") +
#  geom_text(label = "N = 7", x = 3, y = -0.7, size = 4, family = "A") +
#  geom_text(label = "N = 12", x = 4, y = -0.7, size = 4, family = "A") +
#  geom_text(label = "N = 12", x = 5, y = -0.7, size = 4, family = "A") +
#  geom_text(label = "N = 9", x = 6, y = -0.7, size = 4, family = "A") +
  #Add plot title
  geom_label(label = "Day 9-16", x = 3.5, y = 8, size = 4)
  
#Add significant levels 
     ## Water treatment 20 sign. overall (*)
     ## Air treatment 20 sign. overall (**)
     ## Water treatment 20 and air 30 sign. interaction (*)

day9to16_PLOT_PanelC



###########################
# BOXPLOT Day 17-24
# Panel D
##########################

###SET UP THE DATA FOR PLOTTING

#For the boxplot
day17to24_boxplot <- heatrecovery_plotting %>% 
  filter(trial.block =="Day 17-24")

#For the estimates
day17to24_estimates <- PO_modelestimates_byBlock %>% 
  filter(trial.block =="Day 17-24") 

#Look at the model estimates
summary(recovery.glmm)

#Making the boxplot

day17to24_PLOT_PanelD <- day17to24_boxplot %>% 
  ggplot() +
  geom_point(aes(x = Treatment, y = total.number.mussels, 
                 shape = trial.block, color = Treatment),
             position = position_jitterdodge(dodge.width = 1),
             alpha = 0.5) +
  scale_color_manual(name = "Treatment (°C)", 
                     values = c("#9DC7C8", "#197BBD", "#05299E" , "#FC6DAB", "#FC814A", "#8D0801")) +
  scale_shape_manual(name = "Trial block",
                     values = c(8)) +
  LW_theme +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_rect(fill = "grey91")) +
  labs(x = "Treatment (°C)", y = "Total # of mussels consumed per individual") +
  scale_y_continuous(breaks = c(0,1,2,3,4,5,6,7,8), limits = c(-0.6,8))



##ADDING THE ESTIMATES

day17to24_PLOT_PanelD <- day17to24_PLOT_PanelD +
  geom_point(data = day17to24_estimates,
             aes(x = Treatment, y = estimate, group = trial.block, 
                 color = Treatment),
             position = position_dodge(width = 0.75), size = 4) +
  scale_color_manual(name = "Treatment (°C)", 
                     values = c("#9DC7C8", "#197BBD", "#05299E" , "#FC6DAB", "#FC814A", "#8D0801")) +
  #Add confidence interval for the model coefficient
  geom_errorbar(data = day17to24_estimates,
                aes(x = Treatment, y = estimate, group = trial.block,
                    ymin = conf.low, ymax = conf.high, color = Treatment),
                position = position_dodge(width = 0.75), width = 0.5, linewidth = 1) +
  #Add sample sizes
#  geom_text(label = "N = 11", x = 1, y = -0.7, size = 4, family = "A") +
#  geom_text(label = "N = 12", x = 2, y = -0.7, size = 4, family = "A") +
#  geom_text(label = "N = 7", x = 3, y = -0.7, size = 4, family = "A") +
#  geom_text(label = "N = 12", x = 4, y = -0.7, size = 4, family = "A") +
#  geom_text(label = "N = 12", x = 5, y = -0.7, size = 4, family = "A") +
#  geom_text(label = "N = 9", x = 6, y = -0.7, size = 4, family = "A") +
  #Add plot title
  geom_label(label = "Day 17-24", x = 3.5, y = 8, size = 4)

#Add significant levels 
## Water treatment 20 sign. overall (*)
## Air treatment 20 sign. overall (**)
## Water treatment 20 and air 30 sign. interaction (*)

day17to24_PLOT_PanelD



#####################################
##
### Combine the plots together
##
#####################################

#Panel A: heatrecovery_JITTERPLOT.means
#Panel B: day1to8_PLOT_PanelB
#Panel C: day9to16_PLOT_PanelC
#Panel D: day17to24_PLOT_PanelD


MS_TotalMussels_exp1 <- heatrecovery_JITTERPLOT.estimates / (day1to8_PLOT_PanelB + day9to16_PLOT_PanelC + day17to24_PLOT_PanelD)+
  plot_layout(nrow = 2, heights = c(1, 1)) +
  plot_annotation(tag_levels = "a") 

MS_TotalMussels_exp1


#png("MS Figures/Fig.3_TotalMusselsConsumed_Exp1.png", width = 10, height = 9, units = "in", res = 600)
#MS_TotalMussels_exp1
#dev.off()







