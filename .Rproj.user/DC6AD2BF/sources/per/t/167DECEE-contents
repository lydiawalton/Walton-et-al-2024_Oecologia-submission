##################################################################################
# Title: Investigating the impact of air temperature treatments on the feeding
#        activity of juvenile Pisaster ochraceus (Experiment 2)
#
# Author: Lydia Walton
# Last updated: 19-Mar-2024
#################################################################################

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
PO.collection <- read.csv("MS Data/Pisaster_Collections_August2023.csv")
PO.acclimation <- read.csv("MS Data/Pisaster_Acclimation_August2023.csv")
PO.survivorship <- read.csv("MS Data/Pisaster_Mortality-data_August2023.csv")

###Create theme for plotting
#Font change to Times New Roman as "A"
#windowsFonts(A = windowsFont("Times New Roman"))

#Lydia's theme
LW_theme <- theme_classic() +
  theme(#text = element_text(family = "A"),
        axis.title.x = element_text(vjust = -1, size = 12),
        axis.title.y = element_text(vjust = 2, size = 11),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        panel.background = element_rect(fill = NA, color = "black"))


###########################################
## Data cleaning and manipulation
###########################################

#Make labels for plots
food.labels <- c("+mussels", "-mussels")
names(food.labels) <- c("+mussels", "+mussels")

air.labels <- c("25°C air treatment", "30°C air treatment")
names(air.labels) <- c("25", "30")

#Make sure air treatments are factors
PO.acclimation$air.treatment <- as.factor(PO.acclimation$air.treatment)
PO.collection$air.treatment <- as.factor(PO.collection$air.treatment)


#Make a new df for the cumulative consumption of mussels for living sea stars
PO.cumulative.consumption <- PO.acclimation %>% 
  filter(food.treatment == "Fed") %>%  #Only want the fed treatments
  subset(select = c(star.ID, food.treatment, air.treatment, trial.day, mussels.consumed, cum.mussels)) %>% #Subset with the variables we want
  mutate(Treatment = case_when(
    food.treatment =="Fed" & air.treatment=="25" ~ "+mussels/25°C air", 
    food.treatment =="Fed" & air.treatment=="30" ~ "+mussels30°C air", 
    TRUE ~ "other"
  ))

##Merge this df with the survival data
PO.cumulative.consumption <- merge(PO.cumulative.consumption, PO.survivorship,
                                   by = c("star.ID", "food.treatment", "air.treatment"))

#Omit stars after they have died
PO.cumulative.consumption <- na.omit(PO.cumulative.consumption)


#############################
# Cumulative consumption
# ALIVE sea stars
############################

##SAMPLE SIZES FOR PLOT = Alive/25 (n=8), Alive/30 (n=3), Dead/25 (n=0), Dead/30 (n=2)

#Make a df with only the alive sea stars
PO.cum.consumption.ALIVE <- PO.cumulative.consumption %>% 
  filter(status == "1") %>% 
  group_by(Treatment, trial.day) %>% 
  summarise(total.cum.mussels = sum(cum.mussels))

#Make the plot

VW.cumconsumption.ALIVE.plot <- PO.cum.consumption.ALIVE %>% 
  ggplot(aes(x = trial.day, y = total.cum.mussels, color = Treatment)) +
  geom_line() +
  geom_point() +
  scale_color_manual(labels = c("+mussels/25°C air","+mussels/30°C air"),
                     values = c("#8FD5A6","#0C8346")) +
  LW_theme +
  theme(legend.position = "bottom") +
  labs(x = "Day of trial", y = "Cumulative number of mussels consumed") +
  scale_x_continuous(breaks = seq(0, 16, by = 2)) +
  scale_y_continuous(breaks = seq(0, 80, by = 5)) +
  geom_text(label= "n = 8", x=14.5, y=68, size=4, color = "black") +
  geom_text(label= "n = 3", x=15.5, y=20, size=4, color = "black")


VW.cumconsumption.ALIVE.plot


###############################
# Total consumption
# ALIVE sea stars 
###############################

#Set up the data for total number of mussels
PO.total.cons.ALIVE <- PO.cumulative.consumption %>% 
  filter(status == "1") %>% 
  group_by(star.ID, Treatment) %>% 
  summarise(total.num.mussels = sum(mussels.consumed))


#Test for normality and equal variance
##NORMALITY
with(PO.total.cons.ALIVE, shapiro.test(total.num.mussels[Treatment == "+mussels/25°C air"]))
with(PO.total.cons.ALIVE, shapiro.test(total.num.mussels[Treatment == "+mussels/30°C air"]))

##EQUAL VARIANCE
var.test(total.num.mussels ~ Treatment, data = PO.total.cons.ALIVE)

#Assumptions are met so we can now do a t-test
t.test(total.num.mussels ~ Treatment, data = PO.total.cons.ALIVE,
       var.equal = TRUE)
### p = 0.05947 No sign. diff but trend towards lower feeding in 30air (but also way lower sample size)

#Make the plot
set.seed(838)
VW.total.cons.ALIVE.plot <- PO.total.cons.ALIVE %>% 
  ggplot(aes(x = Treatment, y = total.num.mussels, color = Treatment)) +
  geom_jitter(width = 0.15) +
  scale_color_manual(labels = c("+mussels/25°C air","+mussels/30°C air"),
                     values = c("#8FD5A6","#0C8346")) +
  LW_theme +
  theme(legend.position = "none") +
  labs(x = "Treatment", y = "Total number of mussels consumed") +
  #geom_text(label= "N = 8", x=1, y=4, size=4, family = "A", color = "black") +
  #geom_text(label= "N = 3", x=2, y=9, size=4, family = "A", color = "black") +
  scale_y_continuous(breaks = seq(0, 12, by = 1)) +
  geom_label(label = "ns", x=1.5, y=11, size=4, color = "black", fill = "white", family = "A")

VW.total.cons.ALIVE.plot


##################################
#SUPP FIGURE - Feeding activity
#################################

MS_FeedingActivity_exp2 <- VW.cumconsumption.ALIVE.plot + VW.total.cons.ALIVE.plot +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom") 

MS_FeedingActivity_exp2


#png("MS Figures/Fig.S1_FeedingActivity_Exp2.png", width = 8, height = 6, units = "in", res = 600)
#MS_FeedingActivity_exp2
#dev.off()


############################################
# Cumulative consumption per individual
# Alive and Dead
############################################

#Make a df with only the alive sea stars
PO.cum.consumption.ALIVE.ind <- PO.cumulative.consumption %>% 
  filter(status == "1") %>% 
  group_by(star.ID, Treatment, trial.day)

PO.cum.consumption.ALIVE.ind$star.ID <- 
  factor(PO.cum.consumption.ALIVE.ind$star.ID,
         ordered = TRUE,
         levels = c("VW_2023_ind5", "VW_2023_ind16", "VW_2023_ind17", "VW_2023_ind22",
                    "VW_2023_ind32", "VW_2023_ind35", "VW_2023_ind46", "VW_2023_ind50",
                    "VW_2023_ind4", "VW_2023_ind11", "VW_2023_ind12"))

#Make the plot

VW.cumconsumption.ALIVE.plot.ind <- PO.cum.consumption.ALIVE.ind %>% 
  ggplot(aes(x = trial.day, y = cum.mussels, color = Treatment)) +
  geom_line() +
  geom_point() +
  facet_wrap(~star.ID)+
  scale_color_manual(name = "Treatment",
                     values = c("#8FD5A6","#0C8346")) +
  LW_theme +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 8)) +
  labs(x = "Day of trial", y = "Cumulative number of mussels consumed") +
  ggtitle("Survived") +
  scale_x_continuous(breaks = seq(0, 16, by = 2)) +
  scale_y_continuous(breaks = seq(0, 12, by = 1)) 

VW.cumconsumption.ALIVE.plot.ind


#Make a df with only the alive sea stars
PO.cum.consumption.DEAD.ind <- PO.cumulative.consumption %>% 
  filter(status == "2") %>% 
  group_by(star.ID, Treatment, trial.day)

PO.cum.consumption.DEAD.ind$star.ID <- 
  factor(PO.cum.consumption.DEAD.ind$star.ID,
         ordered = TRUE,
         levels = c("VW_2023_ind18", "VW_2023_ind25", "VW_2023_ind27", "VW_2023_ind41",
                    "VW_2023_ind47", "VW_2023_ind53", 
                    "VW_2023_ind2", "VW_2023_ind6","VW_2023_ind10", "VW_2023_ind33",
                    "VW_2023_ind38", "VW_2023_ind39", "VW_2023_ind44", "VW_2023_ind45",
                    "VW_2023_ind48", "VW_2023_ind49", "VW_2023_ind55"))

#Make the plot

VW.cumconsumption.DEAD.plot.ind <- PO.cum.consumption.DEAD.ind %>% 
  ggplot(aes(x = trial.day, y = cum.mussels, color = Treatment)) +
  geom_line() +
  geom_point() +
  facet_wrap(~star.ID)+
  scale_color_manual(name = "Treatment",
                     values = c("#8FD5A6","#0C8346")) +
  LW_theme +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 7)) +
  labs(x = "Day of trial", y = "Cumulative number of mussels consumed") +
  ggtitle("Died") +
  scale_x_continuous(breaks = seq(0, 16, by = 2)) +
  scale_y_continuous(breaks = seq(0, 5, by = 1)) 


VW.cumconsumption.DEAD.plot.ind



##################################
# Supplementary figure
#################################

Supp_Ind.FeedingActivity_exp2 <- VW.cumconsumption.ALIVE.plot.ind + VW.cumconsumption.DEAD.plot.ind +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom") 

Supp_Ind.FeedingActivity_exp2


#png("MS Figures/Fig.S2_FeedingActivity_Exp2.png", width = 10, height = 8, units = "in", res = 600)
#Supp_Ind.FeedingActivity_exp2
#dev.off()


