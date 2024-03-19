###################################################################################
# Title: Comparison of the feeding activity and survivorship of Pisaster ochraceus
#        between the feeding activity experiment and the metabolic rate experiment 
#
# Author: Lydia Walton
# Last updated: 19-Mar-2024
##################################################################################

dev.off() #Clear open plots
rm(list = ls(all=T)) #Clear environment

#Set working directory
library(here)
#Data manipulation
library(dplyr)
library(tidyr)
#Data visualization
library(ggplot2)
library(ggpmisc)
library(tibble)
library(kableExtra)
library(sjPlot)
library(RColorBrewer)
library(extrafont)
library(ggpubr)
library(patchwork)
#Modeling
library(quantreg)
library(AICcmodavg)
library(glmmTMB)
library(EnvStats)
library(broom)
library(DHARMa)
library(emmeans)
#Survival analysis
library(survival)
library(ranger)
library(ggfortify)
library(survminer)

#Load data
PO.collection.exp2 <- read.csv("MS Data/Pisaster_Collections_August2023.csv")
PO.acclimation.exp2 <- read.csv("MS Data/Pisaster_Acclimation_August2023.csv")

PO.collection.exp1 <- read.csv("MS Data/Pisaster_Collections_June2023.csv")
PO.heatstress.exp1 <- read.csv("MS Data/Pisaster_HeatstressPeriod_June2023.csv")

PO.survival.exp2 <- read.csv("MS Data/Pisaster_Mortality-data_August2023.csv")
PO.survival.exp1 <- read.csv("MS Data/Pisaster_Mortality-data_June2023.csv")

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

##Get animal silhouettes from phylopic
library(rphylopic)
seastar <- pick_phylopic(name = "Asteriidae", n = 1)

mussel <- pick_phylopic(name = "Mytilus galloprovincialis", n= 1)

# PhyloPic has a lot of contributors and we should acknowledge 
# their work. You can get data about images using get_attribution

# Get valid uuid
star.cont <- get_uuid(name = "Asteriidae", n=1)
# Get attribution data for uuid
get_attribution(uuid = star.cont)

mussel.cont <- get_uuid(name = "Mytilus galloprovincialis", n= 1)
get_attribution(uuid = mussel.cont)

#######################
# Survival Comparison
#######################

#Add new column to separate data and filter to water20 and air25/30 
##(Treatments found across both experiments)
PO.survival.exp1 <- PO.survival.exp1 %>% 
  add_column(Month = "June") %>% 
  filter(water.treatment == "20") %>% 
  filter(air.treatment == "25" | air.treatment =="30") %>% 
  subset(select = c(star.ID, air.treatment, disc.diameter, collection.weight, time, status, Month))

#New column names to match exp2
colnames(PO.survival.exp1) <- c("star.ID", "air.treatment", "disc.diameter", "wet.weight",
                                "time", "status", "Month")

##Add new column to exp2
PO.survival.exp2.1 <- PO.survival.exp2 %>% 
  add_column(Month = "August") %>% 
  filter(food.treatment =="Fed") %>% 
  subset(select = c(star.ID, air.treatment, disc.diameter, wet.weight, time, status, Month))

#merge the data frames
PO.survival.master <- rbind(PO.survival.exp1, PO.survival.exp2.1)

#Set the months to a specific order
PO.survival.master$Month <- factor(PO.survival.master$Month, order = TRUE,
                                   levels = c("June", "August"))


##----------------KM Curve

#Compute the survival curve (Kaplan-Meier estimate)
km.fit.comparison <- survfit(Surv(time, status) ~ air.treatment + Month, 
                             data = PO.survival.master)
print(km.fit.comparison)

# Summary of survival curves
summary(km.fit.comparison)
# Access to the sort summary table
summary(km.fit.comparison)$table

#Access the values returned by survfit()
survival.comparison <- data.frame(time = km.fit.comparison$time,
                                  n.risk = km.fit.comparison$n.risk,
                                  n.event = km.fit.comparison$n.event,
                                  n.censor = km.fit.comparison$n.censor,
                                  surv = km.fit.comparison$surv,
                                  upper = km.fit.comparison$upper,
                                  lower = km.fit.comparison$lower
)
head(survival.comparison)

#Plotting the KM curve
Fig.kapstar.comparison <- ggsurvplot(km.fit.comparison, data = PO.survival.master,
                                     #conf.int = TRUE, 
                                     #pval = TRUE,
                                     legend.labs = c("25°C air & June", "25°C air & August", 
                                                     "30°C air & June", "30°C air & August"),
                                     #Legend.title = "Water temperature treatment (°C)",
                                     legend = c(0.125,0.25),
                                     palette = c("#FC814A","#8FD5A6", "#8D0801","#0C8346"),
                                     legend.title = "Treatment",
                                     break.time.by = 1,
                                     xlab = "Time (days)")
# risk.table = TRUE,
# tables.theme = theme_cleantable(),


Fig.kapstar.comparison$plot <- Fig.kapstar.comparison$plot +
  #Add theme to plot
  ggplot2::theme(#text = element_text(family = "A"),
                 axis.title.x = element_text(vjust = -1, size = 10),
                 axis.title.y = element_text(vjust = 3, size = 10),
                 axis.text = element_text(size = 8),
                 panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 axis.line = element_line(colour = "black"),
                 legend.position = "none",
                 legend.text = element_text(size = 10)) +
  #Add p-value
  geom_text(label = "p-value < 0.001", x = 1, y = 0.1,
            size = 4, family = "A") +
  #Add animal outline
  add_phylopic(x = 1, y = 0.3, img = seastar, alpha = 1, ysize = 0.3)

Fig.kapstar.comparison

##Cox proportional hazards model for both factors

cox.PO <- coxph(Surv(time, status) ~ air.treatment*Month, 
                data = PO.survival.master)

summary(cox.PO) #p < 0.001

##For air treatment
cox.PO.air <- coxph(Surv(time, status) ~ air.treatment, 
                data = PO.survival.master)

summary(cox.PO.air) #p = 0.04

##For month
cox.PO <- coxph(Surv(time, status) ~ Month, 
                data = PO.survival.master)

summary(cox.PO) #p < 0.001




#####################################
# Feeding activity (total mussels)
#####################################


#Subset the data with the variables we want
PO.collection.exp1.subset <- PO.collection.exp1 %>% #Collection data
  subset(select = c(star.ID, disc.diameter, collection.weight))


PO.heatstress.exp1.subset <- PO.heatstress.exp1 %>% 
  filter(trial.type == "Heat stress") %>% 
  subset(select=c(star.ID, water.treatment, air.treatment, heat.trial.day, total.mussels.consumed, cum.mussels))

#Merge the collection and heat stress df together
exp1.master <- merge(PO.collection.exp1.subset, PO.heatstress.exp1.subset, by = "star.ID")

#Remove the dead sea stars by getting rid of the NAs 
exp1.master <- na.omit(exp1.master)

#Remove #67 (4 armed star) and the sea stars that died later
exp1.master <- exp1.master %>% 
  filter(star.ID != "SB_2023_ind67", star.ID != "SB_2023_ind9", star.ID != "SB_2023_ind11", star.ID != "SB_2023_ind18", star.ID != "SB_2023_ind19", star.ID != "SB_2023_ind21", star.ID != "SB_2023_ind23", star.ID != "SB_2023_ind53", star.ID != "SB_2023_ind58")

#add a new column
exp1.master <- exp1.master %>% 
  add_column(Month = "June") %>% 
  filter(water.treatment == "20") %>% 
  filter(air.treatment == "25" | air.treatment =="30") %>% 
  subset(select = c(star.ID, air.treatment, heat.trial.day, disc.diameter,
                    collection.weight, total.mussels.consumed, cum.mussels, Month))

##----------Load experiment 2 data

#Merge necessary dfs
exp2.merge <- merge(PO.acclimation.exp2, PO.survival.exp2, 
                    by = c("star.ID", "food.treatment", "air.treatment"))

##Subset and filter the data 
exp2.master <- exp2.merge %>% 
  add_column(Month = "August") %>% 
  filter(food.treatment == "Fed") %>%  #Only want the fed treatments
  filter(status == "1") %>% #Only want sea stars that survived
  subset(select = c(star.ID, air.treatment, trial.day, disc.diameter, wet.weight,
                    mussels.consumed, cum.mussels, Month))


##--------Rename columns and merge the dfs
colnames(exp1.master) <- c("star.ID", "air.treatment", "trial.day", "disc.diameter", "wet.weight", "mussels.consumed", "cum.mussels", "Month")

colnames(exp2.master) <- c("star.ID", "air.treatment", "trial.day", "disc.diameter", "wet.weight", "mussels.consumed", "cum.mussels", "Month")

#bind df together
exp.compare.master <- rbind(exp1.master, exp2.master)

#Set the months to a specific order
exp.compare.master$Month <- factor(exp.compare.master$Month, order = TRUE,
                                   levels = c("June", "August"))


##Model the data with a glmm

##RESCALE body size for lm's
exp.compare.master$disc.diameter_z <- arm::rescale(exp.compare.master$disc.diameter)

#Df for modelling
comparison.modelling <- exp.compare.master %>% 
  mutate(Treatment = case_when(
    Month =="June" & air.treatment=="25" ~ "25°C air & June", 
    Month =="June" & air.treatment=="30" ~ "30°C air & June",  
    Month =="August" & air.treatment=="25" ~ "25°C air & August", 
    Month =="August" & air.treatment=="30" ~ "30°C air & August",  
    TRUE ~ "other"
  )) %>%
  group_by(star.ID, disc.diameter_z, air.treatment, Month, Treatment) %>% 
  summarise(total.number.mussels = sum(mussels.consumed))

#Is there a correlation between body size and mussels variables?
cor(comparison.modelling$total.number.mussels, comparison.modelling$disc.diameter_z)
#Weak negative correlation

#----------MODEL

#Look at the total number of mussels consumed during the heat stress period by the different treatment groups 
compare.glmm <- glmmTMB(total.number.mussels ~ air.treatment*Month + (1|disc.diameter_z),
                     data = comparison.modelling, 
                     family = poisson)

summary(compare.glmm) #AIC for this model is 150.1

#Test the model to see if there are issues
compare.dharma <- simulateResiduals(fittedModel = compare.glmm, plot = F)

plot(compare.dharma)

## ------------ANOVA
#Normality
compare.norm <- comparison.modelling %>% 
  group_by(Treatment) 

shapiro.test(compare.norm$total.number.mussels) #NORMAL

#ANOVA
compare.feeding.aov <- aov(total.number.mussels ~ Treatment,
                        data = comparison.modelling)

summary(compare.feeding.aov)

#Look at equal variance
bartlett.test(total.number.mussels ~ Treatment, 
              data = comparison.modelling)


#TukeyHSD

compare.feeding.tukey <- TukeyHSD(compare.feeding.aov)

compare.feeding.tukey

# compact letter display
cld <- multcompView::multcompLetters4(compare.feeding.aov, compare.feeding.tukey)
print(cld)


###############################
# BOXPLOT
#############################

#Make a df with the total mussels per treatment  
compare.total.mussels <- exp.compare.master %>% 
  mutate(Treatment = case_when(
    Month =="June" & air.treatment=="25" ~ "25°C air & June", 
    Month =="June" & air.treatment=="30" ~ "30°C air & June",  
    Month =="August" & air.treatment=="25" ~ "25°C air & August", 
    Month =="August" & air.treatment=="30" ~ "30°C air & August",  
    TRUE ~ "other"
  )) %>%
  group_by(star.ID, air.treatment, Month, Treatment) %>% 
  summarise(total.number.mussels = sum(mussels.consumed))

#Set the months to a specific order
compare.total.mussels$Treatment <- factor(compare.total.mussels$Treatment, order = TRUE,
                                       levels = c("25°C air & June", "30°C air & June",
                                                  "25°C air & August", "30°C air & August"))

#Make the plot
compare.totalmussel.Fig <- compare.total.mussels %>% 
  ggplot() +
  geom_boxplot(aes(x = Month, y = total.number.mussels, fill = Treatment),
               outlier.shape = NA) +
  scale_fill_manual(name = "Treatment",
                    values = c("#FC814A","#8FD5A6", "#8D0801","#0C8346")) +
  LW_theme +
  theme(legend.position = "bottom",
        axis.title.y = element_text(vjust = 2, size = 10)) +
  labs(x = "Month", y = "Total number of mussels consumed") +
  scale_y_continuous(breaks = c(2,3,4,5,6,7,8,9,10,11,12,13), limits = c(2,13.5))+
  ##Sample sizes 
  geom_text(label= "N = 12", x=0.815, y= 9.5, size=4, family = "A") +
  geom_text(label= "N = 9", x=1.19, y=5.5, size=4, family = "A") +
  geom_text(label= "N = 8", x=1.815, y=9.5, size= 4, family = "A") +
  geom_text(label= "N = 3", x=2.19, y=6.5, size= 4, family = "A") +
  #TUKEY results
  geom_text(label= "a", x=0.815, y= 13.6, size=4, color = "#FC814A", family = "A") +
  geom_text(label= "b", x=1.19, y=11.6, size=4, color = "#8FD5A6", family = "A") +
  geom_text(label= "ab", x=1.815, y=11.6, size= 4, color = "#8D0801", family = "A") +
  geom_text(label= "ab", x=2.19, y=8.6, size= 4, color = "#0C8346", family = "A")

compare.totalmussel.Fig


##############
# MS FIGURE
#############


MS_comparison_Fig <- Fig.kapstar.comparison$plot + compare.totalmussel.Fig + 
  plot_layout(nrow = 2, heights = c(1, 1)) +
  plot_annotation(tag_levels = "a") +
  theme(legend.position = "bottom")


MS_comparison_Fig

png("MS Figures/Fig.6_Exp-Mortality-FeedingComparison.png", width = 8, height = 8, units = "in", res = 600)
MS_comparison_Fig
dev.off()







