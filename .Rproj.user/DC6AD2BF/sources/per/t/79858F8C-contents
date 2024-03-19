###################################################################################
# Title: Investigating the impact of contrasting water and air temperatures on the 
#        survivorship of juvenile Pisaster ochraceus 
#
# Author: Lydia Walton
# Last updated: 19-Mar-2024
####################################################################################

dev.off() #Clear open plots
rm(list = ls(all=T)) #Clear environment

#Set working directory
library(here)
#Packages for data manipulation
library(dplyr)
library(tidyr)
#Packages for survival analysis
library(survival)
library(ranger)
library(survminer)
#Packages for data visualization
library(ggplot2)
library(ggfortify)
library(patchwork)
library(extrafont)


##Load data files
seastar.collection <- read.csv("MS Data/Pisaster_Collections_June2023.csv")
##Time - survival time in days
##Status - censoring status; 1 = censored, 2 = dead
all_PO_survival <- read.csv("MS Data/Pisaster_Mortality-data_June2023.csv") #For the sea stars

Mussel_survival <- read.csv("MS Data/Mussel_Mortality-data_June2023.csv") #For the mussels


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

#-------------------------------------------------
# Survival analysis for Pisaster ochraceus
all_PO_survival$water.treatment <- as.factor(all_PO_survival$water.treatment)
all_PO_survival$air.treatment <- as.factor(all_PO_survival$air.treatment)

#Compute the survival curve (Kaplan-Meier estimate)
km.fit.PO <- survfit(Surv(time, status) ~ water.treatment + air.treatment, data = all_PO_survival)
print(km.fit.PO)

# Summary of survival curves
summary(km.fit.PO)
# Access to the sort summary table
summary(km.fit.PO)$table

#Access the values returned by survfit()
survival.PO <- data.frame(time = km.fit.PO$time,
                           n.risk = km.fit.PO$n.risk,
                           n.event = km.fit.PO$n.event,
                           n.censor = km.fit.PO$n.censor,
                           surv = km.fit.PO$surv,
                           upper = km.fit.PO$upper,
                           lower = km.fit.PO$lower)
head(survival.PO)

##Plot the survival curve (Kaplan Meier curve)
Fig.kapstar <- ggsurvplot(km.fit.PO, data = all_PO_survival,
                              # conf.int = TRUE, 
                              #pval = TRUE,
                              legend.labs = c("15°C water/20°C air", "15°C water/25°C air",
                                              "15°C water/30°C air", "20°C water/20°C air",
                                              "20°C water/25°C air", "20°C water/30°C air"),
                              legend = c(0.4,0.25),
                              palette = c("#9DC7C8", "#197BBD", "#05299E" , "#FC6DAB", "#FC814A", "#8D0801"),
                              legend.title = "Temperature treatment",
                              break.time.by = 1,
                              xlab = "Time (days)")
# risk.table = TRUE,
# tables.theme = theme_cleantable()

##You can edit the KM curve in ggplot but you have to specify $plot
Fig.kapstar$plot <- Fig.kapstar$plot +
  #Add theme to plot
  ggplot2::theme(#text = element_text(family = "A"),
                 axis.title.x = element_text(vjust = -1, size = 12),
                 axis.title.y = element_text(vjust = 3, size = 12),
                 axis.text = element_text(size = 8),
                 panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 axis.line = element_line(colour = "black"),
                 legend.position = "none") +
  #Add shading for recovery period
  ggplot2::annotate("rect", fill = "grey", alpha = 0.3, 
                    xmin = 17, xmax = 24,
                    ymin = -Inf, ymax = Inf) +
  #Add p-value
  geom_text(label = "p-value = 0.01", x = 0.75, y = 0.1,
            size = 4) +
  #Add animal outline
  add_phylopic(x = 0.75, y = 0.3, img = seastar, alpha = 1, ysize = 0.3)

Fig.kapstar #Look at the plot


##Cox proportional hazards model for the entire data set

cox.PO <- coxph(Surv(time, status)~water.treatment*air.treatment, 
                  data = all_PO_survival)

summary(cox.PO)

##Cox model for day 1-8
PO.cox.subset <- all_PO_survival %>% 
  filter(time != "16")

cox.PO.day1to8 <- coxph(Surv(time, status)~water.treatment, 
                        data = PO.cox.subset)

summary(cox.PO.day1to8) #No significant difference between 30 air treatments 

##Cox model w/out 30air and 15water

PO.cox.no30air15water <- all_PO_survival %>% 
  filter(water.treatment == "15" | air.treatment == "30")

cox.PO.no30air15water <- coxph(Surv(time, status)~water.treatment*air.treatment, 
                               data = PO.cox.no30air15water)

summary(cox.PO.no30air15water) # p = 0.005 (likelihood ratio) and p = 0.01 (logrank test)

#Cox model w/out 30air and 20water

PO.cox.no30air20water <- all_PO_survival %>% 
  filter(water.treatment == "20" | air.treatment == "30")

cox.PO.no30air20water <- coxph(Surv(time, status)~water.treatment*air.treatment, 
                               data = PO.cox.no30air20water)

summary(cox.PO.no30air20water) # p = 0.005 (likelihood ratio) and p = 0.01 (logrank test)

#-------------------------------------------------
# Survival analysis for Mytilus spp.


#Add temp treatment info
seastarformerging <- seastar.collection %>% 
  subset(select = c(star.ID, water.treatment, air.treatment, group))

Mussel_survival <- merge(Mussel_survival, seastarformerging, by = "star.ID")


#Compute the survival curve (Kaplan-Meier estimate)
km.fit.mussel <- survfit(Surv(time, status) ~ water.treatment + air.treatment, data = Mussel_survival)

print(km.fit.mussel)

# Summary of survival curves
summary(km.fit.mussel)
# Access to the sort summary table
summary(km.fit.mussel)$table

#Access the values returned by survfit()
survival.mussel <- data.frame(time = km.fit.mussel$time,
                              n.risk = km.fit.mussel$n.risk,
                              n.event = km.fit.mussel$n.event,
                              n.censor = km.fit.mussel$n.censor,
                              surv = km.fit.mussel$surv,
                              upper = km.fit.mussel$upper,
                              lower = km.fit.mussel$lower)

head(survival.mussel)


##Plot the survival curve (Kaplan Meier curve)
Fig.kapmussel <- 
  ggsurvplot(km.fit.mussel, data = Mussel_survival,
             # conf.int = TRUE, 
             #pval = TRUE,
             surv.median.line = "hv",
             legend.labs = c("15°C water/20°C air", "15°C water/25°C air",
                             "15°C water/30°C air", "20°C water/20°C air",
                             "20°C water/25°C air", "20°C water/30°C air"),
             legend = c(0.4,0.25),
             palette = c("#9DC7C8", "#197BBD", "#05299E" , "#FC6DAB", "#FC814A", "#8D0801"),
             #palette = c("#FF88DC", "#88CCEE", "#661100", "#117733",  "#332288", 
              #           "#DDCC77"),
             legend.title = "Temperature treatment",
             break.time.by = 1,
             xlab = "Time (days)") 

##You can edit the KM curve in ggplot but you have to specify $plot
Fig.kapmussel$plot <- Fig.kapmussel$plot +
  #Add theme to plot
  ggplot2::theme(#text = element_text(family = "A"),
                 axis.title.x = element_text(vjust = -1, size = 12),
                 axis.title.y = element_text(vjust = 3, size = 12),
                 panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 axis.line = element_line(colour = "black"),
                 legend.position = "bottom",
                 legend.text = element_text(size = 10)) +
  #Add shading for recovery period
  ggplot2::annotate("rect", fill = "grey", alpha = 0.32, 
                    xmin = 17, xmax = 24,
                    ymin = -Inf, ymax = Inf) +
  #Add p-value
  geom_text(label = "p-value = 0.04", x = 1.25, y = 0.075,
            size = 4) +
  #Add animal outline
  add_phylopic(x = 1.5, y = 0.25, img = mussel, alpha = 1, ysize = 0.2)

 Fig.kapmussel

##Cox proportional hazards model
cox.mussel <- coxph(Surv(time, status)~water.treatment*air.treatment, data = Mussel_survival)

summary(cox.mussel)
##exp(coef) describes the hazard ratio = 0.90 - going from 15 to 20 C results in ~10% reduction in hazard
## HR < 1: reduction in hazard (protective)
## HR > 1: increase in hazard
## HR = 1: no effect



#---------------------
## Making the MS figure

MS_Surv_exp1 <- Fig.kapstar$plot + Fig.kapmussel$plot +
  plot_layout(nrow = 2, heights = c(0.75, 1)) +
  plot_annotation(tag_levels = "a") 

MS_Surv_exp1


#png("MS Figures/Fig.2_Mortality_Exp1.png", width = 8, height = 8, units = "in", res = 600)
#MS_Surv_exp1
#dev.off()

