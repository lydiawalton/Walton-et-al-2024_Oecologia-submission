###############################################################################
# Title: Investigating the impact of contrasting heat treatments and food
#        availability on the survivorship of juvenile Pisaster ochraceus 
#
# Author: Lydia Walton
# Last updated: 19-Mar-2024
###############################################################################

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
PO.collection <- read.csv("MS Data/Pisaster_Collections_August2023.csv")
##Time - survival time in days
##Status - censoring status; 1 = censored, 2 = dead
PO.survivorship <- read.csv("MS Data/Pisaster_Mortality-data_August2023.csv")
Mussel.survivorship <- read.csv("MS Data/Mussel_Mortality-data_August2023.csv")


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

#Compute the survival curve (Kaplan-Meier estimate)
km.fit.PO.exp2 <- survfit(Surv(time, status) ~ food.treatment + air.treatment, data = PO.survivorship)
print(km.fit.PO.exp2)

# Summary of survival curves
summary(km.fit.PO.exp2)
# Access to the sort summary table
summary(km.fit.PO.exp2)$table

#Access the values returned by survfit()
survival.PO.exp2 <- data.frame(time = km.fit.PO.exp2$time,
                          n.risk = km.fit.PO.exp2$n.risk,
                          n.event = km.fit.PO.exp2$n.event,
                          n.censor = km.fit.PO.exp2$n.censor,
                          surv = km.fit.PO.exp2$surv,
                          upper = km.fit.PO.exp2$upper,
                          lower = km.fit.PO.exp2$lower
)
head(survival.PO.exp2)


##Plot the survival curve (Kaplan Meier curve)

Fig.kapstar.exp2 <- ggsurvplot(km.fit.PO.exp2, data = PO.survivorship,
                             #conf.int = TRUE, 
                             #pval = TRUE,
                             surv.median.line = "hv",
                             #risk.table = TRUE,
                             legend.labs = c("+mussels/25°C air","+mussels/30°C air", 
                                             "-mussels/25°C air", "-mussels/30°C air"),
                             legend.title = "Treatment",
                             legend = c(0.8,0.85),
                             palette = c("#8FD5A6","#0C8346", "#C490D1", "purple4"),
                             break.time.by = 1,
                             xlab = "Time (days)")


Fig.kapstar.exp2$plot <- Fig.kapstar.exp2$plot +
  #Add theme to plot
  ggplot2::theme(#text = element_text(family = "A"),
                 axis.title.x = element_text(vjust = -1, size = 10),
                 axis.title.y = element_text(vjust = 3, size = 10),
                 axis.text = element_text(size = 8),
                 panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 axis.line = element_line(colour = "black"),
                 legend.position = "top",
                 legend.text = element_text(size = 10)) +
  #Add p-value
  geom_text(label = "p-value = 0.2", x = 0.75, y = 0.1,
            size = 4, family = "A") +
  #Add animal outline
  add_phylopic(x = 0.75, y = 0.25, img = seastar, alpha = 1, ysize = 0.2) 

Fig.kapstar.exp2

#Cox proportional hazards model
cox.PO.exp2 <- coxph(Surv(time, status)~food.treatment*air.treatment, 
                      data = PO.survivorship)
summary(cox.PO.exp2)
## HR < 1: reduction in hazard (protective)
## HR > 1: increase in hazard
## HR = 1: no effect


#-------------------------------------------------
# Survival analysis for Mytilus spp.

#Compute the survival curve (Kaplan-Meier estimate)
km.fit.mussel.exp2 <- survfit(Surv(time, status) ~ food.treatment + air.treatment, 
                           data = Mussel.survivorship)

print(km.fit.mussel.exp2)

# Summary of survival curves
summary(km.fit.mussel.exp2)
# Access to the sort summary table
summary(km.fit.mussel.exp2)$table

#Access the values returned by survfit()
survival.mussel.exp2 <- data.frame(time = km.fit.mussel.exp2$time,
                                n.risk = km.fit.mussel.exp2$n.risk,
                                n.event = km.fit.mussel.exp2$n.event,
                                n.censor = km.fit.mussel.exp2$n.censor,
                                surv = km.fit.mussel.exp2$surv,
                                upper = km.fit.mussel.exp2$upper,
                                lower = km.fit.mussel.exp2$lower
)
head(survival.mussel.exp2)

##Plot the survival curve (Kaplan Meier curve)
Fig.kapmussel.exp2 <- ggsurvplot(km.fit.mussel.exp2, data = Mussel.survivorship,
                               # conf.int = TRUE, 
                               #pval = TRUE,
                               surv.median.line = "hv",
                               #risk.table = TRUE,
                               legend.labs = c("25°C air/+ mussels","30°C air/+ mussels"),
                               legend.title = "Treatment",
                               legend = c(0.8,0.85),
                               palette = c("#8FD5A6","#0C8346"),
                               break.time.by = 1,
                               xlab = "Time (days)")


Fig.kapmussel.exp2$plot <- Fig.kapmussel.exp2$plot +
  #Add theme to plot
  ggplot2::theme(#text = element_text(family = "A"),
                 axis.title.x = element_text(vjust = -1, size = 10),
                 axis.title.y = element_text(vjust = 3, size = 10),
                 panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 axis.line = element_line(colour = "black"),
                 legend.position = "none",
                 legend.text = element_text(size = 10)) +
  #Add p-value
  geom_text(label = "p-value = 0.9", x = 1, y = 0.075,
            size = 4, family = "A") +
  #Add animal outline
  add_phylopic(x = 1, y = 0.25, img = mussel, alpha = 1, ysize = 0.2)

Fig.kapmussel.exp2


#Cox proportional hazards model
cox.mussel.exp2 <- coxph(Surv(time, status)~air.treatment, 
                       data = Mussel.survivorship)
summary(cox.mussel.exp2)


#---------------------
## Making the MS figure

MS_Surv_exp2 <- Fig.kapstar.exp2$plot + Fig.kapmussel.exp2$plot +
  plot_layout(nrow = 2, heights = c(1, 0.75)) +
  plot_annotation(tag_levels = "a") 

MS_Surv_exp2


#png("MS Figures/Fig.5_Mortality_Exp2.png", width = 8, height = 8, units = "in", res = 600)
#MS_Surv_exp2
#dev.off()
