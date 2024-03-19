######################################################################################
# Title: Creating a site map for Pisaster survey & collection locations in Bamfield
#
# Author: Lydia Walton
# Last updated: 19-Mar-2024
######################################################################################

dev.off() #Clear open plots
rm(list = ls(all=T)) #Clear environment

#Set working directory
library(here)
# Libraries needed to run this code
library(raster)
library(maps) 
library(mapdata)
library(maptools)
library(rgeos)
library(rgdal)
library(ggplot2)
library(ggsn)
library(tidyverse)
library(PBSmapping)
library(cowplot)


###Create theme for plotting
#Font change to Times New Roman as "A"
#windowsFonts(A = windowsFont("Times New Roman"))

#Lydia's theme
LW_theme <- theme_classic() +
  theme(#text = element_text(family = "A"),
        #axis.title.x = element_text(vjust = -1, size = 10),
        #axis.title.y = element_text(vjust = 2, size = 10),
        #strip.text = element_text(size = 10),
        panel.background = element_rect(fill = NA, color = "black"))



#################
# Van Isle Map
#################

#Use this data for the Van Isle map (includes USA border)
data(nepacLLhigh)

# crop map to show the west coast of Canada
VanIsle.Map <- ggplot() +
  geom_polygon(data = nepacLLhigh, aes(x = X, y = Y, group = PID),
               col = "black", fill = "grey75", lwd = 0.01) +
  coord_map(xlim = c(-129, -122), ylim = c(47.5, 51.5)) +
  LW_theme +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  
  #labs(x = "Longitude", y = "Latitude") +
  #Add text labels
  annotate("text", x = -123, y = 51, label = "Canada", 
           size = 4, fontface = "bold") +
  annotate("text", x = -126.3, y = 50.1, label = "Vancouver \nIsland", 
           size = 3, fontface = "bold") +
  
  #Add box where Barkley Sound is
  annotate("rect", xmin = -125.005, xmax = -125.2665, ymin = 48.8, ymax = 48.9, 
           col = "red", fill = NA, size = 0.8) 


VanIsle.Map



########################
# PO Site Map 
########################

#Read in sites
POsites <- read.csv("MS Data/Summer2023_SiteMap.csv")

#Hakai basemap
hakai <- readOGR("MS Data/COAST_TEST2.shp")

### chose the lat/long extent you want to show
BamfieldSites <- extent(-125.17, -125.1, 48.8, 48.855)

### crop your shapefile polygons to the extent defined
# takes a moment to run 
hakai.cropped <- crop(hakai,BamfieldSites)

### project and fortify (i.e. turn into a dataframe)
hakai.cropped.df <- fortify(hakai.cropped)


#Map of PO collection sites 

SiteMap <- ggplot() + 
  LW_theme +
  geom_polygon(data= hakai.cropped.df, aes(x=long,y=lat,group= group),
               colour= "black", size=0.1, fill='grey75')+
  coord_cartesian(xlim = c(-125.1665, -125.105), ylim=c(48.81, 48.853)) +
  labs(x = "Longitude", y = "Latitude") +
  
  #Add site locations
  geom_point(data=POsites, aes(x=Lon, y=Lat, color = Site, shape = Site, fill = Site), 
             size=2, stroke=1.5) +
  scale_colour_manual(values = c("black","#DE369D", "#DE369D", "#DE369D")) +
  scale_fill_manual(values = c("black","#DE369D", "#DE369D", "#DE369D")) +
  scale_shape_manual(values = c(25,19,19,19)) +
  theme(legend.position = "none") +
  
  #Add scale bar
  scalebar(data = hakai.cropped.df, location = "bottomright",
           dist = 0.5, transform = TRUE, dist_unit = "km",
           model = "WGS84", height = 0.02, st.size = 3, st.dist = 0.02,
           st.bottom = FALSE, border.size = 0.5, anchor = c(x = -125.105, y = 48.81)) +
  
  #Add north arrow
  north(data = hakai.cropped.df, location = "topright",
        scale = 0.1, symbol = 1, anchor= c(x = -125.102, y = 48.854)) 

SiteMap <- SiteMap +
  #Add site labels
  geom_label(aes(label= "BMSC", x=-125.1352, y= 48.8335), size=3, family = "A") +
  geom_label(aes(label= "EB", x=-125.1431, y= 48.8385), size=3, family = "A") +
  geom_label(aes(label= "GN", x=-125.1109, y= 48.8365), size=3, family = "A") +
  geom_label(aes(label= "SP", x=-125.1296, y= 48.8345), size=3, family = "A") +
  geom_label(aes(label= "Bamfield Inlet", x=-125.1385, y= 48.8265), size=3) +
  geom_text(aes(label= "Bamfield", x=-125.12, y= 48.82), size=4, family = "A", fontface = "bold") 

SiteMap


###################################
# MS MAP (Combined with inset)
####################################

#Combine the plots, adding VI map as an inset in the top left corner
(map_w_inset <- ggdraw(SiteMap) +
   draw_plot(VanIsle.Map, x = -0.047, y = 0.603, width = 0.58, height = 0.38))


####Export file as a .png file

#png(filename = here("", "Fig.1_SiteMap_Summer2023.png"), width = 7.5, height = 5.5, 
#    units = "in", pointsize = 15, res = 600)

#(map_w_inset <- ggdraw(SiteMap) +
#    draw_plot(VanIsle.Map, x = -0.047, y = 0.603, width = 0.58, height = 0.38))

#dev.off()




