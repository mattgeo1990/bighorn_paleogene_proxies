#makes an Aridity Index map with study site lpcatopns (Fig 1 of CZ17O soil water manuscript)
#uses v3 of the global AI-PET database

library(ggplot2) 
library(raster)
#library(rgdal)
library(dplyr)
library(usmap)
library(maps)

library(terra)
library(sf)

#####################MAIN SITES###########
#coordinates of sites of interest
name <-c(
         "Mojave CRS station",
         "Jornada CSAND",
         "EGSR",
         "Reynolds Nancy Gulch")

lats <- c(
          35.10613,
          32.515077,
          42.463680, #this is the north gate entrance lat lon
          43.169714)

lons <- c( 
           -115.55074,
           -106.790915, 
           -84.000862,
           -116.710982)

df <- data.frame("name" = name, "lon" = lons, "lat" = lats)
coords <- data.frame(x=lons,y=lats)

#specific sites to extract
name_v2 <- c("CRS", "JTree", "PJ", "NancyGulch", "CSAND", "ESGR")
lats_v2 <- c(35.10613,35.169108, 35.195799, 43.169714, 32.515077, 42.45703949126451)
lons_v2 <- c(-115.55074,-115.480113,  -115.347572, -116.710982, -106.790915, -83.9978986477466)
df_specificsites <- data.frame("name" = name_v2, "lon" = lons_v2, "lat" = lats_v2)
coords_v2 <- data.frame(x = lons_v2, y = lats_v2)

#################
#AI is from Zomer and Trabucco 2022, PET calculated from World Clim Data
#citation:
  #Zomer, R. J., Trabucco, A. 2022. Version 3 of the “Global Aridity Index and Potential Evapotranspiration (ET0) Database”: Estimation of Penman-Monteith Reference Evapotranspiration. 
  #Available online from the CGIAR-CSI GeoPortal at: https://cgiarcsi.community/2019/01/24/global-aridity-index-and-potential-evapotranspiration-climate-database-v3/

#local directory where file is stored. Download TIF file from Zomer et al. in order to run code and update the directory
dir_v3 <- "~/Library/CloudStorage/OneDrive-IndianaUniversity/Research/TripleO/soils/WorldClimData/Global-AI_ET0_v3_annual/"
setwd(dir_v3)

AI <- raster("ai_v3_yr.tif")

# load United States state map data from maps library
MainStates <- map_data("state")

#turns MainStates into a raster
spg <- MainStates
coordinates(spg)<- ~long+lat
gridded(spg)
rasterDF<-raster(spg)

AI<- crop(AI[[1]], extent(rasterDF))

#convert to data frame for use in ggplot
AI_df <- as.data.frame(AI, xy = TRUE)
AI_df <- AI_df[complete.cases(AI_df), ] #get rid rows with NA values
AI_df$awi_pm_sr_yr = AI_df$awi_pm_sr_yr/10000 #convert units in the input file

#customize colors
cuts=c(0,0.05,0.2,0.5,0.65,0.8,1,1.25,10) #set breaks

#make another variable with the breaks
AI_df$awi_pm_sr_yr_2 = cut(AI_df$awi_pm_sr_yr, breaks = cuts)
my_col = c("#C21500", "#E37800", "#E3AF00","#F6F205", "#90F605", "#6AB505", "#05B56A","#057AB5", "#0325D1")


AIplot <- ggplot()+
#commented out for faster plotting while getting points in
  geom_raster(data = AI_df, aes(x=x, y=y, fill = awi_pm_sr_yr_2))+
  scale_fill_manual(values = my_col, name = "Aridity Index") + 
  
  #state outlines
 geom_polygon( data=MainStates, aes(x=long, y=lat, group=group),
                color="black", fill=NA )+
  
  #Main study sites
  geom_point(data = coords, aes(x = x, y =y), shape = 21, color = "black", fill = "white", size = 4)+
  
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ylab("Latitude (°)")+
  xlab("Longitude (°)")
#  coord_map()
AIplot

path.to.figs = "~/Documents/Soil Water 17O Manuscript/code_figures/MS_refined_figs_v2"
ggsave(filename = "AImap_v3.pdf", plot = AIplot, device = cairo_pdf, width = 6, height = 3.2, path = path.to.figs)


#########AI values at the points of interest##########
AIvalues <- extract(AI,coords)/10000
df<-cbind.data.frame(df,AIvalues)

AIvalues_v2 <- extract(AI,coords_v2)/10000
df_specificsites<-cbind.data.frame(df_specificsites,AIvalues_v2)
