#this script compares the d18O-D17O values generated in our study with other types of obs'd waters

#other data comes from Aron et al. 2021 Chem Geo
#with the addition of soil water data from Beverly et al. 2021 and Alexandre et al. 2025

library(ggplot2)  
library(ggpubr)
library(scales)
library(dplyr)
library(forcats)
library(lubridate)
library(readxl)
library(plotly)
library(stringr)
library(readr) #for read_delim

######## SET UP AND LOAD FILES######

setwd("~/Documents/Soil Water 17O Manuscript/data")
path.to.figs <- "~/Documents/Soil Water 17O Manuscript/code_figures/MS_draft_figs/"
path.to.figs.refined <- "~/Documents/Soil Water 17O Manuscript/code_figures/MS_refined_figs_v2/"

sw <- read.csv(file = "sw.csv", header = T)
mw <- read.csv(file = "mw.csv", header = T)

####organize by sites, set up for figures####
MOJsw <- sw %>% filter(siteID.1 == "MOJ")
JORsw <- sw %>% filter(siteID.1 == "JOR" & siteID.2 == "CSAND") #remove Red lake playa from diff rlated study
REYsw <- sw %>% filter(siteID.1 == "REY")
ESGRsw <- sw %>% filter(siteID.1 == "ESGR")

MOJmetw <- mw %>% filter(siteID.1 == "MOJ")
JORmetw <- mw %>% filter(siteID.1 == "JOR") 
REYmetw <- mw %>% filter(siteID.1 == "REY")
ESGRmetw <- mw %>% filter(siteID.1 == "ESGR") %>% filter(IPL.num.simple != 2018) #removing one sample that is too evaporated to be reasonable and had 75% headspace
#as of nov 2024, this now includes data from Kelson ESGR MS too

#colors from AGU poster
MOJcolor <- "darkorange4"
MOJcolor.minor <- "darkorange"
MOJcolor.minor2 <- "darkorange3"
REYcolor <- "darkolivegreen"
REYcolor.minor <- "darkolivegreen3"
JORcolor <- "#506a98" #this is a dark blue
JORcolor.minor <- "#99b8de" #this is a light blue
ESGRcolor <- "violetred4"
ESGRcolor.minor <- "violetred1"

site.1.colors <- c( "JOR" = JORcolor,
                    "MOJ" = MOJcolor,
                    "REY" = REYcolor,
                    "ESGR" = ESGRcolor)

site.1.colors.minor <- c( "JOR" = JORcolor.minor,
                          "MOJ" = MOJcolor.minor,
                          "REY" = REYcolor.minor,
                          "ESGR" = ESGRcolor.minor)

MStheme <- theme_bw()+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(axis.text.x=element_text(angle=60, hjust=1), axis.title = element_text(size = 12),axis.text = element_text(size = 12))

MStheme_isos <- theme_bw()+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(axis.title = element_text(size = 12),axis.text = element_text(size = 12))

##########BA modeling data######

#soil evaporation model runs from Kelson et al. 2023 https://doi.org/10.1016/j.gca.2023.06.034
#model run used in figure in manuscript is included in github repository, all else commented out
BAdf <- read.csv("~/Documents/Soil Water 17O Manuscript/data/BarnesAllisonModeling_fromKelson2023/BA1983_baserun.csv")
#BAdrydf <- read.csv("~/Documents/Soil Water 17O Manuscript/data/BarnesAllisonModeling_fromKelson2023/BA1983_dryrun.csv")
#BAhumiddf <- read.csv("~/Documents/Soil Water 17O Manuscript/data/BarnesAllisonModeling_fromKelson2023/BA1983_humidrun.csv")

BAdf$d.excess <- BAdf$dDsw -8 *BAdf$d18Osw
#BAdrydf$d.excess <- BAdrydf$dDsw -8 *BAdrydf$d18Osw
#BAhumiddf$d.excess <- BAhumiddf$dDsw -8 *BAhumiddf$d18Osw

###BA trends
BAdf$dp18Osw <- log(BAdf$d18Osw/1000+1)*1000
BAdf$dp17Osw <- (BAdf$D17Osw/1000)+0.528*BAdf$dp18Osw
SWL_BAdf <- lm(BAdf$d17Osw ~ BAdf$d18Osw)
summary(SWL_BAdf)
lambda_empirical_BAdf <- lm(BAdf$dp17Osw~BAdf$dp18Osw)
summary(lambda_empirical_BAdf)

#BAdrydf$dp18Osw <- log(BAdrydf$d18Osw/1000+1)*1000
#BAdrydf$dp17Osw <- (BAdrydf$D17Osw/1000)+0.528*BAdrydf$dp18Osw
#SWL_BAdrydf <- lm(BAdrydf$d17Osw ~ BAdrydf$d18Osw)
#summary(SWL_BAdrydf)
#lambda_empirical_BAdrydf <- lm(BAdrydf$dp17Osw~BAdrydf$dp18Osw)
#summary(lambda_empirical_BAdrydf)

#BAhumiddf$dp18Osw <- log(BAhumiddf$d18Osw/1000+1)*1000
#BAhumiddf$dp17Osw <- (BAhumiddf$D17Osw/1000)+0.528*BAhumiddf$dp18Osw
#lambda_empirical_BAhumiddf <- lm(BAhumiddf$dp17Osw~BAhumiddf$dp18Osw)
#summary(lambda_empirical_BAhumiddf)

#model comparison in isotope space
#ggplot()+
#  geom_point(data = BAdf, aes(x = d18Osw, y = dDsw), size = 2, shape = 21, color = "black")+
#  geom_point(data = BAdrydf, aes(x = d18Osw, y = dDsw), size = 2, shape = 21, color = "red")+
#  geom_point(data = BAhumiddf, aes(x = d18Osw, y = dDsw), size = 2, shape = 21, color = "blue")+

#ggplot()+
#  geom_point(data = BAdf, aes(x = D17Osw, y = d.excess), size = 2, shape = 21, color = "black")+
#  geom_point(data = BAdrydf, aes(x = D17Osw, y = d.excess), size = 2, shape = 21, color = "red")+
#  geom_point(data = BAhumiddf, aes(x = D17Osw, y = d.excess), size = 2, shape = 21, color = "blue")+
#  geom_point(data  = sw, aes(x = D17O_pmg, y = d_excess.mean), size = 2, shape = 22, fill = "black")+
#  geom_point(data  = mw, aes(x = D17O_pmg, y = d_excess.mean), size = 2, shape = 22, fill = "purple")

#########other data sources############

#run this first
#can be downloaded from Aron et al. 2021 chem geo https://doi.org/10.1016/j.chemgeo.2020.120026
#/Users/juliakelson/Library/CloudStorage/OneDrive-IndianaUniversity/Research/TripleO/waterDB_FromPAron/published_17O_water_database.r

#also add Beverly 2021 data, https://doi.org/10.1016/j.epsl.2021.116952
#Table S13: δ18O, δD, and Δ17O water data. Values in ‰VSMOW-SLAP unless otherwise indicated.
Beverly <- read_excel("~/Documents/Soil Water 17O Manuscript/data/Beverly2021_waterisos_supp_forR.xlsx")
Beverly.sw <- Beverly %>% filter(Type == "Soil")

#add Alexandre et al. 2025 data
#https://doi.pangaea.de/10.1594/PANGAEA.984133
Alexandre.sw <-read_delim("~/Documents/Soil Water 17O Manuscript/data/AlexandreA-etal_2025/datasets/6_Soil_water_isotopes.tab", skip = 42,delim = "\t", show_col_types = FALSE)
Alexandre.mw <-read_delim("~/Documents/Soil Water 17O Manuscript/data/AlexandreA-etal_2025/datasets/5_Precipitation_isotopes.tab", skip = 33,delim = "\t", show_col_types = FALSE)
Alexandre.pw <-read_delim("~/Documents/Soil Water 17O Manuscript/data/AlexandreA-etal_2025/datasets/7_Plant_water_isotopes.tab", skip = 61,delim = "\t", show_col_types = FALSE)



######## STATISTIC ON d-excess, D17O correlation for soil water and precipitation, this study only

d18OD17Olm<- lm(sw$D17O_pmg~sw$d18O.mean)
summary(d18OD17Olm)

precip_only <- mw[mw$Water.Type == "Precipitation",]
d18OD17Olm_precip<- lm(precip_only$D17O_pmg~precip_only$d18O.mean)
summary(d18OD17Olm_precip)

DD_sw<- lm(sw$D17O_pmg~sw$d_excess.mean)
summary(DD)

DD_precip <-lm(precip_only$D17O_pmg~precip_only$d_excess.mean)
summary(DD_precip)
########FIGURES #######

d18O_D17O <- ggplot() + 
 geom_point(data=plant.water, aes(x=dp18O,y=(dp17O-0.528*dp18O)*1000), color="darkgreen",size=1.5, shape = 17) + #triangle
  geom_point(data=precipitation, aes(x=dp18O,y=(dp17O-0.528*dp18O)*1000), color="gray",size=1.5, shape = 15) + #rectanlge
  geom_point(data=surface.water, aes(x=dp18O,y=(dp17O-0.528*dp18O)*1000), color="firebrick",size=1.5, shape = 15) + #rectangle
  
  #data from Beverly 2021
  geom_point(aes(x = Beverly.sw$d18O, y = Beverly.sw$avg_D17O_pmg), color = "gold", size = 1.5)+
  
  #data from Alexandre 2025
  geom_point(aes(x = (exp(Alexandre.sw$`δ'18O soil w [‰ SMOW]`/1000)-1)*1000 , y =  Alexandre.sw$`17O xs soil w [per meg]` ), 
             color = "gold", size = 1.5)+
  geom_point(aes(x = (exp(Alexandre.pw$`δ'18O plant water [‰ SMOW] (Observed)`/1000)-1)*1000 , y =  Alexandre.pw$`17O xs plant water [per meg] (Observed)`), 
             color = "darkgreen", size = 1.5,shape = 17)+
  #excluding plant water data because there is no d-excess data reported
  
  #need to unprime the dp18O
  
  #data from this study
  geom_point(aes(x = mw$d18O.mean[mw$Water.Type == "Precipitation"], y= mw$D17O_pmg[mw$Water.Type == "Precipitation"]), fill = "gray",color = "black", shape = 22, size = 3)+
  geom_point(aes(x = mw$d18O.mean[mw$Water.Type != "Precipitation"], y= mw$D17O_pmg[mw$Water.Type != "Precipitation"]), fill = "firebrick",color = "black", shape = 22, size = 3)+
  geom_point(aes(x = sw$d18O.mean, y= sw$D17O_pmg, fill = sw$siteID.1 ), shape = 21, size = 3)+
  scale_fill_manual(name=" ", values = site.1.colors.minor)+

  #scale_x_continuous(limits = c(-76,35), expand = c(0, 0), breaks=seq(-75,25,25)) +
  #scale_y_continuous(limits = c(-290,120), expand = c(0, 0), breaks=seq(-200,100,100)) +
  
  geom_line(data = BAdf[BAdf$z>BAdf$z_ef[1],], aes(x = d18Osw, y = D17Osw), size = 1, color = "black")+
  geom_smooth(data = mw[mw$Water.Type == "Precipitation",], aes(x = d18O.mean, y = D17O_pmg), method = 'lm', se = TRUE, color = "black", linetype = "dotdash", size = 0.5)+
  geom_smooth(data = sw, aes(x = d18O.mean, y = D17O_pmg), method = 'lm', se = TRUE, color = "black", linetype = 2, size = 0.5)+
  
  #geom_abline(slope = d18OD17Olm$coefficients[2], intercept = d18OD17Olm$coefficients[1])+
  theme(legend.position ="right") + 
  labs(x=expression(delta*"'"^"18"*"O (\u2030, VSMOW-SLAP)"), y=expression(Delta*"'"^"17"*"O (per meg, VSMOW-SLAP)")) +
  #ggtitle(expression(bold(delta*"'"^"18"*"O vs. "*Delta*"'"^"17"*"O")))  +
  scale_x_continuous(limits = c(-25,25), expand = c(0, 0), breaks=seq(-20,20,10)) +
  scale_y_continuous(limits = c(-205,110), expand = c(0,0), breaks = seq(-200, 100, 50))  +
  MStheme_isos

plot(d18O_D17O)  

#only plots some of the water types from Aron et al. 2021
dxs_D17O_2 <- ggplot()+
 geom_point(data = plant.water, aes(x=d.excess,y=(dp17O-0.528*dp18O)*1000), color="darkgreen", size=1.5, shape = 17) + 
 geom_point(data = precipitation, aes(x=d.excess,y=(dp17O-0.528*dp18O)*1000) , color="gray", size=1.5, shape = 15) + 
 geom_point(data = surface.water, aes(x=d.excess,y=(dp17O-0.528*dp18O)*1000), color="firebrick",size=1.5, shape = 15) + 
  
 labs(x="d-excess (\u2030, VSMOW-SLAP)", y=expression(Delta*"'"^"17"*"O (per meg, VSMOW-SLAP)")) +
 # ggtitle(expression(bold("d-excess vs. "*Delta*"'"^"17"*"O"))) + 
  #scale_color_manual(name=" ", values = broad.colors) +
  
  #data from Beverly
  geom_point(aes(x = Beverly.sw$dexcess, y = Beverly.sw$avg_D17O_pmg), color = "gold", size = 1.5)+
  
  #data from Alexandre
  geom_point(aes(x = Alexandre.sw$`d xs soil w [‰]`, y = Alexandre.sw$`17O xs soil w [per meg]`), color = "gold", size = 1.5)+
  #geom_point(aes(x = (exp(Alexandre.pw$ , y =  Alexandre.pw$`17O xs plant water [per meg] (Observed)`), 
  #           color = "darkgreen", size = 1.5,shape = 17)+
  #no dexcess data for the plant waters!
  
  #data from this study
  geom_point(aes(x = mw$d_excess.mean[mw$Water.Type == "Precipitation"], y= mw$D17O_pmg[mw$Water.Type == "Precipitation"]), fill = "gray",shape = 22, size = 3)+
  geom_point(aes(x = mw$d_excess.mean[mw$Water.Type != "Precipitation"], y= mw$D17O_pmg[mw$Water.Type != "Precipitation"]), fill = "firebrick",shape = 22, size = 3)+
  geom_point(aes(x = sw$d_excess.mean, y= sw$D17O_pmg, fill = sw$siteID.1 ),shape = 21, size = 3)+
  scale_fill_manual(name=" ", values = site.1.colors.minor)+
  
  geom_line(data = BAdf[BAdf$z>BAdf$z_ef[1],], aes(x = d.excess, y = D17Osw), size = 1,color = "black")+
  geom_smooth(data = mw[mw$Water.Type == "Precipitation",], aes(x = d_excess.mean, y = D17O_pmg), method = 'lm', se = TRUE, color = "black", linetype = "dotdash", size = 0.5)+
  geom_smooth(data = sw, aes(x = d_excess.mean, y = D17O_pmg), method = 'lm', se = TRUE, color = "black", linetype = 2, size = 0.5)+
  
  
  #geom_abline(intercept = DD$coefficients[1], slope = DD$coefficients[2])+
  scale_x_continuous(limits = c(-155,40), expand = c(0, 0), breaks=seq(-150,30,30)) +
  scale_y_continuous(limits = c(-205,110), expand = c(0,0), breaks = seq(-200, 100, 50))  +
  MStheme_isos

plot(dxs_D17O_2)  

isoscompare <- ggarrange(d18O_D17O, dxs_D17O_2, nrow =1, ncol =2, common.legend = TRUE)
isoscompare
path.to.figs.postR1 <- "~/Documents/Soil Water 17O Manuscript/code_figures/MS_figures_postR1/" 

ggsave(filename = "isos_compare_wAlexandre.pdf", plot = isoscompare, device = cairo_pdf, height=6,width=10.75, path = path.to.figs.postR1)


