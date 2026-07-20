#code started Dec 31 2024
#time dependent plots for the manuscript
#p1: d18O vs time, sw and mw
#p2: d-excess vs time, sw and mw
#p3: D17O vs time, sw and mw
#p4: soil moisture
#p5: precipitation (in mm)

#modified from Drever AGU plots

library(ggplot2) 
library(ggpubr)
library(scales)
library(dplyr)
library(forcats)
library(lubridate)
library(readxl)
library(plotly)
library(stringr)

path.to.figs <- "~/Documents/Soil Water 17O Manuscript/code_figures/MS_draft_figs/"
path.to.figs.refined <- "~/Documents/Soil Water 17O Manuscript/code_figures/MS_refined_figs_v2/"

#####read in the isotope data########

#combined Picarro and IRMS data, made in CZ17O_dataviz_waters_forMS.R
mw <-read.csv(file = "~/Documents/Soil Water 17O Manuscript/data/mw.csv") 
sw <-read.csv(file = "~/Documents/Soil Water 17O Manuscript/data/sw.csv") 

mw$date.collection <- as_date(mw$date.collection)
sw$date.collection <- as_date(sw$date.collection)


###### MeanWeightedP values entered from spreadsheet, updated picarro data 12/31/24 #####
mean.weighted.precip.MOJ.d18O <- -7.315 #Weighted_Precip_SMichigan.R in weighted precip code in sw manuscript folder, Dec 31 2024
mean.weighted.precip.MOJ.d2H <- -57.078
mean.weighted.precip.MOJ.dexcess <- mean.weighted.precip.MOJ.d2H - mean.weighted.precip.MOJ.d18O*8

mean.weighted.precip.REY.d18O <- -11.99 #Weighted_Precip_Reynolds.R in weighted precip code in sw manuscript folder, Dec 30 2024
mean.weighted.precip.REY.d2H <- -90.14
mean.weighted.precip.REY.dexcess <- mean.weighted.precip.REY.d2H - mean.weighted.precip.REY.d18O*8

mean.weighted.precip.JOR.d18O <- -5.075 #Weighted_Precip_Reynolds.R in weighted precip code in sw manuscript folder, Dec 30 2024
mean.weighted.precip.JOR.d2H <- -40.78
mean.weighted.precip.JOR.dexcess <- mean.weighted.precip.JOR.d2H - mean.weighted.precip.JOR.d18O*8

mean.weighted.precip.ESGR.d18O <- -8.656341 #Weighted_Precip_SMichigan.R in weighted precip code in sw manuscript folder, Dec 30 2024
mean.weighted.precip.ESGR.d2H <- -54.69655
mean.weighted.precip.ESGR.dexcess <- mean.weighted.precip.ESGR.d2H - mean.weighted.precip.ESGR.d18O*8
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

MStheme_time <- theme_bw()+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 12)) #only formats the y axis so the x axis label can appear only on bottom panel?
  
######### REYNOLDS TIME PLOTS ##########
#load in precipitation data from RCEW, averaged to be daily in the code Drever_Reynolds_Elena_3_18_24.R
REYppt <- read.csv("~/Documents/Soil Water 17O Manuscript/data/climo data/Reynolds climo data/RCEW_098c_daily_ppta.csv", sep = ",")
REYdailyppt_studyperiod_only <- REYppt %>% filter(datetime > "2021-08-01")
REYdailyppt_studyperiod_only$datetime <- as_date(REYdailyppt_studyperiod_only$datetime) #not sure why this is needed

#load in SM data, modified from Drever_Reynolds_Elena_3_18_24.R
path.to.data = "~/Documents/Soil Water 17O Manuscript/data/climo data/Reynolds climo data/Reynolds Huber Complete Climate Data/098"

#SM_2015 <- read.csv(paste(path.to.data,"q098fhydraprobe2015.csv", sep ="/"), header = T) #reading in all the soil moisture data
#SM_2016 <- read.csv(paste(path.to.data,"q098fhydraprobe2016.csv", sep ="/"), header = T)
#SM_2017 <- read.csv(paste(path.to.data,"q098fhydraprobe2017.csv", sep ="/"), header = T)
#SM_2018 <- read.csv(paste(path.to.data,"q098fhydraprobe2018.csv", sep ="/"), header = T)
#SM_2019 <- read.csv(paste(path.to.data,"q098fhydraprobe2019.csv", sep ="/"), header = T)
#SM_2020 <- read.csv(paste(path.to.data,"q098fhydraprobe2020.csv", sep ="/"), header = T)
SM_2021 <- read.csv(paste(path.to.data,"q098fhydraprobe2021.csv", sep ="/"), header = T)
SM_2022 <- read.csv(paste(path.to.data,"q098fhydraprobe2022.csv", sep ="/"), header = T)
SM_2023 <- read.csv(paste(path.to.data,"q098fhydraprobe2023.csv", sep ="/"), header = T)
SM_2024 <- read.csv(paste(path.to.data,"q098fhydraprobe2024.csv", sep ="/"), header = T)

#SM <- rbind(SM_2015, SM_2016, SM_2017, SM_2018, SM_2019, SM_2020, SM_2021, SM_2022, SM_2023, SM_2024) #combining data sets from different years
SM <- rbind(SM_2021, SM_2022, SM_2023, SM_2024) #combining data sets from different years
SM$datetime <- as_date(ymd_hm(SM$datetime)) #lubridate (ymd_hm) parses as a POSISCt format (not sure why), as_date makes it a date
#SM <- filter(SM, year(datetime) >= 2021 & year(datetime) <= 2024)
SM[SM == -999] <- NA #getting rid of the -999 values so they won't mess up calculations later on
SM[SM == 0] <- NA #getting rid of the two 0 values so they won't mess up calculations later on
SM$wat030i[SM$wat030i == 0.263] <- NA #getting rid obviously bad values 

#might use
REYprecip <- REYmetw %>% filter(Water.Type == "Precipitation")
REYprecip <- arrange(REYprecip, date.collection)
REYswdeep <- REYsw %>% filter(bottom.depth > 52)
REYswsurf <- REYsw %>% filter(bottom.depth < 18)

ptimeP <- ggplot()+
  geom_line(data = REYdailyppt_studyperiod_only, aes(x = datetime, y = ppta), color = REYcolor)+
  xlab("")+
  ylab("Daily Precipitation (mm)")+
  scale_y_continuous(limits = c(0, 45), expand = c(0,0))+
  scale_x_date(date_breaks = "2 month", limits = as.Date(c("2021-08-01", "2023-09-30")), expand = c(0,0), labels = NULL)+
  #rectangles for DJF and JJA
  geom_rect(aes(xmin =as.Date("2021-12-01"), xmax = as.Date("2022-02-28"), ymin = -Inf, ymax = Inf),fill = "lightblue",  alpha = 0.6) +
  geom_rect(aes(xmin =as.Date("2022-12-01"), xmax = as.Date("2023-02-28"), ymin = -Inf, ymax = Inf),fill = "lightblue", alpha = 0.6) +
  geom_rect(aes(xmin =as.Date("2022-06-01"), xmax = as.Date("2022-08-30"), ymin = -Inf, ymax = Inf),fill = "coral1", alpha = 0.6) +
  geom_rect(aes(xmin =as.Date("2023-06-01"), xmax = as.Date("2023-08-30"), ymin = -Inf, ymax = Inf),fill = "coral1", alpha = 0.6) +
  MStheme
ptimeP

#ggplotly(ptimeP)
ptimeSM <- ggplot() +
  geom_line(data = SM, aes(x = as_date(datetime), y = as.numeric(wat005i)), color = "black") + #interplant 5 cm
  geom_line(data = SM, aes(x = as_date(datetime), y = as.numeric(wat030i)), color = "gray60") + #interplant 30 cm
  geom_line(data = SM, aes(x = as_date(datetime), y = as.numeric(wat063i)), color = "gray90") + #interplant 63 cm
  xlab("")+
  ylab("Soil Moisture (-)")+
  scale_x_date(date_breaks = "2 month", limits = as.Date(c("2021-08-01", "2023-09-30")), expand = c(0,0), labels = NULL)+
  #rectangles for DJF and JJA
  geom_rect(aes(xmin =as.Date("2021-12-01"), xmax = as.Date("2022-02-28"), ymin = -Inf, ymax = Inf),fill = "lightblue",  alpha = 0.6) +
  geom_rect(aes(xmin =as.Date("2022-12-01"), xmax = as.Date("2023-02-28"), ymin = -Inf, ymax = Inf),fill = "lightblue", alpha = 0.6) +
  geom_rect(aes(xmin =as.Date("2022-06-01"), xmax = as.Date("2022-08-30"), ymin = -Inf, ymax = Inf),fill = "coral1", alpha = 0.6) +
  geom_rect(aes(xmin =as.Date("2023-06-01"), xmax = as.Date("2023-08-30"), ymin = -Inf, ymax = Inf),fill = "coral1", alpha = 0.6) +
  MStheme
ptimeSM

#ggplotly(ptimeSM)
ptime1 <- ggplot()+
  geom_line(data = REYprecip, aes(x = date.collection, y = d18O.mean), size = 0.5, color = REYcolor)+
  geom_point(data = REYprecip, aes(x = date.collection, y = d18O.mean, shape = Water.Type), size = 3, color = REYcolor, fill = REYcolor, shape = 22)+
  geom_line(data = REYswdeep, aes(x = date.collection, y = d18O.mean), size = 0.5, color = "gray90")+
  geom_line(data = REYswsurf, aes(x = date.collection, y = d18O.mean), size = 0.5, color = "black")+
  geom_point(data = REYsw, aes(x = date.collection, y = d18O.mean, color = bottom.depth), size = 3)+
  scale_color_gradient(high = "grey90", low = "black")+
  # geom_line(data = REYswsurf, aes(x = date.collection, y = d18O.mean), size = 1, color = REYcolor.minor)+#doesnt look good
  xlab("")+
  ylab(expression(paste(delta^{18}, "O (‰, SMOW)")))+
  scale_x_date(date_breaks = "2 month", limits = as.Date(c("2021-08-01", "2023-09-30")), expand = c(0,0), labels = NULL)+
  geom_rect(aes(xmin =as.Date("2021-12-01"), xmax = as.Date("2022-02-28"), ymin = -Inf, ymax = Inf),fill = "lightblue",  alpha = 0.6) +
  geom_rect(aes(xmin =as.Date("2022-12-01"), xmax = as.Date("2023-02-28"), ymin = -Inf, ymax = Inf),fill = "lightblue", alpha = 0.6) +
  geom_rect(aes(xmin =as.Date("2022-06-01"), xmax = as.Date("2022-08-30"), ymin = -Inf, ymax = Inf),fill = "coral1", alpha = 0.6) +
  geom_rect(aes(xmin =as.Date("2023-06-01"), xmax = as.Date("2023-08-30"), ymin = -Inf, ymax = Inf),fill = "coral1", alpha = 0.6) +
  
  MStheme
ptime1
#the lines aren't connecting because sample 2000 doesn't have [good] Picarro data

ptime2 <-  ggplot()+
  geom_point(data = REYprecip, aes(x = date.collection, y = d_excess.mean), size = 3, color = REYcolor, fill = REYcolor, shape = 22)+
  geom_line(data = REYprecip, aes(x = date.collection, y = d_excess.mean), size = 0.5, color = REYcolor)+
  geom_line(data = REYswdeep, aes(x = date.collection, y = d_excess.mean), size = 0.5, color = "gray90")+
  geom_line(data = REYswsurf, aes(x = date.collection, y = d_excess.mean), size = 0.5, color = "black")+
  geom_point(data = REYsw, aes(x = date.collection, y = d_excess.mean, color = bottom.depth), size = 3)+
  scale_color_gradient(high = "grey90", low = "black")+
  ylab(expression(paste(italic("d"),"-excess (‰, SMOW)"))) +
  xlab("")+
  scale_x_date(date_breaks = "2 month", limits = as.Date(c("2021-08-01", "2023-09-30")), expand = c(0,0), labels = NULL)+
  geom_rect(aes(xmin =as.Date("2021-12-01"), xmax = as.Date("2022-02-28"), ymin = -Inf, ymax = Inf),fill = "lightblue",  alpha = 0.6) +
  geom_rect(aes(xmin =as.Date("2022-12-01"), xmax = as.Date("2023-02-28"), ymin = -Inf, ymax = Inf),fill = "lightblue", alpha = 0.6) +
  geom_rect(aes(xmin =as.Date("2022-06-01"), xmax = as.Date("2022-08-30"), ymin = -Inf, ymax = Inf),fill = "coral1", alpha = 0.6) +
  geom_rect(aes(xmin =as.Date("2023-06-01"), xmax = as.Date("2023-08-30"), ymin = -Inf, ymax = Inf),fill = "coral1", alpha = 0.6) +
  
  MStheme 
  
#the lines aren't connecting because sample 2000 doesn't have [good] Picarro data
ptime2

ptime3 <- ggplot()+
  geom_line(data = REYprecip, aes(x = date.collection, y = D17O_pmg), size = 0.5, color = REYcolor)+
  geom_line(data = REYswdeep, aes(x = date.collection, y = D17O_pmg), size = 0.5, color = "gray90")+
  geom_line(data = REYswsurf, aes(x = date.collection, y = D17O_pmg), size = 0.5, color = "black")+
  
  geom_errorbar(data  = REYprecip, aes(x = date.collection, y = D17O_pmg, ymin = D17O_pmg - D17O_err, ymax = D17O_pmg + D17O_err), width = 0, color = REYcolor)+
  geom_errorbar(data  = REYsw, aes(x = date.collection, y = D17O_pmg, ymin = D17O_pmg - D17O_err, ymax = D17O_pmg + D17O_err), width = 0, color = "black")+
  
  geom_point(data = REYprecip, aes(x = date.collection, y = D17O_pmg), size = 3, fill = REYcolor, color = REYcolor, shape = 22)+
  geom_point(data = REYsw, aes(x = date.collection, y = D17O_pmg, color = bottom.depth), size = 3)+
  
  scale_color_gradient(high = "grey90", low = "black")+
  
  
  ylab(expression(Delta^"'17"*"O (pmg, SMOW-SLAP)"))+
  xlab("")+
  scale_x_date(date_breaks = "2 month" , date_labels = "%m/%y", limits = as.Date(c("2021-08-01", "2023-09-30")), expand = c(0,0))+
  
  #rectangles for DJF and JJA
  geom_rect(aes(xmin =as.Date("2021-12-01"), xmax = as.Date("2022-02-28"), ymin = -Inf, ymax = Inf),fill = "lightblue",  alpha = 0.6) +
  geom_rect(aes(xmin =as.Date("2022-12-01"), xmax = as.Date("2023-02-28"), ymin = -Inf, ymax = Inf),fill = "lightblue", alpha = 0.6) +
  geom_rect(aes(xmin =as.Date("2022-06-01"), xmax = as.Date("2022-08-30"), ymin = -Inf, ymax = Inf),fill = "coral1", alpha = 0.6) +
  geom_rect(aes(xmin =as.Date("2023-06-01"), xmax = as.Date("2023-08-30"), ymin = -Inf, ymax = Inf),fill = "coral1", alpha = 0.6) +
  MStheme 
  
ptime3

ptime <- ggarrange(ptimeP,ptimeSM, ptime1, ptime2, ptime3, nrow = 5, ncol = 1, common.legend = TRUE)
ptime

ggsave(filename = "ptime_REY.pdf", plot = ptime, device = cairo_pdf, height=8.25,width=10.75, path = path.to.figs.refined)



######### JORNADA TIME PLOTS########
#left off here Dec 31 2024

JORprecip <- JORmetw %>% filter(Water.Type == "Precipitation")
JORswdeep <- JORsw %>% filter(bottom.depth>55)
JORswsurf <- JORsw %>% filter(bottom.depth<20)

#this has data from UTEP theses and the UMIch updated d18O values
PicarroData <- read_excel("~/Documents/Soil Water 17O Manuscript/code_figures/Mean Weighted Precip Code/input files for Jornada/JornadaUTEP_RainCollection_Summary.xlsx")
#isodf <- PicarroData %>% filter(Site == "JER") %>% filter(Repeat_measurement == 0)

#reading in precip data
path.to.data = "~/Documents/Soil Water 17O Manuscript/data/climo data/Jornada climo data/LTER climate data/"

metdf <- read.csv(paste(path.to.data,"knb-lter-jrn.21-437033.29.csv",sep = ""), header = T) #1 hour summary of met data starting in 2013
metdf$date <- ymd_hms(metdf$Date)

#daily climate data
climodf <- read.csv("/Users/juliakelson/Documents/Soil Water 17O Manuscript/data/climo data/Jornada climo data/LTER climate data/knb-lter-jrn.210437048.40.csv")
climodf$DateFormatted <- ymd(climodf$Date)
climodfshort <- climodf %>% filter(DateFormatted >= "2018-01-01 00:00")  #data starts in 2021
climodfshort <- climodfshort %>% mutate(climodfshort, Date = date(DateFormatted))

#soil moisture data
#this version stops in 2022
#SMdf <- read.csv("/Users/juliakelson/Documents/Soil Water 17O Manuscript/data/climo data/Jornada climo data/LTER climate data/knb-lter-jrn.210437078.29.csv", sep = ",", header = T) #30 minute data
SMdf <- read.csv("/Users/juliakelson/Documents/Soil Water 17O Manuscript/data/climo data/Jornada climo data/LTER climate data/knb-lter-jrn.210437093.44.csv", sep = ",", header = T) #daily
SMdf$date <- as_date(SMdf$Date)
SMdf <- SMdf %>% filter(SMdf$date > as_date("2018-01-01") & SMdf$date < as_date("2024-01-01") ) #shorten

ptimePJOR <- ggplot()+
  geom_line(data = climodfshort, aes(x = Date, y = Ppt_mm_Tot), color = JORcolor)+
  xlab("")+
  ylab("Daily Precipitation (mm)")+
  scale_y_continuous(limits = c(0, 30), expand = c(0,0))+
  scale_x_date(date_breaks = "2 month", limits = as.Date(c("2018-10-01", "2023-10-31")), expand = c(0,0), labels = NULL)+
  #rectangles for DJF and JJA
  geom_rect(aes(xmin =as.Date(c("2018-12-01", "2019-12-01", "2020-12-01","2021-12-01", "2022-12-01")), 
                xmax = as.Date(c("2019-02-28", "2020-02-28", "2021-02-28", "2022-02-28", "2023-02-28")), ymin = -Inf, ymax = Inf),fill = "lightblue",  alpha = 0.6) +
  geom_rect(aes(xmin =as.Date(c("2019-06-01","2020-06-01", "2021-06-01", "2022-06-01", "2023-06-01")), 
                xmax = as.Date(c("2019-08-30","2020-08-30", "2021-08-30", "2022-08-30", "2023-08-30")), ymin = -Inf, ymax = Inf),fill = "coral1", alpha = 0.6) +

  MStheme
ptimePJOR

ggplotly(ptimePJOR)
ptimeSMJOR <- ggplot() +
  geom_line(data = SMdf, aes(x = as_date(date), y = VwcCorr_Avg_202_10cm), color = "gray90") + # 10 cm
  geom_line(data = SMdf, aes(x = as_date(date), y = VwcCorr_Avg_202_30cm), color = "gray60") + # 30 cm
  xlab("")+
  ylab("Soil Moisture (-)")+
  scale_x_date(date_breaks = "2 month", limits = as.Date(c("2018-10-01", "2023-10-31")), expand = c(0,0), labels = NULL)+
  #rectangles for DJF and JJA
  geom_rect(aes(xmin =as.Date(c("2018-12-01", "2019-12-01", "2020-12-01","2021-12-01", "2022-12-01")), 
                xmax = as.Date(c("2019-02-28", "2020-02-28", "2021-02-28", "2022-02-28", "2023-02-28")), ymin = -Inf, ymax = Inf),fill = "lightblue",  alpha = 0.6) +
  geom_rect(aes(xmin =as.Date(c("2019-06-01","2020-06-01", "2021-06-01", "2022-06-01", "2023-06-01")), 
                xmax = as.Date(c("2019-08-30","2020-08-30", "2021-08-30", "2022-08-30", "2023-08-30")), ymin = -Inf, ymax = Inf),fill = "coral1", alpha = 0.6) +
  
  MStheme
ptimeSMJOR

ggplotly(ptimeSMJOR)
#time plots with soil and soil water
ptime1JOR <- ggplot()+
  #these are the data points re-measured at UMich
  #geom_line(data = JORprecip, aes(x = date.collection, y = d18O.mean), size = 0.5, color = JORcolor)+
  #geom_point(data = JORprecip, aes(x = date.collection, y = d18O.mean), size = 3, fill = JORcolor, shape = 23)+
  geom_line(data = PicarroData, aes(x = as_date(CollectionDate), y = d18O_updated), size = 0.5, color = JORcolor)+
  geom_point(data = PicarroData, aes(x = as_date(CollectionDate), y = d18O_updated, shape = Site), size = 3, fill = JORcolor)+
  scale_shape_manual(values = c(23,22))+  
  geom_line(data = JORswdeep, aes(x = date.collection, y = d18O.mean), size = 0.5, color = "gray90")+
  geom_line(data = JORswsurf, aes(x = date.collection, y = d18O.mean), size = 0.5, color = "black")+
  geom_point(data = JORsw, aes(x = date.collection, y = d18O.mean, color = bottom.depth), size = 3)+
  scale_color_gradient(high = "grey90", low = "black")+
  scale_x_date(date_breaks = "2 month", limits = as.Date(c("2018-10-01", "2023-10-31")), expand = c(0,0), labels = NULL)+
  xlab("")+
  ylab(expression(paste(delta^{18}, "O (‰, SMOW)")))+
  #rectangles for DJF and JJA
  geom_rect(aes(xmin =as.Date(c("2018-12-01", "2019-12-01", "2020-12-01","2021-12-01", "2022-12-01")), 
                xmax = as.Date(c("2019-02-28", "2020-02-28", "2021-02-28", "2022-02-28", "2023-02-28")), ymin = -Inf, ymax = Inf),fill = "lightblue",  alpha = 0.6) +
  geom_rect(aes(xmin =as.Date(c("2019-06-01","2020-06-01", "2021-06-01", "2022-06-01", "2023-06-01")), 
                xmax = as.Date(c("2019-08-30","2020-08-30", "2021-08-30", "2022-08-30", "2023-08-30")), ymin = -Inf, ymax = Inf),fill = "coral1", alpha = 0.6) +
  
  MStheme
ptime1JOR

ptime2JOR <-  ggplot()+
  #geom_point(data = JORprecip, aes(x = date.collection, y = d_excess.mean), size = 3, fill = JORcolor, shape = 23)+
  #geom_line(data = JORprecip, aes(x = date.collection, y = d_excess.mean), size = 0.5, color = JORcolor)+
  geom_line(data = PicarroData, aes(x = as_date(CollectionDate), y = d2H_updated), size = 0.5, color = JORcolor)+
  geom_point(data = PicarroData, aes(x = as_date(CollectionDate), y = d2H_updated, shape = Site), size = 3, fill = JORcolor)+
  scale_shape_manual(values = c(23,22))+  
  
  geom_line(data = JORswsurf, aes(x = date.collection, y = d_excess.mean), size = 0.5, color = "black")+
  geom_line(data = JORswdeep, aes(x = date.collection, y = d_excess.mean), size = 0.5, color = "gray90")+
  geom_point(data = JORsw, aes(x = date.collection, y = d_excess.mean, color = bottom.depth), size = 3)+
  scale_color_gradient(high = "grey90", low = "black")+
  scale_x_date(date_breaks = "2 month", limits = as.Date(c("2018-10-01", "2023-10-31")), expand = c(0,0), labels = NULL)+
  ylab(expression(paste(italic("d"),"-excess (‰, SMOW)"))) +
  xlab("")+
  #rectangles for DJF and JJA
  geom_rect(aes(xmin =as.Date(c("2018-12-01", "2019-12-01", "2020-12-01","2021-12-01", "2022-12-01")), 
                xmax = as.Date(c("2019-02-28", "2020-02-28", "2021-02-28", "2022-02-28", "2023-02-28")), ymin = -Inf, ymax = Inf),fill = "lightblue",  alpha = 0.6) +
  geom_rect(aes(xmin =as.Date(c("2019-06-01","2020-06-01", "2021-06-01", "2022-06-01", "2023-06-01")), 
                xmax = as.Date(c("2019-08-30","2020-08-30", "2021-08-30", "2022-08-30", "2023-08-30")), ymin = -Inf, ymax = Inf),fill = "coral1", alpha = 0.6) +
  
  MStheme
ptime2JOR

ptime3JOR <- ggplot()+
  geom_line(data = JORprecip, aes(x = date.collection, y = D17O_pmg), size = 0.5, color = JORcolor)+
  
  geom_errorbar(data  = JORprecip, aes(x = date.collection, y = D17O_pmg, ymin = D17O_pmg - D17O_err, ymax = D17O_pmg + D17O_err), width = 0, color = JORcolor)+
    geom_point(data = JORprecip, aes(x = date.collection, y = D17O_pmg, shape = siteID.2), size = 3, fill = JORcolor)+
  scale_shape_manual(values = c(23,22))+  
  
  geom_errorbar(data  = JORsw, aes(x = date.collection, y = D17O_pmg, ymin = D17O_pmg - D17O_err, ymax = D17O_pmg + D17O_err), width = 0, color = "black")+
  
  geom_point(data = JORsw, aes(x = date.collection, y = D17O_pmg, color = bottom.depth), size = 3)+
  scale_color_gradient(high = "grey90", low = "black")+
  ylab(expression(Delta^"'17"*"O (pmg, SMOW-SLAP)"))+
  xlab("Date")+
  scale_x_date(date_breaks = "2 month", date_labels = "%m/%y", limits = as.Date(c("2018-10-01", "2023-10-31")), expand = c(0,0), labels = NULL)+
  #rectangles for DJF and JJA
  geom_rect(aes(xmin =as.Date(c("2018-12-01", "2019-12-01", "2020-12-01","2021-12-01", "2022-12-01")), 
                xmax = as.Date(c("2019-02-28", "2020-02-28", "2021-02-28", "2022-02-28", "2023-02-28")), ymin = -Inf, ymax = Inf),fill = "lightblue",  alpha = 0.6) +
  geom_rect(aes(xmin =as.Date(c("2019-06-01","2020-06-01", "2021-06-01", "2022-06-01", "2023-06-01")), 
                xmax = as.Date(c("2019-08-30","2020-08-30", "2021-08-30", "2022-08-30", "2023-08-30")), ymin = -Inf, ymax = Inf),fill = "coral1", alpha = 0.6) +
  
  MStheme

ptime3JOR

ptimeJOR <- ggarrange(ptimePJOR, ptimeSMJOR, ptime1JOR, ptime2JOR, ptime3JOR, nrow = 5, common.legend = TRUE)
ptimeJOR

ggsave(filename = "ptime_JOR.pdf", plot = ptimeJOR, device = cairo_pdf, height=8.25,width=10.75, path = path.to.figs.refined)


######### MOJAVE TIME PLOTS#########

MOJprecip <- MOJmetw %>% filter(Water.Type == "Precipitation" | Water.Type == "Snow")
MOJswdeep <- MOJsw %>% filter(bottom.depth>55)
MOJswsurf <- MOJsw %>% filter(bottom.depth<20)

MOJmetw.unevap <- MOJmetw %>% filter(Water.Type == "Precipitation" | Water.Type == "Spring" | Water.Type == "Snow")
#ignores well samples, puddle samples

#read in climo data, using CRS for simplicity
CRSclim <- read.table(file = "~/Documents/Soil Water 17O Manuscript/data/climo data/Mojave climo data/Creosote_daily.txt", skip = 7, header = T, sep = ",")
CRSclim$Date <- as_date(CRSclim$Date)

ptimePMOJ <- ggplot()+
  geom_line(data = CRSclim, aes(x = Date, y = Rain_mm_Tot_sum), color = MOJcolor, size = 0.5)+
  xlab("")+
  ylab("Daily Precipitation (mm)")+
  scale_y_continuous(limits = c(0, 40), expand = c(0,0))+
  #scale_x_date(date_breaks = "2 month", limits = as.Date(c("2021-10-01", "2023-09-30")), expand = c(0,0), labels = NULL)+
  #rectangles for DJF and JJA
  geom_rect(aes(xmin =as.Date(c("2021-12-01", "2022-12-01")), 
                xmax = as.Date(c("2022-02-28", "2023-02-28")), ymin = -Inf, ymax = Inf),fill = "lightblue",  alpha = 0.6) +
  geom_rect(aes(xmin =as.Date(c("2022-06-01", "2023-06-01")), 
                xmax = as.Date(c( "2022-08-30", "2023-08-30")), ymin = -Inf, ymax = Inf),fill = "coral1", alpha = 0.6) +
  MStheme

ggplotly(ptimePMOJ)

ptimeSMMOJ <- ggplot()+
  geom_line(data = CRSclim, aes(x = Date, y= VWC_1_5_avg), color = "black")+
  geom_line(data = CRSclim, aes(x = Date, y= VWC_1_25_avg), color = "gray50")+
  geom_line(data = CRSclim, aes(x = Date, y= VWC_1_50_avg), color = "gray70")+
  geom_line(data = CRSclim, aes(x = Date, y= VWC_1_125_avg), color = "gray90")+
  xlab("")+
  ylab("Soil Moisture (-)")+
  scale_y_continuous(limits = c(0, 0.18), expand = c(0,0))+
  scale_x_date(date_breaks = "2 month", limits = as.Date(c("2021-10-01", "2023-09-30")), expand = c(0,0), labels = NULL)+
  #rectangles for DJF and JJA
  geom_rect(aes(xmin =as.Date(c("2021-12-01", "2022-12-01")), 
                xmax = as.Date(c("2022-02-28", "2023-02-28")), ymin = -Inf, ymax = Inf),fill = "lightblue",  alpha = 0.6) +
  geom_rect(aes(xmin =as.Date(c("2022-06-01", "2023-06-01")), 
                xmax = as.Date(c( "2022-08-30", "2023-08-30")), ymin = -Inf, ymax = Inf),fill = "coral1", alpha = 0.6) +
  MStheme
ptimeSMMOJ

ptimeCO2MOJ <- ggplot()+
  geom_line(data = CRSclim, aes(x = Date, y= CO2_1_5_avg), color = "black")+
  geom_line(data = CRSclim, aes(x = Date, y= CO2_1_25_avg), color = "gray50")+
  geom_line(data = CRSclim, aes(x = Date, y= CO2_1_50_avg), color = "gray70")+
  geom_line(data = CRSclim, aes(x = Date, y= CO2_1_125_avg), color = "gray90")+
  xlab("")+
  ylab("pCO2 (ppm)")+
  #rectangles for DJF and JJA
  geom_rect(aes(xmin =as.Date(c("2021-12-01", "2022-12-01")), 
                xmax = as.Date(c("2022-02-28", "2023-02-28")), ymin = -Inf, ymax = Inf),fill = "lightblue",  alpha = 0.6) +
  geom_rect(aes(xmin =as.Date(c("2022-06-01", "2023-06-01")), 
                xmax = as.Date(c( "2022-08-30", "2023-08-30")), ymin = -Inf, ymax = Inf),fill = "coral1", alpha = 0.6) +
  #scale_y_continuous(limits = c(0, 0.18), expand = c(0,0))+
#  scale_x_date(date_breaks = "2 month", limits = as.Date(c("2021-10-01", "2023-09-30")), expand = c(0,0), labels = NULL)+
  MStheme
ptimeCO2MOJ

#time plots with soil and soil water
ptime1MOJ <- ggplot()+
  #geom_point(data = MOJsw, aes(x = date.of.collection, y = d18O, shape = siteID.2), fill = MOJcolor, size = 3)+
  geom_line(data = MOJswsurf %>% filter(siteID.2 == "CRS"), aes(x = date.collection, y = d18O.mean), color = "black", size = 0.5)+
  geom_line(data = MOJswdeep %>% filter(siteID.2 == "CRS"), aes(x = date.collection, y = d18O.mean), color = "gray90", size = 0.5)+
  geom_line(data = MOJprecip, aes(x = date.collection, y = d18O.mean), color = MOJcolor, size = 0.5)+
  geom_point(data = MOJsw, aes(x = date.collection, y = d18O.mean, fill = bottom.depth, shape = siteID.2), size = 3)+
  scale_shape_manual(values = c(21,23, 24, 25))+
  scale_fill_gradient(high = "grey90", low = "black")+
  geom_point(data = MOJprecip, aes(x = date.collection, y = d18O.mean), fill = MOJcolor, size = 3, shape = 22)+
  scale_x_date(date_breaks = "2 month", limits = as.Date(c("2021-10-01", "2023-09-30")), expand = c(0,0), labels = NULL)+
  #rectangles for DJF and JJA
  geom_rect(aes(xmin =as.Date(c("2021-12-01", "2022-12-01")), 
                xmax = as.Date(c("2022-02-28", "2023-02-28")), ymin = -Inf, ymax = Inf),fill = "lightblue",  alpha = 0.6) +
  geom_rect(aes(xmin =as.Date(c("2022-06-01", "2023-06-01")), 
                xmax = as.Date(c( "2022-08-30", "2023-08-30")), ymin = -Inf, ymax = Inf),fill = "coral1", alpha = 0.6) +
  ylab(expression(paste(delta^{18}, "O (‰, SMOW)")))+
  xlab("")+
  MStheme
ptime1MOJ

ptime2MOJ <- ggplot()+
  #geom_point(data = MOJsw, aes(x = date.of.collection, y = d18O, shape = siteID.2), fill = MOJcolor, size = 3)+
  geom_line(data = MOJswsurf %>% filter(siteID.2 == "CRS"), aes(x = date.collection, y = d_excess.mean), color = "black", size = 0.5)+
  geom_line(data = MOJswdeep %>% filter(siteID.2 == "CRS"), aes(x = date.collection, y = d_excess.mean), color = "gray90", size = 0.5)+
  geom_line(data = MOJprecip, aes(x = date.collection, y = d_excess.mean), color = MOJcolor, size = 0.5)+
  geom_point(data = MOJsw, aes(x = date.collection, y = d_excess.mean, fill = bottom.depth, shape = siteID.2), size = 3)+
  scale_shape_manual(values = c(21,23, 24, 25))+
  scale_fill_gradient(high = "grey90", low = "black")+
  geom_point(data = MOJprecip, aes(x = date.collection, y = d_excess.mean), fill = MOJcolor, size = 3, shape = 22)+
  scale_x_date(date_breaks = "2 month", limits = as.Date(c("2021-10-01", "2023-09-30")), expand = c(0,0), labels = NULL)+
  ylab(expression(paste(italic("d"),"-excess (‰, SMOW)"))) +
  xlab("")+
  #rectangles for DJF and JJA
  geom_rect(aes(xmin =as.Date(c("2021-12-01", "2022-12-01")), 
                xmax = as.Date(c("2022-02-28", "2023-02-28")), ymin = -Inf, ymax = Inf),fill = "lightblue",  alpha = 0.6) +
  geom_rect(aes(xmin =as.Date(c("2022-06-01", "2023-06-01")), 
                xmax = as.Date(c( "2022-08-30", "2023-08-30")), ymin = -Inf, ymax = Inf),fill = "coral1", alpha = 0.6) +
  MStheme
ptime2MOJ


ptime3MOJ <- ggplot()+
  geom_line(data = MOJswsurf %>% filter(siteID.2 == "CRS"), aes(x = date.collection, y = D17O_pmg), color = "black", size = 0.5)+
  geom_line(data = MOJswdeep %>% filter(siteID.2 == "CRS"), aes(x = date.collection, y = D17O_pmg), color = "gray90", size = 0.5)+
  geom_line(data = MOJprecip, aes(x = date.collection, y = D17O_pmg), color = MOJcolor, size = 0.5)+
  
  geom_errorbar(data  = MOJsw, aes(x = date.collection, y = D17O_pmg, ymin = D17O_pmg - D17O_err, ymax = D17O_pmg + D17O_err), width = 0, color = "black")+
  
  geom_point(data = MOJsw, aes(x = date.collection, y = D17O_pmg, fill = bottom.depth, shape = siteID.2), size = 3)+
  scale_shape_manual(values = c(21,23, 24, 25))+
  scale_fill_gradient(high = "grey90", low = "black")+
  geom_errorbar(data  = MOJprecip, aes(x = date.collection, y = D17O_pmg, ymin = D17O_pmg - D17O_err, ymax = D17O_pmg + D17O_err), width = 0, color = MOJcolor)+
  geom_point(data = MOJprecip, aes(x = date.collection, y = D17O_pmg), fill = MOJcolor, size = 3, shape = 22)+
  scale_x_date(date_breaks = "2 month", date_labels = "%m/%y", limits = as.Date(c("2021-10-01", "2023-09-30")), expand = c(0,0))+
  ylab(expression(Delta^"'17"*"O (pmg, SMOW-SLAP)"))+
  xlab("")+
  #rectangles for DJF and JJA
  geom_rect(aes(xmin =as.Date(c("2021-12-01", "2022-12-01")), 
                xmax = as.Date(c("2022-02-28", "2023-02-28")), ymin = -Inf, ymax = Inf),fill = "lightblue",  alpha = 0.6) +
  geom_rect(aes(xmin =as.Date(c("2022-06-01", "2023-06-01")), 
                xmax = as.Date(c( "2022-08-30", "2023-08-30")), ymin = -Inf, ymax = Inf),fill = "coral1", alpha = 0.6) +
  MStheme
ptime3MOJ

#attempt to skootch everything closer
timep<- ggarrange(ptimePMOJ,NULL, ptimeSMMOJ, NULL,ptime1MOJ,NULL, ptime2MOJ,NULL, ptime3MOJ, 
                  nrow = 9, ncol = 1, heights = c(1,-0.2, 1, -0.2, 1, -0.2, 1, -0.2, 1), align = "hv", common.legend = TRUE)

#timep<- ggarrange(ptimePMOJ, ptimeSMMOJ, ptime1MOJ, ptime2MOJ, ptime3MOJ, nrow = 5, ncol = 1, align = "hv", common.legend = TRUE)

timep
ggsave(filename = "ptime_MOJ.pdf", plot = timep, device = cairo_pdf, width=10.75,height=8.25, path = path.to.figs.refined)


ptimeSMMOJ <- ggplot()+
  #geom_line(data = CRSclim, aes(x = Date, y= VWC_1_5_avg), color = "black")+
  geom_line(data = CRSclim, aes(x = Date, y= VWC_1_25_avg), color = MOJcolor)+
  geom_line(data = CRSclim, aes(x = Date, y= VWC_1_50_avg), color = MOJcolor.minor2)+
  geom_line(data = CRSclim, aes(x = Date, y= VWC_1_125_avg), color = MOJcolor.minor)+
  xlab("")+
  ylab("Soil Moisture (-)")+
  scale_y_continuous(limits = c(0, 0.18), expand = c(0,0))+
  #scale_x_date(date_breaks = "2 month", limits = as.Date(c("2021-10-01", "2023-09-30")), expand = c(0,0), labels = NULL)+
  scale_x_date(date_breaks = "2 month", date_labels = "%m/%y", limits = as.Date(c("2021-07-01", "2023-09-30")), expand = c(0,0))+

  MStheme+
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 16)) #only formats the y axis so the x axis label can appear only on bottom panel?

ptimeSMMOJ

ggsave(filename = "ptimeSMMOJ.pdf", plot = ptimeSMMOJ, device = cairo_pdf, height=6.44,width=10.31, path = path.to.figs)

ptimeSTMOJ <- ggplot()+
  #geom_line(data = CRSclim, aes(x = Date, y= TS_1_5_avg), color = "black")+
  geom_line(data = CRSclim, aes(x = Date, y= TS_1_25_avg), color = MOJcolor)+
  geom_line(data = CRSclim, aes(x = Date, y= TS_1_50_avg), color = MOJcolor.minor2)+
  geom_line(data = CRSclim, aes(x = Date, y= TS_1_125_avg), color = MOJcolor.minor)+
  xlab("")+
  ylab("Soil Temperature (°C)")+
  #scale_y_continuous(limits = c(0, 0.18), expand = c(0,0))+
  #scale_x_date(date_breaks = "2 month", limits = as.Date(c("2021-10-01", "2023-09-30")), expand = c(0,0), labels = NULL)+
  scale_x_date(date_breaks = "2 month", date_labels = "%m/%y", limits = as.Date(c("2021-07-01", "2023-09-30")), expand = c(0,0))+
  
  MStheme+
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 16)) #only formats the y axis so the x axis label can appear only on bottom panel?

ptimeSTMOJ
ggsave(filename = "ptimeSTMOJ.pdf", plot = ptimeSTMOJ, device = cairo_pdf, height=6.44,width=10.31, path = path.to.figs)



ptimeCO2MOJ <- ggplot()+
  #geom_line(data = CRSclim, aes(x = Date, y= CO2_1_5_avg), color = "black")+
  geom_line(data = CRSclim, aes(x = Date, y= CO2_1_25_avg), color = MOJcolor)+
  geom_line(data = CRSclim, aes(x = Date, y= CO2_1_50_avg), color = MOJcolor.minor2)+
  geom_line(data = CRSclim, aes(x = Date, y= CO2_1_125_avg), color = MOJcolor.minor)+
  xlab("")+
  ylab("pCO2 (ppm)")+
  #scale_y_continuous(limits = c(0, 0.18), expand = c(0,0))+
  scale_x_date(date_breaks = "2 month", date_labels = "%m/%y", limits = as.Date(c("2021-07-01", "2023-09-30")), expand = c(0,0))+
  MStheme+
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 16)) #only formats the y axis so the x axis label can appear only on bottom panel?

ptimeCO2MOJ

ggsave(filename = "ptimeCO2MOJ.pdf", plot = ptimeCO2MOJ, device = cairo_pdf, height=6.44,width=10.31, path = path.to.figs)




######### MOJAVE three sites comparison plots######

#read in climo data, using CRS for simplicity
CRSclim <- read.table(file = "~/Documents/Soil Water 17O Manuscript/data/climo data/Mojave climo data/Creosote_daily.txt", skip = 7, header = T, sep = ",")
CRSclim$Date <- as_date(CRSclim$Date)

JTclim <- read.table(file = "~/Documents/Soil Water 17O Manuscript/data/climo data/Mojave climo data/JTree_daily.txt", skip = 7, header = T, sep = ",")
JTclim$Date <- as_date(JTclim$Date)

PJclim <- read.table(file = "~/Documents/Soil Water 17O Manuscript/data/climo data/Mojave climo data/PJ_daily.txt", skip = 7, header = T, sep = ",")
PJclim$Date <- as_date(PJclim$Date)

#dates of rainfall collection
DSCprecip <- MOJprecip %>% filter(siteID.2 == "DSC")
DSCprecip$Date <- as_date(DSCprecip$date.collection)

ptimePMOJ_compare <- ggplot()+
  #geom_line(data = PJclim, aes(x = Date, y = Rain_mm_Tot_sum), color = "red", size = 0.5)+ #removed because this data is so spotty its not useful to plot
  
  geom_line(data = CRSclim, aes(x = Date, y = Rain_mm_Tot_sum), color = MOJcolor, size = 0.5)+
  geom_line(data = JTclim, aes(x = Date, y = Rain_mm_Tot_sum), color = "orange", size = 0.5)+
  geom_point(data = DSCprecip, aes(x = Date, y= 10), size = 2, alpha = 0.5)+  
  xlab("")+
  ylab("Daily Precipitation (mm)")+
  scale_y_continuous(limits = c(0, 40), expand = c(0,0))+
  scale_x_date(date_breaks = "2 month", limits = as.Date(c("2021-10-01", "2023-09-30")), expand = c(0,0))+
  MStheme

ptimePMOJ_compare

ptimeSMMOJ_compare <- ggplot()+
  geom_line(data = CRSclim, aes(x = Date, y= VWC_1_25_avg), color = MOJcolor)+
  geom_line(data = JTclim, aes(x = Date, y= VWC_1_25_avg), color = "orange")+
  geom_line(data = PJclim, aes(x = Date, y= VWC_1_25_avg), color = "red")+
  geom_point(data = DSCprecip, aes(x = Date, y= 0.10), size = 2, alpha = 0.5)+  
  
  xlab("")+
  ylab("Soil Moisture (-)")+
  scale_y_continuous(limits = c(0, 0.18), expand = c(0,0))+
  scale_x_date(date_breaks = "2 month", limits = as.Date(c("2021-10-01", "2023-09-30")), expand = c(0,0))+
  MStheme
ptimeSMMOJ_compare

ggarrange(ptimePMOJ_compare, ptimeSMMOJ_compare, nrow = 2)

ptimeCO2MOJ_compare <- ggplot()+
  
  geom_line(data = CRSclim, aes(x = Date, y= CO2_1_25_avg), color = MOJcolor)+
  geom_line(data = JTclim, aes(x = Date, y= CO2_1_25_avg), color = "orange")+
  geom_line(data = PJclim, aes(x = Date, y= CO2_1_25_avg), color = "red")+
  
  xlab("")+
  ylab("pCO2 (ppm)")+
  #scale_y_continuous(limits = c(0, 0.18), expand = c(0,0))+
  #  scale_x_date(date_breaks = "2 month", limits = as.Date(c("2021-10-01", "2023-09-30")), expand = c(0,0), labels = NULL)+
  MStheme
ptimeCO2MOJ_compare



######## MICHIGAN TIME PLOTS #########

#precipitation data
NOAAclim <- read.csv(file = "~/Documents/Soil Water 17O Manuscript/data/climo data/ESGR climo data/3886033.csv")
ESGRdailyP <- NOAAclim %>%
  group_by(DATE) %>%
  summarise(PRCP = mean(PRCP, na.rm = TRUE))

ESGRdailyP$YEAR <- year(ESGRdailyP$DATE)

ESGR_AnnualP <- ESGRdailyP %>%
  group_by(YEAR) %>%
  summarise(YEAR_PRCP = sum(PRCP, na.rm = TRUE))

#soil moisutre data, cleaned up in code "ESGR_logger_viz.R", data presented in Kelson 2024 papers
ESGRSM <- read.csv(file = "~/Documents/Soil Water 17O Manuscript/data/climo data/ESGR climo data/ESGR_SM_combined_forR_clean.csv") #some continuous data
#needs a little more clean up
ESGRSM <- add_row(ESGRSM, X = 26024.4, RecordNo = 1, Date_Time_GMT.05 = "10/27/21 16:28", 
                  WaterContent_Deep = NA, WaterContent_Mid = NA, WaterContent_Shallow = NA, Date = "2021-10-27 02:39:00", .before = 23436)

ESGRSM <- add_row(ESGRSM, X = 27147.4, RecordNo = 1, Date_Time_GMT.05 = "12/11/21 12:01", 
                  WaterContent_Deep = NA, WaterContent_Mid = NA, WaterContent_Shallow = NA, Date = "2021-12-11 12:01:00", .before = 24428)

ESGRSM <- add_row(ESGRSM, X = 27147.4, RecordNo = 1, Date_Time_GMT.05 = "07/08/22 16:22", 
                  WaterContent_Deep = NA, WaterContent_Mid = NA, WaterContent_Shallow = NA, Date = "2022-07-08 16:22:00", .before = 25053)

ESGRSM <- add_row(ESGRSM, X = 25053.5, RecordNo = 137, Date_Time_GMT.05 = "08/12/21 0:09", 
                  WaterContent_Deep = NA, WaterContent_Mid = NA, WaterContent_Shallow = NA, Date = "2021-08-12 00:09:00", .before = 22465)


ESGRSM <- ESGRSM %>% filter(as_date(Date) < as_date("2022-07-16 23:37:00")) #remove bad data at the end


ptimePESGR <- ggplot()+
  geom_line(data = ESGRdailyP, aes(x = as_date(DATE), y = PRCP), color = ESGRcolor)+
  xlab("")+
  ylab("Daily Precipitation (mm)")+
  scale_y_continuous(limits = c(0, 41), expand = c(0,0))+
  scale_x_date(date_breaks = "2 month", limits = as_date(c("2020-10-01", "2022-11-30")), expand = c(0,0), labels = NULL)+
  geom_rect(aes(xmin =as.Date(c("2020-12-01", "2021-12-01")), xmax = as.Date(c("2021-02-28", "2022-02-28")), ymin = -Inf, ymax = Inf),fill = "lightblue",  alpha = 0.6) +
  geom_rect(aes(xmin =as.Date(c("2021-06-01", "2022-06-01")), xmax = as.Date(c("2021-08-30", "2022-08-30")), ymin = -Inf, ymax = Inf),fill = "coral1", alpha = 0.6) +
  MStheme
ptimePESGR


ptimeSMESGR <- ggplot()+
  geom_line(data = ESGRSM, aes(x = as_date(Date), y = WaterContent_Shallow), color = "black") + #10 cm
  geom_line(data = ESGRSM, aes(x = as_date(Date), y = WaterContent_Mid), color = "gray60") + # 30 cm
  geom_line(data = ESGRSM, aes(x = as_date(Date), y = WaterContent_Deep), color = "gray90") + #60 cm
  xlab("")+
  ylab("Soil Moisture (-)")+
  scale_x_date(date_breaks = "2 month", limits = as_date(c("2020-10-01", "2022-11-30")), expand = c(0,0), labels = NULL)+
  geom_rect(aes(xmin =as.Date(c("2020-12-01", "2021-12-01")), xmax = as.Date(c("2021-02-28", "2022-02-28")), ymin = -Inf, ymax = Inf),fill = "lightblue",  alpha = 0.6) +
  geom_rect(aes(xmin =as.Date(c("2021-06-01", "2022-06-01")), xmax = as.Date(c("2021-08-30", "2022-08-30")), ymin = -Inf, ymax = Inf),fill = "coral1", alpha = 0.6) +
  MStheme
ptimeSMESGR


ESGRprecip <- ESGRmetw %>% filter(Water.Type == "Precipitation")
ESGRswsurf <- ESGRsw %>% filter(bottom.depth < 15)
ESGRswdeep <- ESGRsw %>% filter(bottom.depth > 49)

ptime1ESGR <- ggplot()+
  geom_line(data = ESGRprecip, aes(x = date.collection, y = d18O.mean), size = 0.5, color = ESGRcolor)+
  geom_point(data = ESGRprecip, aes(x = date.collection, y = d18O.mean), size = 3, fill = ESGRcolor, shape  = 22)+
  geom_line(data = ESGRswdeep, aes(x = date.collection, y = d18O.mean), size = 0.5, color = "gray90")+
  geom_line(data = ESGRswsurf, aes(x = date.collection, y = d18O.mean), size = 0.5, color = "black")+
  geom_point(data = ESGRsw, aes(x = date.collection, y = d18O.mean, color = bottom.depth), size = 3)+
  scale_color_gradient(high = "grey90", low = "black")+
  # geom_line(data = REYswsurf, aes(x = date.collection, y = d18O.mean), size = 1, color = REYcolor.minor)+#doesnt look good
  xlab("")+
  ylab(expression(paste(delta^{18}, "O (‰, SMOW)")))+
  scale_x_date(date_breaks = "2 month", limits = as_date(c("2020-10-01", "2022-11-30")), expand = c(0,0), labels = NULL)+
  MStheme
ptime1ESGR
#the lines aren't connecting because sample 2000 doesn't have [good] Picarro data

ptime2ESGR <-  ggplot()+
  geom_line(data = ESGRprecip, aes(x = date.collection, y = d_excess.mean), size = 0.5, color = ESGRcolor)+
  geom_point(data = ESGRprecip, aes(x = date.collection, y = d_excess.mean), size = 3, fill = ESGRcolor, shape  = 22)+
  geom_line(data = ESGRswdeep, aes(x = date.collection, y = d_excess.mean), size = 0.5, color = "gray90")+
  geom_line(data = ESGRswsurf, aes(x = date.collection, y = d_excess.mean), size = 0.5, color = "black")+
  geom_point(data = ESGRsw, aes(x = date.collection, y = d_excess.mean, color = bottom.depth), size = 3)+
  scale_color_gradient(high = "grey90", low = "black")+
  ylab(expression(paste(italic("d"),"-excess (‰, SMOW)"))) +
  xlab("")+
  scale_x_date(date_breaks = "2 month", limits = as_date(c("2020-10-01", "2022-11-30")), expand = c(0,0), labels = NULL)+
  geom_rect(aes(xmin =as.Date(c("2020-12-01", "2021-12-01")), xmax = as.Date(c("2021-02-28", "2022-02-28")), ymin = -Inf, ymax = Inf),fill = "lightblue",  alpha = 0.6) +
  geom_rect(aes(xmin =as.Date(c("2021-06-01", "2022-06-01")), xmax = as.Date(c("2021-08-30", "2022-08-30")), ymin = -Inf, ymax = Inf),fill = "coral1", alpha = 0.6) +
  MStheme 

#the lines aren't connecting because sample 2000 doesn't have [good] Picarro data
ptime2ESGR

ptime3ESGR <- ggplot()+
  geom_line(data = ESGRprecip, aes(x = date.collection, y = D17O_pmg), size = 0.5, color = ESGRcolor)+
  geom_errorbar(data  = ESGRprecip, aes(x = date.collection, y = D17O_pmg, ymin = D17O_pmg - D17O_err, ymax = D17O_pmg + D17O_err), width = 0, color = ESGRcolor)+
  geom_point(data = ESGRprecip, aes(x = date.collection, y = D17O_pmg), size = 3, fill = ESGRcolor, shape  = 22)+
  geom_line(data = ESGRswdeep, aes(x = date.collection, y = D17O_pmg), size = 0.5, color = "gray90")+
  geom_line(data = ESGRswsurf, aes(x = date.collection, y = D17O_pmg), size = 0.5, color = "black")+
  geom_errorbar(data  = ESGRsw, aes(x = date.collection, y = D17O_pmg, ymin = D17O_pmg - D17O_err, ymax = D17O_pmg + D17O_err, color = bottom.depth), width = 0)+
  geom_point(data = ESGRsw, aes(x = date.collection, y = D17O_pmg, color = bottom.depth), size = 3)+
  scale_color_gradient(high = "grey90", low = "black")+
  ylab(expression(Delta^"'17"*"O (pmg, SMOW-SLAP)"))+
  xlab("")+
  scale_x_date(date_breaks = "2 month", date_labels = "%m/%y", limits = as_date(c("2020-10-01", "2022-11-30")), expand = c(0,0))+
  #rectangles for DJF and JJA
  geom_rect(aes(xmin =as.Date(c("2020-12-01", "2021-12-01")), xmax = as.Date(c("2021-02-28", "2022-02-28")), ymin = -Inf, ymax = Inf),fill = "lightblue",  alpha = 0.6) +
  geom_rect(aes(xmin =as.Date(c("2021-06-01", "2022-06-01")), xmax = as.Date(c("2021-08-30", "2022-08-30")), ymin = -Inf, ymax = Inf),fill = "coral1", alpha = 0.6) +
  MStheme 

ptime3ESGR

ptimeESGR <- ggarrange(ptimePESGR,ptimeSMESGR, ptime1ESGR, ptime2ESGR, ptime3ESGR, nrow = 5, ncol = 1, common.legend = TRUE)
ptimeESGR
ggsave(filename = "ptime_ESGR.pdf", plot = ptimeESGR, device = cairo_pdf, width=10.75,height=8.25, path = path.to.figs.refined)
