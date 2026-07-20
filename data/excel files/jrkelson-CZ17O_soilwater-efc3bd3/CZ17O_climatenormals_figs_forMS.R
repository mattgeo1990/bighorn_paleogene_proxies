#Make climate normals figures for the 4 study sites considered
#started JRK January 2025

#to include: monthly precipitation, relative humidity of overlying air, air temp (min, mean, max), 
#for each site make a dataframe with climate normals.

library(ggplot2)
library(dplyr)
library(lubridate)
library(ggpubr)
library(plotly)

#####set up plotting theme and site colors########
MStheme <- theme_bw()+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(axis.title = element_text(size = 12),axis.text = element_text(size = 12))

path.to.figs <- "~/Documents/Soil Water 17O Manuscript/code_figures/MS_draft_figs/"
path.to.figs.climate.normals <- "~/Documents/Soil Water 17O Manuscript/code_figures/MS_draft_figs/climate_normals_drafts/"

path.to.figs.refined <- "~/Documents/Soil Water 17O Manuscript/code_figures/MS_refined_figs/"

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
####Mojave########
CRSclim <- read.table(file = "~/Documents/Soil Water 17O Manuscript/data/climo data/Mojave climo data/Creosote_daily.txt", skip = 7, header = T, sep = ",")
CRSclim$Date <- as_date(CRSclim$Date)
CRSclim$Month <- month(CRSclim$Date)
CRSclim$Year <- year(CRSclim$Date)

CRSclim_monthly <- read.table(file = "~/Documents/Soil Water 17O Manuscript/data/climo data/Mojave climo data/Creosote_monthly.txt", skip = 7, header = T, sep = ",")
CRSclim_monthly$Month <- month(CRSclim_monthly$Date)
CRSclim_monthly$Year <- year(CRSclim_monthly$Date)

JTclim_monthly <- read.table(file = "~/Documents/Soil Water 17O Manuscript/data/climo data/Mojave climo data/JTree_monthly.txt", skip = 7, header = T, sep = ",")
JTclim_monthly$Month <- month(JTclim_monthly$Date)
JTclim_monthly$Year <- year(JTclim_monthly$Date)

PJclim_monthly <- read.table(file = "~/Documents/Soil Water 17O Manuscript/data/climo data/Mojave climo data/PJ_monthly.txt", skip = 7, header = T, sep = ",")
PJclim_monthly$Month <- month(PJclim_monthly$Date)
PJclim_monthly$Year <- year(PJclim_monthly$Date)

CRSMonthlyT <- CRSclim_monthly %>%
  group_by(Month) %>%
  summarise(meanAirT = mean(AirTC_Avg_avg, na.rm = TRUE)) 


CRSMinMonthlyT <- CRSclim_monthly %>%
  group_by(Month) %>%
  summarise(minAirT = mean(AirTC_Min_min, na.rm = TRUE)) 

CRSMaxMonthlyT <- CRSclim_monthly %>%
  group_by(Month) %>%
  summarise(maxAirT = max(AirTC_Max_max, na.rm = TRUE)) 

mean(CRSclim$AirTC_Avg_avg, na.rm = TRUE) #mAAT over many years

TempPlotCRS <- ggplot()+
  geom_point(data = CRSMonthlyT, aes(x = Month,y = meanAirT), fill = MOJcolor, color = MOJcolor, shape = 21 )+
  geom_line(data = CRSMonthlyT, aes(x = Month,y = meanAirT), color = MOJcolor)+
  
  scale_x_discrete(limits=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"))+
  scale_y_continuous(limits = c(-6,32), expand = c(0, 0), breaks=seq(-5,30,5)) +
  
  ylab("Air Temp (°C)")+
  xlab("Month")+
  MStheme

TempPlotCRS

CRSMonthlyRH <- CRSclim_monthly %>%
  group_by(Month) %>%
  summarise(meanRH = mean(RH_Avg_avg, na.rm = TRUE))
mean(CRSclim$RH_Avg_avg, na.rm = TRUE) #mean over many years

JTMonthlyRH <- JTclim_monthly %>%
  group_by(Month) %>%
  summarise(meanRH = mean(RH_Avg_avg, na.rm = TRUE))
mean(JTclim_monthly$RH_Avg_avg, na.rm = TRUE) #mean over many years

PJMonthlyRH <- PJclim_monthly %>%
  group_by(Month) %>%
  summarise(meanRH = mean(RH_Avg_avg, na.rm = TRUE))
mean(PJclim_monthly$RH_Avg_avg, na.rm = TRUE) #mean over many years


RHPlotCRS <- ggplot()+
  geom_point(data = CRSMonthlyRH, aes(x = Month, y = meanRH), fill = MOJcolor, color = MOJcolor, shape = 21)+
  geom_line(data = CRSMonthlyRH, aes(x = Month, y = meanRH), color = MOJcolor)+
  
  geom_line(data = JTMonthlyRH, aes(x = Month, y = meanRH), color = "red")+
  geom_line(data = PJMonthlyRH, aes(x = Month, y = meanRH), color = "yellow")+
  
  scale_x_discrete(limits=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"), labels = NULL)+
  scale_y_continuous(limits = c(10,80), expand = c(0, 0), breaks=seq(0,80,10)) +
  xlab("")+
  ylab("RH (%)")+
  MStheme
RHPlotCRS


CRSMonthlyP <- CRSclim_monthly%>%
  group_by(Month) %>%
  filter(!is.na(AirTC_Avg_avg)) %>%
  summarise(monthlyP = mean(Rain_mm_Tot_sum))  #can use mean because the monthly precip is already summed

sum(CRSMonthlyP$monthlyP) #this should be close to the MAP in mm/year

CRSAnnualP <- CRSclim%>%
  group_by(Year) %>%
  filter(!is.na(AirTC_Avg_avg)) %>%
  summarise(MAP = sum(Rain_mm_Tot_sum))  #can use mean because the monthly precip is already summed

mean(CRSAnnualP$MAP[CRSAnnualP$Year == 2018 | CRSAnnualP$Year == 2019 | CRSAnnualP$Year == 2022 | CRSAnnualP$Year == 2023]) #mean of mostly complete years, another estimate of MAP
sd(CRSAnnualP$MAP[CRSAnnualP$Year == 2018 | CRSAnnualP$Year == 2019 | CRSAnnualP$Year == 2022 | CRSAnnualP$Year == 2023]) #mean of mostly complete years, another estimate of MAP


MAPPlotCRS <- ggplot()+
  geom_bar(data = CRSMonthlyP, aes(x = Month,y = monthlyP), stat = "identity", color = MOJcolor, fill = MOJcolor)+
  scale_x_discrete(limits=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"), labels = NULL)+
  scale_y_continuous(limits = c(0,50), expand = c(0, 0), breaks=seq(0,50,10)) +
  xlab("")+
  ylab("Monthly Precip (mm)")+
  MStheme

MAPPlotCRS


ClimateNormalsCRS <- ggarrange(MAPPlotCRS, RHPlotCRS, TempPlotCRS,nrow = 3, ncol = 1)
ggsave(filename = "ClimateNormalsCRS.pdf", plot = ClimateNormalsCRS, device = cairo_pdf, height=4,width=4, path = path.to.figs.climate.normals)


####Reynolds######
#data for Reynolds creek is downloaded from USDA in their publicly available folder
#https://ars-usda.app.box.com/s/4jwgmxyxb8vacosvp1t5sdlibdc5qxoe

#files used in the code are uploaded in github for convenience

REYclim <- read.csv(file = "~/Documents/Soil Water 17O Manuscript/data/climo data/Reynolds climo data/reynolds-creek-098c-climate-l1.csv")
#climate data from Nancy GUlch met site, 098. metadata is in the header of the .dat file with the same name
REYclim$Month <- month(REYclim$datetime)
REYclim$Year <- year(REYclim$datetime)

REYclim$hum3[which(REYclim$hum3 == -999)] <- NA
REYclim$tmp3[which(REYclim$tmp3 == -999)] <- NA

REYMonthlyRH <- REYclim %>%
  mutate(Month = month(datetime)) %>%
  group_by(Month) %>%
  summarise(meanRH = mean(hum3, na.rm = TRUE))
mean(REYclim$hum3, na.rm = TRUE) #mean over many years

REYMonthlyT <- REYclim %>%
  mutate(Month = month(datetime)) %>%
  group_by(Month) %>%
  summarise(meanAirT = mean(tmp3, na.rm = TRUE)) 

REYMinMonthlyT <- REYclim %>%
  mutate(Month = month(datetime)) %>%
  group_by(Month) %>%
  summarise(minAirT = min(tmp3, na.rm = TRUE)) 

REYMaxMonthlyT <- REYclim %>%
  mutate(Month = month(datetime)) %>%
  group_by(Month) %>%
  summarise(maxAirT = max(tmp3, na.rm = TRUE)) 

mean(REYclim$tmp3, na.rm = TRUE) #mean over many years


RHPlotREY <- ggplot()+
  geom_point(data = REYMonthlyRH, aes(x = Month, y = meanRH), fill = REYcolor, color = REYcolor, shape = 21)+
  geom_line(data = REYMonthlyRH, aes(x = Month, y = meanRH), color = REYcolor)+
  scale_x_discrete(limits=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"), labels = NULL)+
  scale_y_continuous(limits = c(10,80), expand = c(0, 0), breaks=seq(0,80,10)) +
  xlab("")+
  ylab("RH (%)")+
  MStheme
RHPlotREY

TempPlotREY <- ggplot()+
  geom_point(data = REYMonthlyT, aes(x = Month,y = meanAirT), fill = REYcolor, color = REYcolor, shape = 21 )+
  geom_line(data = REYMonthlyT, aes(x = Month,y = meanAirT), color = REYcolor)+
  
#  geom_point(data = REYMinMonthlyT, aes(x = Month,y = minAirT), color = REYcolor, fill = REYcolor, shape = 25)+
#  geom_line(data = REYMinMonthlyT, aes(x = Month,y = minAirT), color = REYcolor)+
  
 # geom_point(data = REYMaxMonthlyT, aes(x = Month,y = maxAirT), color = REYcolor, fill = REYcolor, shape = 24)+
#  geom_line(data = REYMaxMonthlyT, aes(x = Month,y = maxAirT), color = REYcolor)+
  
  scale_x_discrete(limits=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"))+
  scale_y_continuous(limits = c(-6,32), expand = c(0, 0), breaks=seq(-5,30,5)) +
  
  ylab("Air Temp (°C)")+
  xlab("Month")+
  MStheme

TempPlotREY

#load in precipitation data from RCEW, averaged to be daily in the code Drever_Reynolds_Elena_3_18_24.R, starts in 2015 
REYppt <- read.csv("~/Documents/Soil Water 17O Manuscript/data/climo data/Reynolds climo data/RCEW_098c_daily_ppta.csv", sep = ",")
REYppt$Month <- month(REYppt$datetime)
REYppt$Year <- year(REYppt$datetime)
unique(REYppt$Year)

REYMonthlyP <- REYppt %>%
  mutate(Month = month(datetime)) %>%
  group_by(Month) %>%
  summarise(monthlyP = sum(ppta/9, na.rm = TRUE)) %>% #divide by 9 because there are 9 years in the record, otherwise its the sum of all Jans for example
  na.omit() #not sure why this is needed, removes a month that is all NAs

REYAnnualP <- REYppt%>%
  group_by(Year) %>%
  #filter(!is.na(AirTC_Avg_avg)) %>%
  summarise(MAP = sum(ppta))  #can use mean because the monthly precip is already summed

mean(REYAnnualP$MAP, na.rm = TRUE)
sd(REYAnnualP$MAP, na.rm = TRUE)

MAPPlotREY <- ggplot()+
  geom_bar(data = REYMonthlyP, aes(x = Month,y = monthlyP), stat = "identity", color = REYcolor, fill = REYcolor)+
  scale_x_discrete(limits=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"), labels = NULL)+
  scale_y_continuous(limits = c(0,50), expand = c(0, 0), breaks=seq(0,50,10)) +
  xlab("")+
  ylab("Monthly Precip (mm)")+
  MStheme

MAPPlotREY

ClimateNormalsREY <- ggarrange(MAPPlotREY, RHPlotREY, TempPlotREY,nrow = 3, ncol = 1)
ggsave(filename = "ClimateNormalsREY.pdf", plot = ClimateNormalsREY, device = cairo_pdf, height=4,width=4, path = path.to.figs.climate.normals)
                    
MAP <- REYppt %>%
  group_by(Year) %>%
  summarise(total = sum(ppta)) %>% #pretty sure this is precip in mm/year
  na.omit()

####Jornada ####
#daily records at CSAND NPP, downloaded from the LTER website using R script knb-lter-jrn.210437048.40.R with the following citation
  #Anderson, J. 2025. Jornada Basin LTER: Wireless meteorological station at NPP C-SAND site: 30-minute summary data: 2013 - ongoing ver 41. 
  #Environmental Data Initiative. https://doi.org/10.6073/pasta/85ca137e355d2b905133d45223650310 (Accessed 2026-01-06).

JORclim <-read.csv(file = "~/Documents/Soil Water 17O Manuscript/data/climo data/Jornada climo data/LTER climate data/knb-lter-jrn.210437048.40.csv")

#2013 record starts in Nov, and the 2024 record is only a few days in January, remove those years
JORclim <- JORclim %>% filter(Year != 2013) %>% filter(Year != 2024)
unique(JORclim$Year)

JORMonthlyT <- JORclim %>%
  mutate(Month = month(Date)) %>%
  group_by(Month) %>%
  summarise(meanAirT = mean(Air_TempC_Avg, na.rm = TRUE)) 

mean(JORclim$Air_TempC_Avg, na.rm = TRUE) #mean over many years

JORMinMonthlyT <- JORclim %>%
  mutate(Month = month(Date)) %>%
  group_by(Month) %>%
  summarise(minAirT = min(Air_TempC_Min, na.rm = TRUE)) 

JORMaxMonthlyT <- JORclim %>%
  mutate(Month = month(Date)) %>%
  group_by(Month) %>%
  summarise(maxAirT = max(Air_TempC_Max, na.rm = TRUE)) 

TempPlotJOR <- ggplot()+
  geom_point(data = JORMonthlyT, aes(x = Month,y = meanAirT), fill = JORcolor, color = JORcolor, shape = 21 )+
  geom_line(data = JORMonthlyT, aes(x = Month,y = meanAirT), color = JORcolor)+
  
 # geom_point(data = JORMinMonthlyT, aes(x = Month,y = minAirT), color = JORcolor, fill = JORcolor, shape = 25)+
#  geom_line(data = JORMinMonthlyT, aes(x = Month,y = minAirT), color = JORcolor)+
  
#  geom_point(data = JORMaxMonthlyT, aes(x = Month,y = maxAirT), color = JORcolor, fill = JORcolor, shape = 24)+
#  geom_line(data = JORMaxMonthlyT, aes(x = Month,y = maxAirT), color = JORcolor)+
  
  scale_x_discrete(limits=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"))+
  scale_y_continuous(limits = c(-6,32), expand = c(0, 0), breaks=seq(-5,30,5)) +
  
  ylab("Air Temp (°C)")+
  xlab("Month")+
  MStheme

TempPlotJOR

JORMonthlyRH <- JORclim %>%
  mutate(Month = month(Date)) %>%
  group_by(Month) %>%
  summarise(meanRH = mean(Relative_Humidity_Avg, na.rm = TRUE))
mean(JORclim$Relative_Humidity_Avg, na.rm = TRUE) #mean over many years

RHPlotJOR <- ggplot()+
  geom_point(data = JORMonthlyRH, aes(x = Month, y = meanRH), fill = JORcolor, color = JORcolor, shape = 21)+
  geom_line(data = JORMonthlyRH, aes(x = Month, y = meanRH),color = JORcolor)+
  scale_x_discrete(limits=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"), labels = NULL)+
  scale_y_continuous(limits = c(10,80), expand = c(0, 0), breaks=seq(0,80,10)) +
  xlab("")+
  ylab("RH (%)")+
  MStheme
RHPlotJOR

JORMonthlyP <- JORclim %>%
  mutate(Month = month(Date)) %>%
  group_by(Month) %>%
  summarise(monthlyP = sum(Ppt_mm_Tot/10, na.rm = TRUE)) %>% #divide by 10 because there are 10 years in the record, otherwise its the sum of all Jans for example
  na.omit() #not sure why this is needed, removes a month that is all NAs

JORYearlyP <- JORclim %>%
  group_by(Year) %>%
  summarise(YearP = sum(Ppt_mm_Tot, na.rm = TRUE)) %>% #divide by 10 because there are 10 years in the record, otherwise its the sum of all Jans for example
  na.omit() #not sure why this is needed, removes a month that is all NAs

mean(JORYearlyP$YearP) #this is teh mean annual precip for the years of 2014-2023
sum(JORMonthlyP$monthlyP) #this should be the same as the MAP

mean(REYAnnualP$MAP, na.rm = TRUE)
sd(REYAnnualP$MAP, na.rm = TRUE)


MAPPlotJOR <- ggplot()+
  geom_bar(data = JORMonthlyP, aes(x = Month,y = monthlyP), stat = "identity", color = JORcolor, fill = JORcolor)+
  scale_x_discrete(limits=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"), labels = NULL)+
  scale_y_continuous(limits = c(0,50), expand = c(0, 0), breaks=seq(0,50,10)) +
  xlab("")+
  ylab("Monthly Precip (mm)")+
  MStheme
MAPPlotJOR

ClimateNormalsJOR <- ggarrange(MAPPlotJOR, RHPlotJOR, TempPlotJOR,nrow = 3, ncol = 1)
ggsave(filename = "ClimateNormalsJOR.pdf", plot = ClimateNormalsJOR, device = cairo_pdf, height=4,width=4, path = path.to.figs.climate.normals)


##### Southern Michigan####
#downloaded temp and precip normals 1981-2010 from https://www.ncdc.noaa.gov/cdo-web/ as a pdf, then made simple csv files with the data
#using Lat: 42.3263° N Lon: 84.0133° W
#relative humidity took a little more leg work, downloaded hourly data from the same SE AA weather station but from a different site

AAtemp_normals <- read.csv("~/Documents/Soil Water 17O Manuscript/data/climo data/ESGR climo data/AnnArborTemperatureNormals.csv")
AAtemp_normals$Mean_C <- (AAtemp_normals$Mean_F-32)*5/9
AAtemp_normals$DailyMax_C <- (AAtemp_normals$DailyMax_F-32)*5/9
AAtemp_normals$DailyMin_C <- (AAtemp_normals$DailyMin_F-32)*5/9

mean(AAtemp_normals$Mean_C) #mean annual temperature

TempPlotESGR <- ggplot()+
  geom_point(data = AAtemp_normals, aes(x = Month,y = Mean_C), fill = ESGRcolor, color = ESGRcolor, shape = 21 )+
  geom_line(data = AAtemp_normals, aes(x = Month,y = Mean_C), color = ESGRcolor)+
  
#  geom_point(data = AAtemp_normals, aes(x = Month,y = DailyMin_C), color = ESGRcolor, fill = ESGRcolor, shape = 25)+
 # geom_line(data = AAtemp_normals, aes(x = Month,y = DailyMin_C), color = ESGRcolor)+
  
#  geom_point(data = AAtemp_normals, aes(x = Month,y = DailyMax_C), color = ESGRcolor, fill = ESGRcolor, shape = 24)+
#  geom_line(data = AAtemp_normals, aes(x = Month,y = DailyMax_C), color = ESGRcolor)+
  
  scale_x_discrete(limits=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"))+
  scale_y_continuous(limits = c(-6,32), expand = c(0, 0), breaks=seq(-5,30,5)) +
  
  ylab("Air Temp (°C)")+
  xlab("Month")+
  MStheme

TempPlotESGR

AAprecip_normals <- read.csv("~/Documents/Soil Water 17O Manuscript/data/climo data/ESGR climo data/AnnArborPrecipNormals.csv")
AAprecip_normals$MeanPrecip_mm <- AAprecip_normals$MeanPrecip_inches*25.4

MAPPlotESGR <- ggplot()+
  geom_bar(data = AAprecip_normals, aes(x = Month,y = MeanPrecip_mm), stat = "identity", color = ESGRcolor, fill = ESGRcolor)+
  scale_x_discrete(limits=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"), labels = NULL)+
  xlab("")+
  ylab("Monthly Precip (mm)")+
  scale_y_continuous(limits = c(0,91), expand = c(0, 0), breaks=seq(0,90,10)) +
  
  MStheme
MAPPlotESGR

sum(AAprecip_normals$MeanPrecip_mm) #approx MAP in mm

#read in file with RH and other data, hourly data from 2015 to 2024
AAhourly <- read.csv(file = "~/Documents/Soil Water 17O Manuscript/data/climo data/ESGR climo data/3920494.csv")
AAhourly <- AAhourly %>% filter(REPORT_TYPE == "FM-15") #remove teh SOD (summary of date) rows
AAhourly$Date <- date(AAhourly$DATE)
AAhourly$Month <- month(AAhourly$Date)
AAhourly$HourlyRelativeHumidity <- as.numeric(AAhourly$HourlyRelativeHumidity)

ggplot()+
  geom_line(data = AAhourly, aes(x = Date, y = HourlyRelativeHumidity))

#this is VERY slow
AAmonthlyRH <- AAhourly %>%
  group_by(Month) %>%
  summarise(meanRH = mean(HourlyRelativeHumidity, na.rm = TRUE))

RHPlotESGR <- ggplot()+
  geom_point(data = AAmonthlyRH, aes(x = Month, y = meanRH), fill = ESGRcolor, color = ESGRcolor, shape = 21)+
  geom_line(data = AAmonthlyRH, aes(x = Month, y = meanRH), color = ESGRcolor)+
  scale_x_discrete(limits=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"), labels = NULL)+
  scale_y_continuous(limits = c(10,80), expand = c(0, 0), breaks=seq(0,80,10)) +
  xlab("")+
  ylab("RH (%)")+
  MStheme
RHPlotESGR

mean(AAmonthlyRH$meanRH)

ClimateNormalsESGR <- ggarrange(MAPPlotESGR, RHPlotESGR, TempPlotESGR,nrow = 3, ncol = 1)
ggsave(filename = "ClimateNormalsESGR.pdf", plot = ClimateNormalsESGR, device = cairo_pdf, height=4,width=4, path = path.to.figs.climate.normals)

###########combine everything############

ClimateNormalsAll <- ggarrange(ClimateNormalsCRS, ClimateNormalsJOR, ClimateNormalsREY, ClimateNormalsESGR,
                               nrow = 2, ncol = 2)

ggsave(filename = "ClimateNormalsAll.pdf", plot = ClimateNormalsAll, device = cairo_pdf, height=8,width=8, path = path.to.figs.climate.normals)

#try plotting climate variables all together

TempPlotAll <- ggplot()+
  geom_point(data = AAtemp_normals, aes(x = Month,y = Mean_C), fill = ESGRcolor, color = ESGRcolor, shape = 21 )+
  geom_line(data = AAtemp_normals, aes(x = Month,y = Mean_C), color = ESGRcolor)+
  
  geom_point(data = REYMonthlyT, aes(x = Month,y = meanAirT), fill = REYcolor, color = REYcolor, shape = 21 )+
  geom_line(data = REYMonthlyT, aes(x = Month,y = meanAirT), color = REYcolor)+
  
  geom_point(data = CRSMonthlyT, aes(x = Month,y = meanAirT), fill = MOJcolor, color = MOJcolor, shape = 21 )+
  geom_line(data = CRSMonthlyT, aes(x = Month,y = meanAirT), color = MOJcolor)+
  
  geom_point(data = JORMonthlyT, aes(x = Month,y = meanAirT), fill = JORcolor, color = JORcolor, shape = 21 )+
  geom_line(data = JORMonthlyT, aes(x = Month,y = meanAirT), color = JORcolor)+
  
  scale_x_discrete(limits=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"))+
  scale_y_continuous(limits = c(-6,32), expand = c(0, 0), breaks=seq(-5,30,5)) +
  
  ylab("Air Temp (°C)")+
  xlab("Month")+
  MStheme

TempPlotAll

RHPlotAll <- ggplot()+
  geom_point(data = CRSMonthlyRH, aes(x = Month, y = meanRH), fill = MOJcolor, color = MOJcolor, shape = 23)+
  geom_line(data = CRSMonthlyRH, aes(x = Month, y = meanRH), color = MOJcolor)+
  
  geom_point(data = JORMonthlyRH, aes(x = Month, y = meanRH), fill = JORcolor, color = JORcolor, shape = 23)+
  geom_line(data = JORMonthlyRH, aes(x = Month, y = meanRH),color = JORcolor)+
  
  geom_point(data = REYMonthlyRH, aes(x = Month, y = meanRH), fill = REYcolor, color = REYcolor, shape = 23)+
  geom_line(data = REYMonthlyRH, aes(x = Month, y = meanRH), color = REYcolor)+
  
  geom_point(data = AAmonthlyRH, aes(x = Month, y = meanRH), fill = ESGRcolor, color = ESGRcolor, shape = 23)+
  geom_line(data = AAmonthlyRH, aes(x = Month, y = meanRH), color = ESGRcolor)+
  
  scale_x_discrete(limits=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"), labels = NULL)+
  scale_y_continuous(limits = c(10,80), expand = c(0, 0), breaks=seq(0,80,10)) +
  xlab("")+
  ylab("RH (%)")+
  MStheme
RHPlotAll

MAPPlotAll <- ggplot()+
  geom_point(data = CRSMonthlyP, aes(x = Month,y = monthlyP),color = MOJcolor, fill = MOJcolor, shape = 21)+  
  geom_line(data = CRSMonthlyP, aes(x = Month,y = monthlyP),color = MOJcolor)+
  
  geom_point(data = REYMonthlyP, aes(x = Month,y = monthlyP),color = REYcolor, fill = REYcolor, shape = 21)+  
  geom_line(data = REYMonthlyP, aes(x = Month,y = monthlyP),color = REYcolor)+
  
  geom_point(data = JORMonthlyP, aes(x = Month,y = monthlyP),color = JORcolor, fill = JORcolor, shape = 21)+  
  geom_line(data = JORMonthlyP, aes(x = Month,y = monthlyP),color = JORcolor)+
  
  geom_point(data = AAprecip_normals, aes(x = Month,y = MeanPrecip_mm), color = ESGRcolor, fill = ESGRcolor, shape  = 21)+
  geom_line(data = AAprecip_normals, aes(x = Month,y = MeanPrecip_mm),color = ESGRcolor)+

  scale_x_discrete(limits=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"), labels = NULL)+
  scale_y_continuous(limits = c(0,92), expand = c(0, 0), breaks=seq(0,90,10)) +
  xlab("")+
  ylab("Monthly Precip (mm)")+
  MStheme

MAPPlotAll


ClimateNormalsAll2 <- ggarrange(MAPPlotAll, RHPlotAll, TempPlotAll,
                               nrow = 3, ncol = 1)

ggsave(filename = "ClimateNormalsAll2.pdf", plot = ClimateNormalsAll2, device = cairo_pdf, height=8,width=8, path = path.to.figs.climate.normals)


