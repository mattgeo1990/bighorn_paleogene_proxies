#plotting CZ17O  WATERS data into a summary panel (Fig 3) and depth profiles for each site (Figs 5,6,7,8)
#also calculates lambda_soil and lambda_met value
#includes d18O and d2H data from Picarro/CRDS
#and triple O data from IRMS
#starts with mw.csv and sw.csv, files that have combined isotope data and metadata

library(ggplot2)  
library(ggpubr)
library(scales)
library(dplyr)
library(forcats)
library(lubridate)
library(readxl)
library(plotly)
library(stringr)

########STEP 1: read in mw and sw data, merged picarro and IRMS data ###########

mw <- read.csv("~/Documents/Soil Water 17O Manuscript/data/mw.csv")
sw <- read.csv("~/Documents/Soil Water 17O Manuscript/data/sw.csv")

###### MeanWeightedP values calculated in other code#####
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
ESGRmetw <- mw %>% filter(siteID.1 == "ESGR") 

precip <- mw %>% filter(Water.Type == "Precipitation")

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


#####identify duplicates in the sw dataset#####
swdups <- sw$SampleIDLong[duplicated(sw$SampleIDLong)]
#this only says TRUE for the second duplicate

swduplicates <- sw %>%
  group_by(SampleIDLong) %>%
  summarise(d18O.mean.dups = mean(d18O.mean), d18O.sd = sd(d18O.mean), d18O_diff = abs(d18O.mean[1]-d18O.mean[2]), total_count = n()) %>%
  filter(total_count>1)

max(abs(swduplicates$d18O_diff))
min(abs(swduplicates$d18O_diff))
median(abs(swduplicates$d18O_diff))

mean(abs(swduplicates$d18O_diff))


#####calculate correlations in precipitation values#######
heavyprecip <- precip %>% filter(d18O.mean > -5)
line1 <- lm(heavyprecip$d_excess.mean~ heavyprecip$d18O.mean)
summary(line1)
line2 <- lm(heavyprecip$D17O_pmg~ heavyprecip$d18O.mean)
summary(line2)

lightprecip <- precip %>% filter(d18O.mean <= -5)
line1 <- lm(lightprecip$d_excess.mean~ lightprecip$d18O.mean)
summary(line1)
line2 <- lm(lightprecip$D17O_pmg~ lightprecip$d18O.mean)
summary(line2)

line1 <- lm(precip$d_excess.mean~ precip$d18O.mean)
summary(line1)
line2 <- lm(precip$D17O_pmg~ precip$d18O.mean)
summary(line2)
line3 <- lm(precip$d_excess.mean~ precip$D17O_pmg)
summary(line3)
line4 <- lm(mw$d_excess.mean~ mw$D17O_pmg)
summary(line4)

line1 <- lm(sw$d_excess.mean~ sw$d18O.mean)
summary(line1)
line2 <- lm(sw$D17O_pmg~ sw$d18O.mean)
summary(line2)
line3 <- lm(sw$d_excess.mean~ sw$D17O_pmg)
summary(line3)
line4 <- lm(sw$D17O_pmg~ sw$d_excess.mean)
summary(line4)

#ggplot(data = sw, aes(x = d_excess.mean, y= D17O_pmg ))+
#  geom_point(fill = "orange",shape = 21, size = 3)+
#  geom_smooth( method='lm')+
#  geom_abline(slope = line4$coefficients[2], intercept = line4$coefficients[1])


###########seasonality analysis for supporting information figure############
#average analyses
mean(precip$D17O_pmg, na.rm = TRUE)

precip$month <- month(precip$date.collection)
precip$season[precip$month == 6 | precip$month == 7 | precip$month == 8 ] <- "summer"
precip$season[precip$month == 9 | precip$month == 10 | precip$month == 11 ] <- "fall"
precip$season[precip$month == 12 | precip$month == 1 | precip$month == 2 ] <- "winter"
precip$season[precip$month == 3 | precip$month == 4 | precip$month == 5 ] <- "spring"

sw$month <- month(sw$date.collection)
sw$season[sw$month == 6 | sw$month == 7 | sw$month == 8 ] <- "summer"
sw$season[sw$month == 9 | sw$month == 10 | sw$month == 11 ] <- "fall"
sw$season[sw$month == 12 | sw$month == 1 | sw$month == 2 ] <- "winter"
sw$season[sw$month == 3 | sw$month == 4 | sw$month == 5 ] <- "spring"

t.test(x = precip$D17O_pmg[precip$season == "summer"], y = precip$D17O_pmg[precip$season == "winter"])
t.test(x = precip$D17O_pmg[precip$season == "winter"], y = precip$D17O_pmg[precip$season == "summer"])

t.test(x = precip$d18O.mean[precip$season == "summer"], y = precip$d18O.mean[precip$season == "winter"])

t.test(x = precip$d_excess.mean[precip$season == "summer"], y = precip$d_excess.mean[precip$season == "winter"])

SEAmet1 <- ggplot()+
  geom_point(data = precip, aes(x = d18O.mean, y = D17O_pmg, fill = season, shape = siteID.1), size = 3)+
  scale_fill_manual(values = c("tomato", "purple", "chartreuse", "dodgerblue" ))+
  scale_shape_manual(values = c(21, 22, 23, 24))+
  xlab(expression(paste(delta^{18}, "O (‰, SMOW)")))+
  ylab(expression(Delta^"'17"*"O (pmg, SMOW-SLAP)"))+
    #scale_fill_manual(values = c("tomato", "purple", "chartreuse", "dodgerblue" ))+
  MStheme_isos
SEAmet1

SEAmet2 <- ggplot()+
  geom_point(data = precip, aes(x = d18O.mean, y = d_excess.mean, fill = season, shape = siteID.1), size = 3)+
  scale_fill_manual(values = c("tomato", "purple", "chartreuse", "dodgerblue" ))+
  scale_shape_manual(values = c(21, 22, 23, 24))+
  xlab(expression(paste(delta^{18}, "O (‰, SMOW)")))+
  ylab(expression(paste(italic("d"),"-excess (‰, SMOW)"))) +
    #scale_fill_manual(values = c("tomato", "purple", "chartreuse", "dodgerblue" ))+
  MStheme_isos

SEAmet2

SEApart1 <- ggarrange(SEAmet1, SEAmet2, nrow = 1, ncol = 2, common.legend = TRUE)

SEASdataforboxes <- data.frame(d18O.mean = c(precip$d18O.mean, sw$d18O.mean), d2H.mean = c(precip$d2H.mean, sw$d2H.mean), 
                              d_excess.mean  = c(precip$d_excess.mean, sw$d_excess.mean), 
                              D17O_pmg = c(precip$D17O_pmg, sw$D17O_pmg), season = c(precip$season, sw$season) )

counts <- SEASdataforboxes %>% count(season)
print(counts)

counts_D17O <- SEASdataforboxes %>% filter(!is.na(D17O_pmg)) %>% count(season)
print(counts_D17O)


SEASboxd18O <- ggplot()+
  geom_boxplot(data = precip, aes(x = season, y = d18O.mean, color = season), outlier.shape = NA)+
  scale_color_manual(values = c("tomato", "purple", "chartreuse", "dodgerblue" ))+
  geom_jitter(data = precip, aes(x = season, y = d18O.mean, shape = siteID.1, color = season, fill = season), size = 2)+
  scale_shape_manual(values = c(21, 22, 23, 24))+
  scale_fill_manual(values = c("tomato", "purple", "chartreuse", "dodgerblue" ))+
  
  #scale_y_continuous(limits = c(-20,15), expand = c(0,0), breaks = seq(-20, 15, 5))+
  xlab("")+
  ylab(expression(paste(delta^{18}, "O (‰, SMOW)")))+
  MStheme_isos #+ theme(axis.text.y = element_blank() )
SEASboxd18O

SEASboxdxs <- ggplot()+
  geom_boxplot(data = precip, aes(x = season, y = d_excess.mean, color = season), outlier.shape = NA)+
  scale_color_manual(values = c("tomato", "purple", "chartreuse", "dodgerblue" ))+
  
  geom_jitter(data = precip, aes(x = season, y = d_excess.mean, shape = siteID.1, color = season, fill = season), size = 2)+
  scale_shape_manual(values = c(21, 22, 23, 24))+
  scale_fill_manual(values = c("tomato", "purple", "chartreuse", "dodgerblue" ))+
  #scale_y_continuous(limits = c(-20,15), expand = c(0,0), breaks = seq(-20, 15, 5))+
  xlab("")+
  ylab(expression(paste(italic("d"),"-excess (‰, SMOW)"))) +
  MStheme_isos #+ theme(axis.text.y = element_blank() )

SEASboxdxs

SEASboxD17O <- ggplot()+
  geom_boxplot(data = precip, aes(x = season, y = D17O_pmg, color = season), outlier.shape = NA)+
  scale_color_manual(values = c("tomato", "purple", "chartreuse", "dodgerblue" ))+
  
  geom_jitter(data = precip, aes(x = season, y = D17O_pmg, shape = siteID.1, color = season, fill = season), size = 2)+
  scale_shape_manual(values = c(21, 22, 23, 24))+
  scale_fill_manual(values = c("tomato", "purple", "chartreuse", "dodgerblue" ))+
  #scale_y_continuous(limits = c(-20,15), expand = c(0,0), breaks = seq(-20, 15, 5))+
  xlab("")+
  ylab(expression(Delta^"'17"*"O (pmg, SMOW-SLAP)"))+
  MStheme_isos #+ theme(axis.text.y = element_blank() )
SEASboxD17O

SEAboxes <- ggarrange(SEASboxd18O, SEASboxdxs, SEASboxD17O, nrow = 3, ncol = 1, common.legend = TRUE)
SEAboxes

SEAtogether <- ggarrange(SEApart1, SEAboxes, ncol = 1, nrow = 2)
SEAtogether
ggsave(filename = "SeasonalitySummary.pdf", plot = SEAtogether, device = cairo_pdf, height=10,width=8, path = path.to.figs.refined)
#needs lots of clean up and adjusting of aspect ratios


############ loads to BA1983 modeling ############
BAdf <- read.csv("/Users/juliakelson/Documents/TripleO/soils/BarnesAllisonFigs/ForMSv8/BA1983_baserun.csv")
BAdrydf <- read.csv("/Users/juliakelson/Documents/TripleO/soils/BarnesAllisonFigs/ForMSv8/BA1983_dryrun.csv")
BAhumiddf <- read.csv("/Users/juliakelson/Documents/TripleO/soils/BarnesAllisonFigs/ForMSv8/BA1983_humidrun.csv")

BAdfMOJ <- read.csv("/Users/juliakelson/Documents/CZ17O_R/figs/Gold_figs/BA1983_baserunMOJ.csv")

BAdf$d.excess <- BAdf$dDsw -8 *BAdf$d18Osw
BAdrydf$d.excess <- BAdrydf$dDsw -8 *BAdrydf$d18Osw
BAhumiddf$d.excess <- BAhumiddf$dDsw -8 *BAhumiddf$d18Osw
BAdfMOJ$d.excess <- BAdfMOJ$dDsw -8 *BAdfMOJ$d18Osw

##BA trends

BAdf$dp18Osw <- log(BAdf$d18Osw/1000+1)*1000
BAdf$dp17Osw <- (BAdf$D17Osw/1000)+0.528*BAdf$dp18Osw
SWL_BAdf <- lm(BAdf$dDsw ~ BAdf$d18Osw)
summary(SWL_BAdf)
lambda_empirical_BAdf <- lm(BAdf$dp17Osw~BAdf$dp18Osw)
summary(lambda_empirical_BAdf)

BAdrydf$dp18Osw <- log(BAdrydf$d18Osw/1000+1)*1000
BAdrydf$dp17Osw <- (BAdrydf$D17Osw/1000)+0.528*BAdrydf$dp18Osw
SWL_BAdrydf <- lm(BAdrydf$dDsw ~ BAdrydf$d18Osw)
summary(SWL_BAdrydf)
lambda_empirical_BAdrydf <- lm(BAdrydf$dp17Osw~BAdrydf$dp18Osw)
summary(lambda_empirical_BAdrydf)

BAhumiddf$dp18Osw <- log(BAhumiddf$d18Osw/1000+1)*1000
BAhumiddf$dp17Osw <- (BAhumiddf$D17Osw/1000)+0.528*BAhumiddf$dp18Osw
lambda_empirical_BAhumiddf <- lm(BAhumiddf$dp17Osw~BAhumiddf$dp18Osw)
summary(lambda_empirical_BAhumiddf)
SWL_BAhumiddf <- lm(BAhumiddf$dDsw ~ BAhumiddf$d18Osw)
summary(SWL_BAhumiddf)

#model comparison in isotope space
ggplot()+
  geom_point(data = BAdf, aes(x = d18Osw, y = dDsw), size = 2, shape = 21, color = "black")+
  geom_point(data = BAdrydf, aes(x = d18Osw, y = dDsw), size = 2, shape = 21, color = "red")+
  geom_point(data = BAhumiddf, aes(x = d18Osw, y = dDsw), size = 2, shape = 21, color = "blue")+
  geom_point(data = BAdfMOJ, aes(x = d18Osw, y = dDsw), size = 2, shape = 21, color = "gray")

ggplot()+
  geom_point(data = BAdf, aes(x = D17Osw, y = d.excess), size = 2, shape = 21, color = "black")+
  geom_point(data = BAdrydf, aes(x = D17Osw, y = d.excess), size = 2, shape = 21, color = "red")+
  geom_point(data = BAhumiddf, aes(x = D17Osw, y = d.excess), size = 2, shape = 21, color = "blue")+
  geom_point(data = BAdfMOJ, aes(x = D17Osw, y = d.excess), size = 2, shape = 21, color = "gray")+
  geom_point(data  = sw, aes(x = D17O_pmg, y = d_excess.mean), size = 2, shape = 22, fill = "black")+
  geom_point(data  = mw, aes(x = D17O_pmg, y = d_excess.mean), size = 2, shape = 22, fill = "purple")



######### explore soil waters, calculate lambda values#####

##calculate the d'18O'd'17O relationship (lambda) for soil waters
lambda_empirical_sw <- lm(sw$dp17O~sw$dp18O) 
summary(lambda_empirical_sw)

lambda_empirical_swREY <- lm(REYsw$dp17O~REYsw$dp18O)
lambda_empirical_swJOR <- lm(JORsw$dp17O~JORsw$dp18O)
lambda_empirical_swMOJ <- lm(MOJsw$dp17O~MOJsw$dp18O)
lambda_empirical_swESGR<- lm(ESGRsw$dp17O~ESGRsw$dp18O)

summary(lambda_empirical_swREY)
summary(lambda_empirical_swJOR)
summary(lambda_empirical_swMOJ)
summary(lambda_empirical_swESGR)

#break MOJ down by site, not used
CRSsw <- MOJsw %>% filter(siteID.2 == "CRS")
JTsw <- MOJsw %>% filter(siteID.2 == "JT")
PJsw <- MOJsw %>% filter(siteID.2 == "PJ")

lambda_empirical_swCRS <- lm(CRSsw$dp17O~CRSsw$dp18O)
lambda_empirical_swJT <- lm(JTsw$dp17O~JTsw$dp18O)
lambda_empirical_swPJ <- lm(PJsw$dp17O~PJsw$dp18O)
summary(lambda_empirical_swCRS)
summary(lambda_empirical_swJT)
summary(lambda_empirical_swPJ)

trendALL <- lm(sw$D17O_pmg ~sw$d18O_IRMS)
summary(trendALL)

##calculate the d'18O'd'17O relationship (lambda) for  ALL met waters
lambda_empirical_mw <- lm(mw$dp17O~mw$dp18O)
summary(lambda_empirical_mw)

lambda_empirical_mwREY <- lm(REYmetw$dp17O~REYmetw$dp18O)
lambda_empirical_mwJOR <- lm(JORmetw$dp17O~JORmetw$dp18O)
lambda_empirical_mwMOJ <- lm(MOJmetw$dp17O~MOJmetw$dp18O)
lambda_empirical_mwESGR<- lm(ESGRmetw$dp17O~ESGRmetw$dp18O)

summary(lambda_empirical_mwREY)
summary(lambda_empirical_mwJOR)
summary(lambda_empirical_mwMOJ)
summary(lambda_empirical_mwESGR)

trendALL <- lm(sw$D17O_pmg ~sw$d18O_IRMS)
summary(trendALL)

##calculate the d18O-d2H relationships for "unevap" met waters (LMWLs)
mw.unevap <- mw %>% filter(Water.Type == "Marsh" | Water.Type == "Creek" |Water.Type == "River" | Water.Type == "Well" | Water.Type == "Precipitation" |
                             Water.Type == "spring" | Water.Type == "Spring") %>%
                  filter(d18O.mean<10)
empirical_LMWL <- lm(mw.unevap$d2H.mean~mw.unevap$d18O.mean)
summary(empirical_LMWL)

empirical_LMWL_REY <- lm(REYmetw$d2H.mean~REYmetw$d18O.mean)
empirical_LMWL_JOR <- lm(JORmetw$d2H.mean~JORmetw$d18O.mean)
empirical_LMWL_JOR_2 <- lm(JORmetw$d2H.mean[JORmetw$d18O.mean<10]~JORmetw$d18O.mean[JORmetw$d18O.mean<10])
empirical_LMWL_MOJ <- lm(MOJmetw$d2H.mean[MOJmetw$Water.Type != "Puddle"]~MOJmetw$d18O.mean[MOJmetw$Water.Type != "Puddle"])
empirical_LMWL_ESGR <- lm(ESGRmetw$d2H.mean~ESGRmetw$d18O.mean) #INCOMPLETE

summary(empirical_LMWL_REY)
summary(empirical_LMWL_JOR)
summary(empirical_LMWL_JOR_2)
summary(empirical_LMWL_MOJ)
summary(empirical_LMWL_ESGR)

#LMWL: precipitation only
mw.precip <- mw %>% filter(Water.Type == "Precipitation")
REYprecip <- mw.precip %>% filter(siteID.1 == "REY")
MOJprecip <- mw.precip %>% filter(siteID.1 == "MOJ")
JORprecip <- mw.precip %>% filter(siteID.1 == "JOR")
ESGRprecip <- mw.precip %>% filter(siteID.1 == "ESGR")

mean(mw.precip$D17O_pmg, na.rm = TRUE)
mean(ESGRprecip$D17O_pmg, na.rm = TRUE)
mean(mw.precip$D17O_pmg[mw.precip$d18O.mean<0], na.rm = TRUE)


MWL_precip <- lm(mw.precip$d2H.mean~mw.precip$d18O.mean) #this doesn'y really mean much
summary(MWL_precip)

empirical_LMWL_precip_REY <- lm(REYprecip$d2H.mean~REYprecip$d18O.mean)
empirical_LMWL_precip_JOR <- lm(JORprecip$d2H.mean~JORprecip$d18O.mean)
empirical_LMWL_precip_JOR_2 <- lm(JORprecip$d2H.mean[JORprecip$d18O.mean<10]~JORprecip$d18O.mean[JORprecip$d18O.mean<10])
empirical_LMWL_precip_MOJ <- lm(MOJprecip$d2H.mean~MOJprecip$d18O.mean)
empirical_LMWL_precip_ESGR <- lm(ESGRprecip$d2H.mean~ESGRprecip$d18O.mean) #INCOMPLETE

summary(empirical_LMWL_precip_REY)
summary(empirical_LMWL_precip_JOR)
summary(empirical_LMWL_precip_JOR_2)
summary(empirical_LMWL_precip_MOJ)
summary(empirical_LMWL_precip_ESGR)


lambda_empirical_precipREY <- lm(REYprecip$dp17O~REYprecip$dp18O)
lambda_empirical_precipJOR <- lm(JORprecip$dp17O~JORprecip$dp18O)
#lambda_empirical_precipJOR_2 <- lm(JORprecip$dp17O[JORprecip$d18O.mean<10]~JORprecip$dp18O[JORprecip$d18O.mean<10]) #this just removes a sample with NA
lambda_empirical_precipMOJ <- lm(MOJprecip$dp17O~MOJprecip$dp18O)
lambda_empirical_precipESGR<- lm(ESGRprecip$dp17O~ESGRprecip$dp18O)

summary(lambda_empirical_precipREY)
summary(lambda_empirical_precipJOR)
#summary(lambda_empirical_precipJOR_2)
summary(lambda_empirical_precipMOJ)
summary(lambda_empirical_precipESGR)

lambda_empirical_precipALL <- lm(mw.precip$dp17O~mw.precip$dp18O)
summary(lambda_empirical_precipALL)

#MWLplot <- ggplot(data = mw.precip)+
#  geom_point(aes(x = d18O.mean, y = d2H.mean, color = siteID.1))+
#  geom_abline(slope = MWL_precip$coefficients[2], intercept = MWL_precip$coefficients[1])
#ggplotly(MWLplot)

#empirical SWLS
empirical_SWL <- lm(sw$d2H.mean~sw$d18O.mean) #this doesn'y really mean much
summary(empirical_SWL)

empirical_SWL_REY <- lm(REYsw$d2H.mean~REYsw$d18O.mean)
empirical_SWL_JOR <- lm(JORsw$d2H.mean~JORsw$d18O.mean)
empirical_SWL_MOJ <- lm(MOJsw$d2H.mean~MOJsw$d18O.mean)
empirical_SWL_ESGR <- lm(ESGRsw$d2H.mean~ESGRsw$d18O.mean)

summary(empirical_SWL_REY)
summary(empirical_SWL_JOR)
summary(empirical_SWL_MOJ)
summary(empirical_SWL_ESGR)




################# Reynolds ###########
#citation for published reynolds data:
#Lohse, Kathleen; Schegel, Melissa; Souza, Jennifer; Warix, Sara; MacNeille, Ruth; Murray, Erin; Radke, Anna; Godsey, Sarah E.; Seyfried, Mark S.; and Finney, Bruce. (2022). Dataset for Stable Isotopes of Precipitation, Surface Water, Spring Water, and Subsurface Waters at the Reynolds Creek Experimental Watershed and Critical Zone Observatory, Southwestern Idaho, USA [Data set]. Retrieved from 10.18122/reynoldscreek.29.boisestate
path.to.data <- "~/Documents/Soil Water 17O Manuscript/data/"

RCEWmet <- read.csv(paste(path.to.data, "RCEWWaterIsotopesFINAL_Lohse2022.csv", sep = ""))
RCEWmet.precip <- RCEWmet %>% filter(SubType == "Precipitation")
LMWL.REY.all <- lm(c(REYprecip$d2H.mean,RCEWmet.precip$dD_SMOW ) ~ c(REYprecip$d18O.mean,RCEWmet.precip$d18O_SMOW) )
summary(LMWL.REY.all)

##depth profiles for Reynolds

tryprof <- ggplot(data = REYsw)+
  #geom_path(aes(x = d18O.mean, y = -bottom.depth, group = date.collection), color = REYcolor, size = 0.5)+
  geom_point(aes(x = d18O.mean, y = -bottom.depth, fill = as.factor(date.collection)), size = 4, shape = 21)+
  scale_fill_manual(values = c("tomato", "dodgerblue", "purple", "chartreuse", "tomato3", "dodgerblue4", "purple4", "chartreuse4", "tomato4"))+
  #geom_point(aes(x = mean.weighted.precip.REY.d18O, y = 0), size = 5, shape = 25, fill = "black")+
  geom_vline(xintercept = mean.weighted.precip.REY.d18O, linetype = "longdash")+
  scale_y_continuous(limits = c(-66,1), expand = c(0, 0), breaks=seq(-60,0,10)) +
  scale_x_continuous(limits = c(-18,1), expand = c(0,0), breaks = seq(-15, 0, 5))+
    xlab(expression(paste(delta^{18}, "O (‰, SMOW)")))+
  #ylab("Depth (cm)")+
  ylab("")+
  MStheme_isos +  theme(axis.text.y = element_blank() )

tryprof

tryprof2 <- ggplot(data = REYsw)+
  #geom_path(aes(x = d_excess.mean, y = -bottom.depth, group = date.collection), color = REYcolor, size = 0.5)+
  geom_point(aes(x = d_excess.mean, y = -bottom.depth, fill = as.factor(date.collection)), size = 4, shape = 21)+
  scale_fill_manual(values = c("tomato", "dodgerblue", "purple", "chartreuse", "tomato3", "dodgerblue4", "purple4", "chartreuse4", "tomato4"))+
  
  geom_vline(xintercept = mean.weighted.precip.REY.dexcess, linetype = "longdash")+
  scale_y_continuous(limits = c(-66,1), expand = c(0, 0), breaks=seq(-60,0,10)) +
  scale_x_continuous(limits = c(-90,14), expand = c(0,0), breaks = seq(-90, 10, 10))+
  xlab(expression(paste(italic("d"),"-excess (‰, SMOW)"))) +
  ylab("")+
  MStheme_isos+
  theme(axis.text.y = element_blank() )

tryprof2

tryprof3 <- ggplot(data = REYsw)+  
  #geom_path(aes(x = D17O_pmg, y = -bottom.depth, group = date.collection), color = REYcolor, size = 0.5)+
  geom_errorbarh(aes(xmin = D17O_pmg - D17O_err, xmax = D17O_pmg+ D17O_err, y = -bottom.depth), height = 0)+
  geom_point(aes(x = D17O_pmg, y = -bottom.depth, fill = as.factor(date.collection)), size = 4, shape = 21)+
  scale_fill_manual(values = c("tomato", "dodgerblue", "purple", "chartreuse", "tomato3", "dodgerblue4", "purple4", "chartreuse4", "tomato4"))+
  geom_vline(xintercept = mean(REYprecip$D17O_pmg, na.rm = TRUE), linetype = "longdash")+
  
  scale_y_continuous(limits = c(-66,1), expand = c(0, 0), breaks=seq(-60,0,10)) +
  scale_x_continuous(limits = c(-55,50), expand = c(0,0), breaks = seq(-50, 50, 20))+
  xlab(expression(Delta^"'17"*"O (pmg, SMOW-SLAP)"))+
  ylab("")+
  MStheme_isos+
  theme(axis.text.y = element_blank() )

tryprof3

profsummREY<- ggarrange(tryprof, tryprof2, tryprof3, ncol = 3, legend = "none")

#add box and whisker plots
REYsw$depth_cat <- NA
REYsw$depth_cat[REYsw$bottom.depth <= 5] <- "f__0 to 5"
REYsw$depth_cat[REYsw$bottom.depth > 5 & REYsw$bottom.depth <= 25] <- "e_5 to 25"
REYsw$depth_cat[REYsw$bottom.depth > 25 & REYsw$bottom.depth <= 40] <- "d_25 to 40"
REYsw$depth_cat[REYsw$bottom.depth > 40] <- "c_40 to 65"
REYsw$soilormet <- "soil"

REYmetw$depth_cat[REYmetw$Water.Type == "Precipitation"] <- "b_Precipitation"
REYmetw$depth_cat[REYmetw$Water.Type == "River_or_stream"]  <- "a_River_or_stream"
REYmetw$soilormet <- "met"

REYdataforboxes <- data.frame(d18O.mean = c(REYsw$d18O.mean, REYmetw$d18O.mean), d2H.mean = c(REYsw$d2H.mean, REYmetw$d2H.mean), 
                              d_excess.mean  = c(REYsw$d_excess.mean, REYmetw$d_excess.mean), 
                           D17O_pmg = c(REYsw$D17O_pmg, REYmetw$D17O_pmg), depth_cat = c(REYsw$depth_cat, REYmetw$depth_cat), 
                           soilormet =c(REYsw$soilormet, REYmetw$soilormet) )

counts <- REYdataforboxes %>% count(depth_cat)
print(counts)

counts_D17O <- REYdataforboxes %>% filter(!is.na(D17O_pmg)) %>% count(depth_cat)
print(counts_D17O)

REYboxd18O <- ggplot()+
  geom_boxplot(data = REYdataforboxes, aes(x = depth_cat, y = d18O.mean, fill= soilormet), outlier.shape = NA)+
  geom_jitter(data = REYdataforboxes, aes(x = depth_cat, y = d18O.mean, fill = soilormet, shape = soilormet), color = "black")+
  scale_shape_manual(values = c(22, 21))+
  coord_flip()+
  scale_fill_manual(values = c(REYcolor, REYcolor.minor))+
  #scale_x_discrete(labels = c("River or Stream","Precipitation",   "40-65 cm","25-40 cm","5-25 cm", "0-5 cm"))+
  #scale_x_discrete(labels = c("","",   "","","", ""))+
  
  scale_y_continuous(limits = c(-18,1), expand = c(0,0), breaks = seq(-15, 0, 5))+
  xlab("")+
  ylab(expression(paste(delta^{18}, "O (‰, SMOW)")))+
  MStheme_isos
REYboxd18O

REYboxdxs <- ggplot()+
  geom_boxplot(data = REYdataforboxes, aes(x = depth_cat, y = d_excess.mean, fill= soilormet), outlier.shape = NA)+
  geom_jitter(data = REYdataforboxes, aes(x = depth_cat, y = d_excess.mean, fill = soilormet, shape = soilormet), color = "black")+
  scale_shape_manual(values = c(22, 21))+
  scale_fill_manual(values = c(REYcolor, REYcolor.minor))+
  scale_x_discrete(labels = c("","",   "","","", ""))+
  #scale_x_discrete(labels = c("River or Stream","Precipitation",   "40-65 cm","25-40 cm","5-25 cm", "0-5 cm"))+
  scale_y_continuous(limits = c(-90,14), expand = c(0,0), breaks = seq(-90, 10, 10))+
  coord_flip()+
  xlab("")+
  ylab(expression(paste(italic("d"),"-excess (‰, SMOW)"))) +
  MStheme_isos
REYboxdxs

REYboxD17O <- ggplot()+
  geom_boxplot(data = REYdataforboxes, aes(x = depth_cat, y = D17O_pmg, fill= soilormet), outlier.shape = NA)+
  geom_jitter(data = REYdataforboxes, aes(x = depth_cat, y = D17O_pmg, fill = soilormet, shape = soilormet), color = "black")+
  scale_shape_manual(values = c(22, 21))+
  scale_fill_manual(values = c(REYcolor, REYcolor.minor))+
  scale_x_discrete(labels = c("","",   "","","", ""))+
  
  #scale_x_discrete(labels = c("River or Stream","Precipitation",   "40-65 cm","25-40 cm","5-25 cm", "0-5 cm"))+
  #scale_y_continuous(limits = c(-50,45), expand = c(0,0), breaks = seq(-50, 40, 10))+
  scale_y_continuous(limits = c(-55,50), expand = c(0,0), breaks = seq(-50, 50, 20))+
  
  coord_flip()+
  ylab(expression(Delta^"'17"*"O (pmg, SMOW-SLAP)"))+
  xlab("")+
  MStheme_isos
REYboxD17O

REYboxes <- ggarrange(REYboxd18O, REYboxdxs, REYboxD17O, nrow =1, ncol = 3, legend = "none")
REYboxes

ggsave(filename = "REYboxes.pdf", plot = REYboxes, device = cairo_pdf, height=3,width=12, path = path.to.figs)

REYprofandboxes <- ggarrange(profsummREY, REYboxes, nrow = 2, ncol = 1, heights = c(1, 0.5))
ggsave(filename = "REYprofandboxes.pdf", plot = REYprofandboxes, device = cairo_pdf, height=8.5,width=11, path = path.to.figs)

#calculate means Reynolds
mean(REYprecip$D17O_pmg, na.rm = TRUE) - mean(REYdataforboxes$D17O_pmg[REYdataforboxes$depth_cat == "c_40 to 65"], na.rm = TRUE)
mean.weighted.precip.REY.d18O - mean(REYdataforboxes$d18O.mean[REYdataforboxes$depth_cat == "c_40 to 65"], na.rm = TRUE)
mean.weighted.precip.REY.d2H - mean(REYdataforboxes$d2H.mean[REYdataforboxes$depth_cat == "c_40 to 65"], na.rm = TRUE)
mean.weighted.precip.REY.dexcess- mean(REYdataforboxes$d_excess.mean[REYdataforboxes$depth_cat == "c_40 to 65"], na.rm = TRUE)


######### Mojave ###########
#depth profiles for MOJ########

tryprofCRS <- ggplot(data = MOJsw %>% filter (siteID.2 == "CRS"))+
  #geom_path(aes(x = d18O.mean, y = -bottom.depth, group = date.collection), color = MOJcolor, size = 0.5)+
  geom_point(aes(x = d18O.mean, y = -bottom.depth, fill = as.factor(date.collection)), size = 4, shape = 21)+
  scale_fill_manual(values = c("tomato", "dodgerblue",  "chartreuse", "tomato3", "dodgerblue4", "chartreuse4", "tomato4"))+
  
  #scale_fill_manual(values = c("darkorange1", "deepskyblue", "darkolivegreen1", "darkorange3", "deepskyblue3", "darkolivegreen3", "darkorange4"))+
  geom_vline(xintercept = mean.weighted.precip.MOJ.d18O, linetype = "longdash")+
  scale_y_continuous(limits = c(-95,0), expand = c(0, 0), breaks=seq(-90,0,10)) +
  scale_x_continuous(limits = c(-13,11), expand = c(0,0), breaks = seq(-12, 10, 2))+
  
  xlab(expression(paste(delta^{18}, "O (‰, SMOW)")))+
  #ylab("Depth (cm)")+
  ylab("")+
  MStheme_isos + theme(axis.text.y = element_blank() )
tryprofCRS

tryprof2CRS <- ggplot(data = MOJsw %>% filter (siteID.2 == "CRS"))+
  geom_point(aes(x = d_excess.mean, y = -bottom.depth, fill = as.factor(date.collection)), size = 4, shape = 21)+
  #geom_line(aes(x = d_excess.mean, y = -bottom.depth, group = date.collection), color = MOJcolor, size = 0.5)+
  scale_fill_manual(values = c("tomato", "dodgerblue",  "chartreuse", "tomato3", "dodgerblue4", "chartreuse4", "tomato4"))+
  geom_vline(xintercept = mean.weighted.precip.MOJ.dexcess, linetype = "longdash")+
  scale_x_continuous(limits = c(-125,15), expand = c(0,0), breaks = seq(-120, 10, 20))+
  scale_y_continuous(limits = c(-95,0), expand = c(0, 0), breaks=seq(-90,0,10)) +
  xlab(expression(paste(italic("d"),"-excess (‰, SMOW)"))) +
  ylab("")+
  MStheme_isos+   theme(axis.text.y = element_blank() )

tryprof2CRS

tryprof3CRS <- ggplot(data = MOJsw %>% filter (siteID.2 == "CRS"))+
  
  geom_errorbarh(aes(xmin = D17O_pmg - D17O_err, xmax = D17O_pmg + D17O_err, y = -bottom.depth), height = 0)+
  geom_point(aes(x = D17O_pmg, y = -bottom.depth, fill = as.factor(date.collection)), size = 4, shape = 21)+
  #geom_line(aes(x = D17O_pmg, y = -bottom.depth, group = date.collection), color = MOJcolor, size = 0.5)+
  scale_fill_manual(values = c("tomato", "dodgerblue",  "chartreuse", "tomato3", "dodgerblue4", "chartreuse4", "tomato4"))+
  geom_vline(xintercept = mean(MOJmetw$D17O_pmg[MOJmetw$Water.Type == "Precipitation"], na.rm = TRUE), linetype = "longdash")+
  scale_y_continuous(limits = c(-95,0), expand = c(0, 0), breaks=seq(-90,0,10)) +
  scale_x_continuous(limits = c(-150,40), expand = c(0,0), breaks = seq(-150, 40, 20))+
  xlab(expression(Delta^"'17"*"O (pmg, SMOW-SLAP)"))+
  ylab("")+
  MStheme_isos+   theme(axis.text.y = element_blank() )

tryprof3CRS

profsummMOJ_CRS<- ggarrange(tryprofCRS, tryprof2CRS, tryprof3CRS, ncol = 3, legend = "none", align = "hv")

profsummMOJ_CRS #not sure why common legend isnt working
ggsave(filename = "MOJ_CRSprofs.pdf", plot = profsummMOJ_CRS, device = cairo_pdf, height=6,width=12, path = path.to.figs)


tryprofJT <- ggplot(data = MOJsw %>% filter (siteID.2 == "JT"))+
  #geom_path(aes(x = d18O.mean, y = -bottom.depth, group = date.collection), color = MOJcolor, size = 0.5)+
  geom_point(aes(x = d18O.mean, y = -bottom.depth, fill = as.factor(date.collection)), size = 4, shape = 21)+
  scale_fill_manual(values = c("tomato", "dodgerblue",  "chartreuse", "tomato3", "dodgerblue4", "chartreuse4", "tomato4"))+
  geom_vline(xintercept = mean.weighted.precip.MOJ.d18O, linetype = "longdash")+
  scale_y_continuous(limits = c(-95,0), expand = c(0, 0), breaks=seq(-90,0,10)) +
  scale_x_continuous(limits = c(-13,11), expand = c(0,0), breaks = seq(-12, 10, 2))+
  xlab(expression(paste(delta^{18}, "O (‰, SMOW)")))+
  #ylab("Depth (cm)")+
  ylab("")+
  MStheme_isos+
  theme(axis.text.y = element_blank() )
tryprofJT

tryprof2JT <- ggplot(data = MOJsw %>% filter (siteID.2 == "JT"))+
  geom_point(aes(x = d_excess.mean, y = -bottom.depth, fill = as.factor(date.collection)), size = 4, shape = 21)+
  scale_fill_manual(values = c("tomato", "dodgerblue",  "chartreuse", "tomato3", "dodgerblue4", "chartreuse4", "tomato4"))+
  #geom_line(aes(x = d_excess.mean, y = -bottom.depth, group = date.collection), color = MOJcolor, size = 0.5)+
  geom_vline(xintercept = mean.weighted.precip.MOJ.dexcess, linetype = "longdash")+
  scale_x_continuous(limits = c(-125,15), expand = c(0,0), breaks = seq(-120, 10, 20))+
  scale_y_continuous(limits = c(-95,0), expand = c(0, 0), breaks=seq(-90,0,10)) +
  xlab(expression(paste(italic("d"),"-excess (‰, SMOW)"))) +
  ylab("")+
  MStheme_isos+   theme(axis.text.y = element_blank() )

tryprof2JT

tryprof3JT <- ggplot(data = MOJsw %>% filter (siteID.2 == "JT"))+
  geom_errorbarh(aes(xmin = D17O_pmg - D17O_err, xmax = D17O_pmg + D17O_err, y = -bottom.depth), height = 0)+
  geom_point(aes(x = D17O_pmg, y = -bottom.depth, fill = as.factor(date.collection)), size = 4, shape = 21)+
  #geom_line(aes(x = D17O_pmg, y = -bottom.depth, group = date.collection), color = MOJcolor, size = 0.5)+
  scale_fill_manual(values = c("tomato", "dodgerblue",  "chartreuse", "tomato3", "dodgerblue4", "chartreuse4", "tomato4"))+
  scale_y_continuous(limits = c(-95,0), expand = c(0, 0), breaks=seq(-90,0,10)) +
  scale_x_continuous(limits = c(-150,40), expand = c(0,0), breaks = seq(-150, 40, 20))+
  geom_vline(xintercept = mean(MOJmetw$D17O_pmg[MOJmetw$Water.Type == "Precipitation"], na.rm = TRUE), linetype = "longdash")+
  xlab(expression(Delta^"'17"*"O (pmg, SMOW-SLAP)"))+
  ylab("")+
  MStheme_isos+   theme(axis.text.y = element_blank() )
tryprof3JT

profsummMOJ_JT<- ggarrange(tryprofJT, tryprof2JT, tryprof3JT, ncol = 3, legend = "none", align = "hv")

profsummMOJ_JT #not sure why common legend isnt working
ggsave(filename = "MOJ_JTprofs.pdf", plot = profsummMOJ_JT, device = cairo_pdf, height=6,width=12, path = path.to.figs)


tryprofPJ <- ggplot(data = MOJsw %>% filter (siteID.2 == "PJ"))+
  #geom_path(aes(x = d18O.mean, y = -bottom.depth, group = date.collection), color = MOJcolor, size = 0.5)+
  geom_point(aes(x = d18O.mean, y = -bottom.depth, fill = as.factor(date.collection)), size = 4, shape = 21)+
  scale_fill_manual(values = c( "dodgerblue",  "chartreuse", "tomato3", "dodgerblue4", "chartreuse4"))+
  
  #scale_fill_manual(values = c("deepskyblue", "darkolivegreen1", "darkorange3", "deepskyblue3", "darkolivegreen3"))+
  
  geom_vline(xintercept = mean.weighted.precip.MOJ.d18O, linetype = "longdash")+
  scale_y_continuous(limits = c(-95,0), expand = c(0, 0), breaks=seq(-90,0,10)) +
  scale_x_continuous(limits = c(-13,11), expand = c(0,0), breaks = seq(-12, 10, 2))+
  
  xlab(expression(paste(delta^{18}, "O (‰, SMOW)")))+
  #ylab("Depth (cm)")+
  ylab("")+
  MStheme_isos+
  theme(axis.text.y = element_blank() )
tryprofPJ

tryprof2PJ <- ggplot(data = MOJsw %>% filter (siteID.2 == "PJ"))+
  #geom_path(aes(x = d_excess.mean, y = -bottom.depth, group = date.collection), color = MOJcolor, size = 0.5)+
  geom_point(aes(x = d_excess.mean, y = -bottom.depth, fill = as.factor(date.collection)), size = 4, shape = 21)+
  scale_fill_manual(values = c( "dodgerblue",  "chartreuse", "tomato3", "dodgerblue4", "chartreuse4"))+
  geom_vline(xintercept = mean.weighted.precip.MOJ.dexcess, linetype = "longdash")+
  scale_x_continuous(limits = c(-125,15), expand = c(0,0), breaks = seq(-120, 10, 20))+
  scale_y_continuous(limits = c(-95,0), expand = c(0, 0), breaks=seq(-90,0,10)) +
  xlab(expression(paste(italic("d"),"-excess (‰, SMOW)"))) +
  ylab("")+
  MStheme_isos+   theme(axis.text.y = element_blank() )

tryprof2PJ

tryprof3PJ <- ggplot(data = MOJsw %>% filter (siteID.2 == "PJ"))+
 # geom_path(aes(x = D17O_pmg, y = -bottom.depth, group = date.collection), color = MOJcolor, size = 0.5)+
  geom_errorbarh(aes(xmin = D17O_pmg - D17O_err, xmax = D17O_pmg + D17O_err, y = -bottom.depth), height = 0)+
    geom_point(aes(x = D17O_pmg, y = -bottom.depth, fill = as.factor(date.collection)), size = 4, shape = 21)+
  scale_fill_manual(values = c( "dodgerblue",  "chartreuse", "tomato3", "dodgerblue4", "chartreuse4"))+
  geom_vline(xintercept = mean(MOJmetw$D17O_pmg[MOJmetw$Water.Type == "Precipitation"], na.rm = TRUE), linetype = "longdash")+
  scale_y_continuous(limits = c(-95,0), expand = c(0, 0), breaks=seq(-90,0,10)) +
  scale_x_continuous(limits = c(-150,40), expand = c(0,0), breaks = seq(-150, 40, 20))+
  xlab(expression(Delta^"'17"*"O (pmg, SMOW-SLAP)"))+
  ylab("")+
  MStheme_isos+   theme(axis.text.y = element_blank() )
tryprof3PJ

profsummMOJ_PJ<- ggarrange(tryprofPJ, tryprof2PJ, tryprof3PJ, ncol = 3, legend = "none", align = "hv")

profsummMOJ_PJ #not sure why common legend isnt working
ggsave(filename = "MOJ_PJprofs.pdf", plot = profsummMOJ_PJ, device = cairo_pdf, height=6,width=12, path = path.to.figs)

#Mojave box and whisker plots
MOJsw$depth_cat <- NA
MOJsw$depth_cat[MOJsw$bottom.depth <= 5] <- "k_0 to 5"
MOJsw$depth_cat[MOJsw$bottom.depth > 5 & MOJsw$bottom.depth <= 25] <- "j_5 to 25"
MOJsw$depth_cat[MOJsw$bottom.depth > 25 & MOJsw$bottom.depth <= 40] <- "i_25 to 40"
MOJsw$depth_cat[MOJsw$bottom.depth > 40 & MOJsw$bottom.depth <= 55] <- "h_40 to 55"
MOJsw$depth_cat[MOJsw$bottom.depth > 55 & MOJsw$bottom.depth <= 70] <- "g_55 to 70"
MOJsw$depth_cat[MOJsw$bottom.depth > 70] <- "f_70 to 100"

MOJsw$soilormet <- "soil"

MOJmetw$depth_cat[MOJmetw$Water.Type == "Precipitation"] <- "e_Precipitation"
#MOJmetw$depth_cat[MOJmetw$Water.Type == "Puddle"]  <- "a_Assorted Surface"
MOJmetw$depth_cat[MOJmetw$Water.Type == "Well"]  <- "b_Well"
#MOJmetw$depth_cat[MOJmetw$Water.Type == "Snow"]  <- "a_Assorted Surface"
MOJmetw$depth_cat[MOJmetw$Water.Type == "Spring"]  <- "a_Spring"

MOJmetw$soilormet <- "met"

MOJdataforboxes <- data.frame(d18O.mean = c(MOJsw$d18O.mean, MOJmetw$d18O.mean), d2H.mean = c(MOJsw$d2H.mean, MOJmetw$d2H.mean), 
                              d_excess.mean  = c(MOJsw$d_excess.mean, MOJmetw$d_excess.mean), 
                              D17O_pmg = c(MOJsw$D17O_pmg, MOJmetw$D17O_pmg), depth_cat = c(MOJsw$depth_cat, MOJmetw$depth_cat), 
                              soilormet =c(MOJsw$soilormet, MOJmetw$soilormet) )

counts <- MOJdataforboxes %>% count(depth_cat)
print(counts)

counts_D17O <- MOJdataforboxes %>% filter(!is.na(D17O_pmg)) %>% count(depth_cat)
print(counts_D17O)

MOJboxd18O <- ggplot()+
  geom_boxplot(data = MOJdataforboxes, aes(x = depth_cat, y = d18O.mean, fill= soilormet), outlier.shape = NA)+
  geom_jitter(data = MOJdataforboxes, aes(x = depth_cat, y = d18O.mean, fill = soilormet, shape = soilormet), color = "black")+
  scale_shape_manual(values = c(22, 21))+
  coord_flip()+
  scale_fill_manual(values = c(MOJcolor, MOJcolor.minor))+
  #scale_x_discrete(labels = c("River or Stream","Precipitation",   "40-65 cm","25-40 cm","5-25 cm", "0-5 cm"))+
  #scale_x_discrete(labels = c("","",   "","","", ""))+
  scale_y_continuous(limits = c(-13,11), expand = c(0,0), breaks = seq(-12, 10, 2))+
  xlab("")+
  ylab(expression(paste(delta^{18}, "O (‰, SMOW)")))+
  MStheme_isos + theme(axis.text.y = element_blank() )
MOJboxd18O

MOJboxdxs <- ggplot()+
  geom_boxplot(data = MOJdataforboxes, aes(x = depth_cat, y = d_excess.mean, fill= soilormet), outlier.shape = NA)+
  geom_jitter(data = MOJdataforboxes, aes(x = depth_cat, y = d_excess.mean, fill = soilormet, shape = soilormet), color = "black")+
  scale_shape_manual(values = c(22, 21))+
  scale_fill_manual(values = c(MOJcolor, MOJcolor.minor))+
  scale_x_discrete(labels = c("","",   "","","", ""))+
  #scale_x_discrete(labels = c("River or Stream","Precipitation",   "40-65 cm","25-40 cm","5-25 cm", "0-5 cm"))+
  scale_y_continuous(limits = c(-125,15), expand = c(0,0), breaks = seq(-120, 10, 20))+
  coord_flip()+
  xlab("")+
  ylab(expression(paste(italic("d"),"-excess (‰, SMOW)"))) +
  MStheme_isos + theme(axis.text.y = element_blank() )
MOJboxdxs

MOJboxD17O <- ggplot()+
  geom_boxplot(data = MOJdataforboxes, aes(x = depth_cat, y = D17O_pmg, fill= soilormet), outlier.shape = NA)+
  geom_jitter(data = MOJdataforboxes, aes(x = depth_cat, y = D17O_pmg, fill = soilormet, shape = soilormet), color = "black")+
  scale_shape_manual(values = c(22, 21))+
  scale_fill_manual(values = c(MOJcolor, MOJcolor.minor))+
  scale_x_discrete(labels = c("","",   "","","", ""))+
  scale_y_continuous(limits = c(-150,40), expand = c(0,0), breaks = seq(-150, 40, 20))+
  coord_flip()+
  ylab(expression(Delta^"'17"*"O (pmg, SMOW-SLAP)"))+
  xlab("")+
  MStheme_isos + theme(axis.text.y = element_blank() )
MOJboxD17O

MOJboxes <- ggarrange(MOJboxd18O, MOJboxdxs, MOJboxD17O, nrow =1, ncol = 3, legend = "none")
MOJboxes

ggsave(filename = "MOJboxes.pdf", plot = MOJboxes, device = cairo_pdf, height=3,width=12, path = path.to.figs)

MOJprofandboxes <- ggarrange(profsummMOJ_CRS, profsummMOJ_JT, profsummMOJ_PJ, MOJboxes, nrow = 4, ncol =1, heights  = c(1,1,1,0.75))
MOJprofandboxes

ggsave(file = "MOJprofandboxes.pdf", plot = MOJprofandboxes, width = 11.5, height = 10, device = cairo_pdf, path.to.figs)

#calculate DD17O 
mean(MOJprecip$D17O_pmg, na.rm = TRUE) - mean(MOJsw$D17O_pmg[MOJsw$bottom.depth > 40], na.rm = TRUE)
mean(MOJprecip$d18O.mean, na.rm = TRUE) - mean(MOJsw$d18O.mean[MOJsw$bottom.depth > 40], na.rm = TRUE)
mean(MOJprecip$d2H.mean, na.rm = TRUE) - mean(MOJsw$d2H.mean[MOJsw$bottom.depth > 40], na.rm = TRUE)
mean(MOJprecip$d_excess.mean, na.rm = TRUE) - mean(MOJsw$d_excess.mean[MOJsw$bottom.depth > 40], na.rm = TRUE)

mean.weighted.precip.MOJ.d18O - mean(MOJsw$d18O.mean[MOJsw$bottom.depth > 40], na.rm = TRUE)
mean.weighted.precip.MOJ.d2H  - mean(MOJsw$d2H.mean[MOJsw$bottom.depth > 40], na.rm = TRUE)
mean.weighted.precip.MOJ.dexcess - mean(MOJsw$d_excess.mean[MOJsw$bottom.depth > 40], na.rm = TRUE)


#Jornada ###########

#depth profiles for Jornada
tryprofJOR <- ggplot(data = JORsw)+
  #geom_path(aes(x = d18O.mean, y = -bottom.depth, group = date.collection), color = MOJcolor, size = 0.5)+
  geom_point(aes(x = d18O.mean, y = -bottom.depth, fill = as.factor(date.collection)), size = 4, shape = 21)+
  scale_fill_manual(values = c("tomato", "purple", "chartreuse", "tomato3", "dodgerblue4", "purple4", "tomato4"))+
  
  geom_vline(xintercept = mean.weighted.precip.JOR.d18O, linetype = "longdash")+
  scale_y_continuous(limits = c(-90,0), expand = c(0, 0), breaks=seq(-90,0,10)) +
  scale_x_continuous(limits = c(-20,15), expand = c(0,0), breaks = seq(-20, 15, 5))+

  xlab(expression(paste(delta^{18}, "O (‰, SMOW)")))+
  ylab("")+
  MStheme_isos+  theme(axis.text.y = element_blank() )

  
tryprofJOR

tryprof2JOR <- ggplot(data = JORsw)+
  #geom_path(aes(x = d_excess.mean, y = -bottom.depth, group = date.collection), color = MOJcolor, size = 0.5)+
  geom_point(aes(x = d_excess.mean, y = -bottom.depth, fill = as.factor(date.collection)), size = 4, shape = 21)+
  scale_fill_manual(values = c("tomato", "purple", "chartreuse", "tomato3", "dodgerblue4", "purple4", "tomato4"))+
  
  geom_vline(xintercept = mean.weighted.precip.JOR.dexcess, linetype = "longdash")+
  scale_y_continuous(limits = c(-90,0), expand = c(0, 0), breaks=seq(-90,0,10)) +
  scale_x_continuous(limits = c(-125,25), expand = c(0,0), breaks = seq(-120, 25, 20))+
  ylab("")+
  xlab(expression(paste(italic("d"),"-excess (‰, SMOW)"))) +
  MStheme_isos +   theme(axis.text.y = element_blank() )
tryprof2JOR

tryprof3JOR <- ggplot(data = JORsw)+
  #geom_path(aes(x = D17O_pmg, y = -bottom.depth, group = date.collection), color = MOJcolor, size = 0.5)+
  geom_point(aes(x = D17O_pmg, y = -bottom.depth, fill = as.factor(date.collection)), size = 4, shape = 21)+
  scale_fill_manual(values = c("tomato", "purple", "chartreuse", "tomato3", "dodgerblue4", "purple4", "tomato4"))+
  
  geom_vline(xintercept = mean(JORmetw$D17O_pmg[JORmetw$Water.Type == "Precipitation"], na.rm = TRUE), linetype = "longdash")+
  scale_y_continuous(limits = c(-90,0), expand = c(0, 0), breaks=seq(-90,0,10)) +
  scale_x_continuous(limits = c(-150,50), expand = c(0,0), breaks = seq(-150, 50, 20))+
  xlab(expression(Delta^"'17"*"O (pmg, SMOW-SLAP)"))+
  ylab("")+
  MStheme_isos+   theme(axis.text.y = element_blank() )
tryprof3JOR

profsummJOR<- ggarrange(tryprofJOR, tryprof2JOR, tryprof3JOR, ncol = 3, legend = "none")

profsummJOR #not sure why common legend isnt working
ggsave(filename = "JORprofs.pdf", plot = profsummJOR, device = cairo_pdf, height=6,width=12, path = path.to.figs)


#Jornada box and whisker plots
JORsw$depth_cat <- NA
JORsw$depth_cat[JORsw$bottom.depth <= 10] <- "k_0 to 10"
JORsw$depth_cat[JORsw$bottom.depth > 10 & JORsw$bottom.depth <= 25] <- "j_10 to 25"
JORsw$depth_cat[JORsw$bottom.depth > 25 & JORsw$bottom.depth <= 40] <- "i_25 to 40"
JORsw$depth_cat[JORsw$bottom.depth > 40 & JORsw$bottom.depth <= 55] <- "h_40 to 55"
JORsw$depth_cat[JORsw$bottom.depth > 55 & JORsw$bottom.depth <= 70] <- "g_55 to 70"
JORsw$depth_cat[JORsw$bottom.depth > 70] <- "f_70 to 100"

JORsw$soilormet <- "soil"

JORmetw$depth_cat[JORmetw$Water.Type == "Precipitation"] <- "b_Precipitation"
JORmetw$depth_cat[JORmetw$Water.Type == "Spring"]  <- "a_Spring"

JORmetw$soilormet <- "met"

JORdataforboxes <- data.frame(d18O.mean = c(JORsw$d18O.mean, JORmetw$d18O.mean), d2H.mean = c(JORsw$d2H.mean, JORmetw$d2H.mean), 
                              d_excess.mean  = c(JORsw$d_excess.mean, JORmetw$d_excess.mean), 
                              D17O_pmg = c(JORsw$D17O_pmg, JORmetw$D17O_pmg), depth_cat = c(JORsw$depth_cat, JORmetw$depth_cat), 
                              soilormet =c(JORsw$soilormet, JORmetw$soilormet) )

counts <- JORdataforboxes %>% count(depth_cat)
print(counts)

counts_D17O <- JORdataforboxes %>% filter(!is.na(D17O_pmg)) %>% count(depth_cat)
print(counts_D17O)

JORboxd18O <- ggplot()+
  geom_boxplot(data = JORdataforboxes, aes(x = depth_cat, y = d18O.mean, fill= soilormet), outlier.shape = NA)+
  geom_jitter(data = JORdataforboxes, aes(x = depth_cat, y = d18O.mean, fill = soilormet, shape = soilormet), color = "black")+
  scale_shape_manual(values = c(22, 21))+
  coord_flip()+
  scale_fill_manual(values = c(JORcolor, JORcolor.minor))+
  #scale_x_discrete(labels = c("River or Stream","Precipitation",   "40-65 cm","25-40 cm","5-25 cm", "0-5 cm"))+
  #scale_x_discrete(labels = c("","",   "","","", ""))+
  scale_y_continuous(limits = c(-20,15), expand = c(0,0), breaks = seq(-20, 15, 5))+
  xlab("")+
  ylab(expression(paste(delta^{18}, "O (‰, SMOW)")))+
  MStheme_isos + theme(axis.text.y = element_blank() )
JORboxd18O

JORboxdxs <- ggplot()+
  geom_boxplot(data = JORdataforboxes, aes(x = depth_cat, y = d_excess.mean, fill= soilormet), outlier.shape = NA)+
  geom_jitter(data = JORdataforboxes, aes(x = depth_cat, y = d_excess.mean, fill = soilormet, shape = soilormet), color = "black")+
  scale_shape_manual(values = c(22, 21))+
  scale_fill_manual(values = c(JORcolor, JORcolor.minor))+
  scale_x_discrete(labels = c("","",   "","","", ""))+
  #scale_x_discrete(labels = c("River or Stream","Precipitation",   "40-65 cm","25-40 cm","5-25 cm", "0-5 cm"))+
  scale_y_continuous(limits = c(-125,25), expand = c(0,0), breaks = seq(-120, 25, 20))+
  coord_flip()+
  xlab("")+
  ylab(expression(paste(italic("d"),"-excess (‰, SMOW)"))) +
  MStheme_isos +   theme(axis.text.y = element_blank() )

JORboxdxs

JORboxD17O <- ggplot()+
  geom_boxplot(data = JORdataforboxes, aes(x = depth_cat, y = D17O_pmg, fill= soilormet), outlier.shape = NA)+
  geom_jitter(data = JORdataforboxes, aes(x = depth_cat, y = D17O_pmg, fill = soilormet, shape = soilormet), color = "black")+
  scale_shape_manual(values = c(22, 21))+
  scale_fill_manual(values = c(JORcolor, JORcolor.minor))+
  scale_x_discrete(labels = c("","",   "","","", ""))+
  scale_y_continuous(limits = c(-150,50), expand = c(0,0), breaks = seq(-150, 50, 20))+
  coord_flip()+
  ylab(expression(Delta^"'17"*"O (pmg, SMOW-SLAP)"))+
  xlab("")+
  MStheme_isos+  theme(axis.text.y = element_blank() )
JORboxD17O

JORboxes <- ggarrange(JORboxd18O, JORboxdxs, JORboxD17O, nrow =1, ncol = 3, legend = "none")
JORboxes

ggsave(filename = "JORboxes.pdf", plot = JORboxes, device = cairo_pdf, height=3,width=12, path = path.to.figs)

JORprofandboxes <- ggarrange(profsummJOR, JORboxes, nrow = 2, ncol =1, heights  = c(1,0.5))
JORprofandboxes

ggsave(file = "JORprofandboxes.pdf", plot = JORprofandboxes, width = 11.5, height = 8, device = cairo_pdf, path.to.figs)


mean(JORprecip$D17O_pmg, na.rm = TRUE) - mean(JORsw$D17O_pmg[JORsw$bottom.depth > 40], na.rm = TRUE)
#mean(JORprecip$d18O.mean, na.rm = TRUE) - mean(JORsw$d18O.mean[JORsw$bottom.depth > 40], na.rm = TRUE)
#mean(JORprecip$d2H.mean, na.rm = TRUE) - mean(JORsw$d2H.mean[JORsw$bottom.depth > 40], na.rm = TRUE)
#mean(JORprecip$d_excess.mean, na.rm = TRUE) - mean(JORsw$d_excess.mean[JORsw$bottom.depth > 40], na.rm = TRUE)

mean.weighted.precip.JOR.d18O - mean(JORsw$d18O.mean[JORsw$bottom.depth > 40], na.rm = TRUE)
mean.weighted.precip.JOR.d2H - mean(JORsw$d2H.mean[JORsw$bottom.depth > 40], na.rm = TRUE)
mean.weighted.precip.JOR.dexcess - mean(JORsw$d_excess.mean[JORsw$bottom.depth > 40], na.rm = TRUE)


######### ESGR ####
tryprofESGR <- ggplot(data = ESGRsw)+
  #geom_path(aes(x = d18O.mean, y = -bottom.depth, group = date.collection), color = ESGRcolor, size = 0.5)+
    geom_point(aes(x = d18O.mean, y = -bottom.depth, fill = as.factor(date.collection)), size = 4, shape = 21)+
  scale_fill_manual(values = c("tomato", "tomato3", "purple", "chartreuse", "tomato4","red", "dodgerblue", "chartreuse3","chartreuse4","lightsalmon", "lightsalmon3", "pink"))+
  geom_vline(xintercept = mean.weighted.precip.ESGR.d18O, linetype = "longdash")+
  scale_y_continuous(limits = c(-95,0), expand = c(0, 0), breaks=seq(-90,0,10)) +
  scale_x_continuous(limits = c(-25,1), expand = c(0,0), breaks = seq(-25, 0, 5))+
  ylab("")+
  xlab(expression(paste(delta^{18}, "O (‰, SMOW)")))+
  MStheme_isos+ theme(axis.text.y = element_blank() )
tryprofESGR

tryprof2ESGR <- ggplot(data = ESGRsw)+
  #geom_path(aes(x = d_excess.mean, y = -bottom.depth, group = date.collection), color = ESGRcolor, size = 0.5)+
  geom_point(aes(x = d_excess.mean, y = -bottom.depth, fill = as.factor(date.collection)), size = 4, shape = 21)+
  scale_fill_manual(values = c("tomato", "tomato3", "purple", "chartreuse", "tomato4","red", "dodgerblue", "chartreuse3","chartreuse4","lightsalmon", "lightsalmon3", "pink"))+
  geom_vline(xintercept = mean.weighted.precip.ESGR.dexcess, linetype = "longdash")+
  scale_y_continuous(limits = c(-95,0), expand = c(0, 0), breaks=seq(-90,0,10)) +
  scale_x_continuous(limits = c(-10,27), expand = c(0,0), breaks = seq(-10, 25, 5))+
  ylab("")+
  xlab(expression(paste(italic("d"),"-excess (‰, SMOW)"))) +
  MStheme_isos + theme(axis.text.y = element_blank() )
tryprof2ESGR

tryprof3ESGR <- ggplot(data = ESGRsw)+
  #geom_path(aes(x = D17O_pmg, y = -bottom.depth, group = date.collection), color = ESGRcolor, size = 0.5)+
  geom_point(aes(x = D17O_pmg, y = -bottom.depth, fill = as.factor(date.collection)), size = 4, shape = 21)+
  scale_fill_manual(values = c("tomato", "tomato3", "purple", "chartreuse", "tomato4","red", "dodgerblue", "chartreuse3","chartreuse4","lightsalmon", "lightsalmon3", "pink"))+
  geom_vline(xintercept = mean(ESGRmetw$D17O_pmg[ESGRmetw$Water.Type == "Precipitation"], na.rm = TRUE), linetype = "longdash")+
  
  scale_y_continuous(limits = c(-95,0), expand = c(0, 0), breaks=seq(-90,0,10)) +
  scale_x_continuous(limits = c(-10,52), expand = c(0,0), breaks = seq(-10, 50, 10))+
  xlab(expression(Delta^"'17"*"O (pmg, SMOW-SLAP)"))+
  ylab("")+
  MStheme_isos + theme(axis.text.y = element_blank() )
tryprof3ESGR

profsummESGR<- ggarrange(tryprofESGR, tryprof2ESGR, tryprof3ESGR, nrow = 1, ncol = 3,legend = "none")

profsummESGR #not sure why common legend isnt working
ggsave(filename = "ESGRprofs.pdf", plot = profsummESGR, device = cairo_pdf, height=6,width=12, path = path.to.figs)

#box and whisker plots
ESGRsw$depth_cat <- NA
ESGRsw$depth_cat[ESGRsw$bottom.depth <= 10] <- "k_0 to 10"
ESGRsw$depth_cat[ESGRsw$bottom.depth > 10 & ESGRsw$bottom.depth <= 30] <- "j_10 to 30"
ESGRsw$depth_cat[ESGRsw$bottom.depth > 30 & ESGRsw$bottom.depth <= 60] <- "i_30 to 60"
ESGRsw$depth_cat[ESGRsw$bottom.depth > 60] <- "f_60 to 90"

ESGRsw$soilormet <- "soil"

ESGRmetw$depth_cat[ESGRmetw$Water.Type == "Precipitation"] <- "c_Precipitation"
ESGRmetw$depth_cat[ESGRmetw$Water.Type == "Marsh"]  <- "a_Marsh"
ESGRmetw$depth_cat[ESGRmetw$Water.Type == "River_or_stream"]  <- "b_River"
ESGRmetw$soilormet <- "met"

ESGRdataforboxes <- data.frame(d18O.mean = c(ESGRsw$d18O.mean, ESGRmetw$d18O.mean), d2H.mean = c(ESGRsw$d2H.mean, ESGRmetw$d2H.mean), 
                              d_excess.mean  = c(ESGRsw$d_excess.mean, ESGRmetw$d_excess.mean), 
                              D17O_pmg = c(ESGRsw$D17O_pmg, ESGRmetw$D17O_pmg), depth_cat = c(ESGRsw$depth_cat, ESGRmetw$depth_cat), 
                              soilormet =c(ESGRsw$soilormet, ESGRmetw$soilormet) )

counts <- ESGRdataforboxes %>% count(depth_cat)
print(counts)

n = counts_D17O <- ESGRdataforboxes %>% filter(!is.na(D17O_pmg)) %>% count(depth_cat)
print(counts_D17O)

ESGRboxd18O <- ggplot()+
  geom_boxplot(data = ESGRdataforboxes, aes(x = depth_cat, y = d18O.mean, fill= soilormet), outlier.shape = NA)+
  geom_jitter(data = ESGRdataforboxes, aes(x = depth_cat, y = d18O.mean, fill = soilormet, shape = soilormet), color = "black")+
  scale_shape_manual(values = c(22, 21))+
  coord_flip()+
  scale_fill_manual(values = c(ESGRcolor, ESGRcolor.minor))+
  #scale_x_discrete(labels = c("River or Stream","Precipitation",   "40-65 cm","25-40 cm","5-25 cm", "0-5 cm"))+
  #scale_x_discrete(labels = c("","",   "","","", ""))+
  scale_y_continuous(limits = c(-25,1), expand = c(0,0), breaks = seq(-25, 0, 5))+
  xlab("")+
  ylab(expression(paste(delta^{18}, "O (‰, SMOW)")))+
  MStheme_isos + theme(axis.text.y = element_blank() )
ESGRboxd18O 

ESGRboxdxs <- ggplot()+
  geom_boxplot(data = ESGRdataforboxes, aes(x = depth_cat, y = d_excess.mean, fill= soilormet), outlier.shape = NA)+
  geom_jitter(data = ESGRdataforboxes, aes(x = depth_cat, y = d_excess.mean, fill = soilormet, shape = soilormet), color = "black")+
  scale_shape_manual(values = c(22, 21))+
  scale_fill_manual(values = c(ESGRcolor, ESGRcolor.minor))+
  scale_x_discrete(labels = c("","",   "","","", ""))+
  #scale_x_discrete(labels = c("River or Stream","Precipitation",   "40-65 cm","25-40 cm","5-25 cm", "0-5 cm"))+
  scale_y_continuous(limits = c(-10,27), expand = c(0,0), breaks = seq(-10, 25, 5))+
  coord_flip()+
  xlab("")+
  ylab(expression(paste(italic("d"),"-excess (‰, SMOW)"))) +
  MStheme_isos + theme(axis.text.y = element_blank() )

ESGRboxdxs

ESGRboxD17O <- ggplot()+
  geom_boxplot(data = ESGRdataforboxes, aes(x = depth_cat, y = D17O_pmg, fill= soilormet), outlier.shape = NA)+
  geom_jitter(data = ESGRdataforboxes, aes(x = depth_cat, y = D17O_pmg, fill = soilormet, shape = soilormet), color = "black")+
  scale_shape_manual(values = c(22, 21))+
  scale_fill_manual(values = c(ESGRcolor, ESGRcolor.minor))+
  scale_x_discrete(labels = c("","",   "","","", ""))+
  scale_y_continuous(limits = c(-10,52), expand = c(0,0), breaks = seq(-10, 50, 10))+
  coord_flip()+
  ylab(expression(Delta^"'17"*"O (pmg, SMOW-SLAP)"))+
  xlab("")+
  MStheme_isos + theme(axis.text.y = element_blank() )
ESGRboxD17O

ESGRboxes <- ggarrange(ESGRboxd18O, ESGRboxdxs, ESGRboxD17O, nrow =1, ncol = 3, legend = "none")
ESGRboxes

ggsave(filename = "ESGRboxes.pdf", plot = ESGRboxes, device = cairo_pdf, height=3,width=12, path = path.to.figs)

ESGRprofandboxes <- ggarrange(profsummESGR, ESGRboxes, nrow = 2, ncol =1, heights  = c(1,0.5))
ESGRprofandboxes

ggsave(file = "ESGRprofandboxes.pdf", plot = ESGRprofandboxes, width = 11.5, height = 8, device = cairo_pdf, path.to.figs)


mean(ESGRprecip$D17O_pmg, na.rm = TRUE) - mean(ESGRsw$D17O_pmg[ESGRsw$bottom.depth > 40], na.rm = TRUE)
mean(ESGRprecip$d18O.mean, na.rm = TRUE) - mean(ESGRsw$d18O.mean[ESGRsw$bottom.depth > 40], na.rm = TRUE)
mean(ESGRprecip$d2H.mean, na.rm = TRUE) - mean(ESGRsw$d2H.mean[ESGRsw$bottom.depth > 40], na.rm = TRUE)
mean(ESGRprecip$d_excess.mean, na.rm = TRUE) - mean(ESGRsw$d_excess.mean[ESGRsw$bottom.depth > 40], na.rm = TRUE)

mean.weighted.precip.ESGR.d18O - mean(ESGRsw$d18O.mean[ESGRsw$bottom.depth > 40], na.rm = TRUE)
mean.weighted.precip.ESGR.d2H - mean(ESGRsw$d2H.mean[ESGRsw$bottom.depth > 40], na.rm = TRUE)
mean.weighted.precip.ESGR.dexcess - mean(ESGRsw$d_excess.mean[ESGRsw$bottom.depth > 40], na.rm = TRUE)


###############Figure 3 isotope summary figures: d18O vs dD, d18O vs d-excess, d18O vs D17O #######
#REYNOLDS three plots
REYp1 <- ggplot()+
  geom_point(data = RCEWmet.precip, aes(x = d18O_SMOW, y= dD_SMOW),fill = REYcolor, size = 2, shape = 22, alpha = 0.5, color = REYcolor)+ #data from Lohse et al 2022
  #geom_point(data = REYmetw, aes(x = d18O.mean, y= d2H.mean),fill = REYcolor, size = 3, shape = 22)+
  geom_point(data = REYmetw[REYmetw$Water.Type != "Precipitation",], aes(x = d18O.mean, y= d2H.mean),fill = REYcolor, size = 3, shape = 23)+
  geom_point(data = REYprecip, aes(x = d18O.mean, y= d2H.mean),fill = REYcolor, size = 3, shape = 22)+
  
  geom_point(data = REYsw, aes(x = d18O.mean, y= d2H.mean, fill = -bottom.depth), size = 3, shape= 21)+
  geom_point(aes(x = mean.weighted.precip.REY.d18O, y= mean.weighted.precip.REY.d2H),fill = "yellow", size = 5, shape = 22)+
  geom_abline(slope = 8, intercept = 10, color = "black", linetype = 2)+ #GMWL for comparison
  
  scale_fill_gradient(low = REYcolor, high = REYcolor.minor)+  
  #lines
  #geom_smooth(data = REYmetw,  aes(x = d18O.mean, y= d2H.mean), method = "lm", se = TRUE, color = REYcolor)+
  geom_smooth(data = REYprecip,  aes(x = d18O.mean, y= d2H.mean), method = "lm", se = TRUE, color = REYcolor)+
  geom_smooth(data = REYsw,  aes(x = d18O.mean, y= d2H.mean), method = "lm", se = TRUE, color = REYcolor.minor)+
  MStheme_isos+
  xlab(expression(delta^"18"*"O (\u2030)"))+
  ylab(expression(delta^"2"*"H (\u2030)"))+
  #RAY specific scales
  scale_x_continuous(limits = c(-21,5), expand = c(0, 0), breaks=seq(-20,5,5)) +
  scale_y_continuous(limits = c(-170,-20), expand = c(0, 0), breaks=seq(-170,-20,30)) 
  #universal scales
  #scale_x_continuous(limits = c(-25,15), expand = c(0, 0), breaks=seq(-25,15,5)) +
  #scale_y_continuous(limits = c(-175,34), expand = c(0, 0), breaks=seq(-170,30,30)) 
  
REYp1

REYp2 <- ggplot()+
  geom_point(data = RCEWmet.precip, aes(x = d18O_SMOW, y= d_excess),fill = REYcolor, size = 2, shape = 22, alpha = 0.5, color = REYcolor)+ #data from Lohse et al 2022
  geom_point(data = REYmetw[REYmetw$Water.Type != "Precipitation",], aes(x = d18O.mean, y= d_excess.mean),fill = REYcolor, size = 3, shape = 23)+
  geom_point(data = REYprecip, aes(x = d18O.mean, y= d_excess.mean),fill = REYcolor, size = 3, shape = 22)+
  geom_point(data = REYsw, aes(x = d18O.mean, y= d_excess.mean, fill = -bottom.depth), size = 3, shape= 21  )+
  geom_point(aes(x = mean.weighted.precip.REY.d18O, y= mean.weighted.precip.REY.dexcess),fill = "yellow", size = 5, shape = 22)+
  
  scale_fill_gradient(low = REYcolor, high = REYcolor.minor)+  
  
  #lines
  geom_smooth(data = REYmetw,  aes(x = d18O.mean, y= d_excess.mean), method = "lm", se = TRUE, color = REYcolor)+
  geom_smooth(data = REYsw,  aes(x = d18O.mean, y= d_excess.mean), method = "lm", se = TRUE, color = REYcolor.minor)+
  MStheme_isos+
  xlab(expression(delta^"18"*"O (\u2030)"))+
  ylab(expression(italic(d)*"-excess"*" (\u2030)"))+
  
  scale_x_continuous(limits = c(-21,5), expand = c(0, 0), breaks=seq(-20,5,5)) +
  scale_y_continuous(limits = c(-90,20), expand = c(0, 0), breaks=seq(-90,15,30)) 
  #scale_x_continuous(limits = c(-25,15), expand = c(0, 0), breaks=seq(-25,15,5)) +
  #scale_y_continuous(limits = c(-122,26), expand = c(0, 0), breaks=seq(-120,20,20)) 

REYp2

REYp3 <- ggplot()+
  
  geom_errorbar(data = REYsw, aes(x = d18O.mean, ymax = D17O_pmg + D17O_err, ymin = D17O_pmg - D17O_err), width = 0)+
  geom_errorbar(data = REYmetw, aes(x = d18O.mean, ymax = D17O_pmg + D17O_err, ymin = D17O_pmg - D17O_err), width = 0)+
  
  geom_point(data = REYmetw[REYmetw$Water.Type != "Precipitation",], aes(x = d18O.mean, y= D17O_pmg),fill = REYcolor, size = 3, shape = 23)+
  geom_point(data = REYmetw[REYmetw$Water.Type == "Precipitation",], aes(x = d18O.mean, y= D17O_pmg),fill = REYcolor, size = 3, shape = 22)+
  geom_point(data = REYsw, aes(x = d18O.mean, y= D17O_pmg, fill = -bottom.depth), size = 3, shape= 21)+
  scale_fill_gradient(low = REYcolor, high = REYcolor.minor)+  
  geom_point(aes(x = mean.weighted.precip.REY.d18O, y= mean(REYmetw$D17O_pmg, na.rm = TRUE)),fill = "yellow", size = 5, shape = 22)+
  
  #lines
  geom_smooth(data = REYmetw,  aes(x = d18O.mean, y= D17O_pmg), method = "lm", se = TRUE, color = REYcolor)+
  geom_smooth(data = REYsw,  aes(x = d18O.mean, y= D17O_pmg), method = "lm", se = TRUE, color = REYcolor.minor)+
  MStheme_isos+
  xlab(expression(delta^"18"*"O (\u2030)"))+
  ylab(expression(Delta^"'17"*"O (per meg)"))+
  scale_x_continuous(limits = c(-21,5), expand = c(0, 0), breaks=seq(-20,5,5)) +
  scale_y_continuous(limits = c(-55,53), expand = c(0, 0), breaks=seq(-50,50,20)) 
  
  #scale_x_continuous(limits = c(-25,15), expand = c(0, 0), breaks=seq(-25,15,5)) +
  #scale_y_continuous(limits = c(-155,52), expand = c(0,0), breaks = seq(-150, 50, 50))  

REYp3

#REYsumm <- ggarrange(REYp1, REYp2, REYp3, ncol = 3, nrow = 1,common.legend = TRUE, legend = "right")
REYsumm <- ggarrange(REYp1, REYp2, REYp3, ncol = 3, nrow = 1,common.legend = TRUE, legend = "none")

REYsumm

#JORNADA three plots
#in revision adding comparison to Arellano 2026 precip data from Odessa TX

Arellano <- read.csv("/Users/juliakelson/Documents/Soil Water 17O Manuscript/data/Arellano_2026_S3_forR.csv")
Odessa <- subset(Arellano, Site == "Odessa")

JORp1 <- ggplot()+
  geom_point(data = JORprecip, aes(x = d18O.mean, y= d2H.mean),fill = JORcolor, size = 3, shape = 22)+
  geom_point(data = JORmetw[JORmetw$Water.Type != "Precipitation",], aes(x = d18O.mean, y= d2H.mean),fill = JORcolor, size = 3, shape = 23)+
  geom_point(data = JORsw, aes(x = d18O.mean, y= d2H.mean, fill = -bottom.depth), size = 3, shape= 21)+
  geom_point(aes(x = mean.weighted.precip.JOR.d18O, y= mean.weighted.precip.JOR.d2H),fill = "yellow", size = 5, shape = 22)+
  
  geom_point(data = Odessa, aes(x = d18O, y = dD), fill = JORcolor, size =2, shape = 22, alpha = 0.5)+
  
  scale_fill_gradient(low = JORcolor, high = JORcolor.minor)+  
  #lines
  geom_smooth(data = JORmetw,  aes(x = d18O.mean, y= d2H.mean), method = "lm", se = TRUE, color = JORcolor)+
  geom_smooth(data = JORsw,  aes(x = d18O.mean, y= d2H.mean), method = "lm", se = TRUE, color = JORcolor.minor)+
  geom_abline(slope = 8, intercept = 10, color = "black", linetype = 2)+
  
  MStheme_isos+
  xlab(expression(delta^"18"*"O (\u2030)"))+
  ylab(expression(delta^"2"*"H (\u2030)"))+
  scale_x_continuous(limits = c(-22,15), expand = c(0, 0), breaks=seq(-25,15,5)) +
  scale_y_continuous(limits = c(-155,35), expand = c(0, 0), breaks=seq(-150,35,30)) 
  #scale_x_continuous(limits = c(-25,15), expand = c(0, 0), breaks=seq(-25,15,5)) +
  #scale_y_continuous(limits = c(-175,34), expand = c(0, 0), breaks=seq(-170,30,30)) 

JORp1

JORp2 <- ggplot()+
  geom_point(data = JORprecip, aes(x = d18O.mean, y= d_excess.mean),fill = JORcolor, size = 3, shape = 22)+
  geom_point(data = JORmetw[JORmetw$Water.Type != "Precipitation",], aes(x = d18O.mean, y= d_excess.mean),fill = JORcolor, size = 3, shape = 23)+
  geom_point(data = JORsw, aes(x = d18O.mean, y= d_excess.mean, fill = -bottom.depth), size = 3, shape= 21  )+
  geom_point(aes(x = mean.weighted.precip.JOR.d18O, y= mean.weighted.precip.JOR.dexcess),fill = "yellow", size = 5, shape = 22)+
  
  geom_point(data = Odessa, aes(x = d18O, y = d_excess), fill = JORcolor, size =2, shape = 22, alpha = 0.5)+
  
  scale_fill_gradient(low = JORcolor, high = JORcolor.minor)+  
  
  #lines
  geom_smooth(data = JORprecip,  aes(x = d18O.mean, y= d_excess.mean), method = "lm", se = TRUE, color = JORcolor)+
  geom_smooth(data = JORsw,  aes(x = d18O.mean, y= d_excess.mean), method = "lm", se = TRUE, color = JORcolor.minor)+
  MStheme_isos+
  xlab(expression(delta^"18"*"O (\u2030)"))+
  ylab(expression(italic(d)*"-excess"*" (\u2030)"))+
  scale_x_continuous(limits = c(-22,15), expand = c(0, 0), breaks=seq(-25,15,5)) +
  scale_y_continuous(limits = c(-125,32), expand = c(0, 0), breaks=seq(-120,30,30)) 
  #scale_x_continuous(limits = c(-25,15), expand = c(0, 0), breaks=seq(-25,15,5)) +
  #scale_y_continuous(limits = c(-122,26), expand = c(0, 0), breaks=seq(-120,20,20)) 

JORp2

JORp3 <- ggplot()+
  
  geom_errorbar(data = JORsw, aes(x = d18O.mean, ymax = D17O_pmg + D17O_err, ymin = D17O_pmg - D17O_err), width = 0)+
  geom_errorbar(data = JORmetw, aes(x = d18O.mean, ymax = D17O_pmg + D17O_err, ymin = D17O_pmg - D17O_err), width = 0)+
  
  geom_point(data = JORprecip, aes(x = d18O.mean, y= D17O_pmg),fill = JORcolor, size = 3, shape = 22)+
  geom_point(data = JORmetw[JORmetw$Water.Type != "Precipitation",], aes(x = d18O.mean, y= D17O_pmg),fill = JORcolor, size = 3, shape = 23)+
  geom_point(data = JORsw, aes(x = d18O.mean, y= D17O_pmg, fill = -bottom.depth), size = 3, shape= 21)+
  geom_point(aes(x = mean.weighted.precip.JOR.d18O, y= mean(JORmetw$D17O_pmg, na.rm = TRUE)),fill = "yellow", size = 5, shape = 22)+
  
  geom_point(data = Odessa, aes(x = d18O, y = D17O_permeg), fill = JORcolor, size =2, shape = 22, alpha = 0.5)+
  
  
  scale_fill_gradient(low = JORcolor, high = JORcolor.minor)+  
  #lines
  geom_smooth(data = JORprecip,  aes(x = d18O.mean, y= D17O_pmg), method = "lm", se = TRUE, color = JORcolor)+
  geom_smooth(data = JORsw,  aes(x = d18O.mean, y= D17O_pmg), method = "lm", se = TRUE, color = JORcolor.minor)+
  MStheme_isos+
  xlab(expression(delta^"18"*"O (\u2030)"))+
  ylab(expression(Delta^"'17"*"O (per meg)"))+
  scale_x_continuous(limits = c(-22,15), expand = c(0, 0), breaks=seq(-25,15,5)) +
  scale_y_continuous(limits = c(-155,52), expand = c(0, 0), breaks=seq(-150,50,30)) 
  #scale_x_continuous(limits = c(-25,15), expand = c(0, 0), breaks=seq(-25,15,5)) +
  #scale_y_continuous(limits = c(-155,52), expand = c(0,0), breaks = seq(-150, 50, 50))  
JORp3


#JORsumm <- ggarrange(JORp1, JORp2, JORp3, ncol = 3, common.legend = TRUE, legend = "right")
JORsumm <- ggarrange(JORp1, JORp2, JORp3, ncol = 3, common.legend = TRUE, legend = "none")

JORsumm

#MOJAVE three plots
MOJp1 <- ggplot()+
  geom_point(data = MOJprecip, aes(x = d18O.mean, y= d2H.mean),fill = MOJcolor, size = 3, shape = 22)+
  geom_point(data = MOJmetw[MOJmetw$Water.Type != "Precipitation",], aes(x = d18O.mean, y= d2H.mean),fill = MOJcolor, size = 3, shape = 23)+
  geom_point(data = MOJsw, aes(x = d18O.mean, y= d2H.mean, fill = -bottom.depth, shape = siteID.2), size = 3)+
  #geom_label(data = MOJsw, aes(x = d18O.mean, y= d2H.mean, label = Identifier_1 ), nudge = 2, size = 1.5)+
  
  scale_fill_gradient(low = MOJcolor, high = MOJcolor.minor)+  
  scale_shape_manual(values = c(21,24, 25))+
  geom_point(aes(x = mean.weighted.precip.MOJ.d18O, y= mean.weighted.precip.MOJ.d2H),fill = "yellow", size = 5, shape = 22)+
  
  #lines
  geom_smooth(data = MOJprecip,  aes(x = d18O.mean, y= d2H.mean), method = "lm", se = TRUE, color = MOJcolor)+
  geom_smooth(data = MOJsw,  aes(x = d18O.mean, y= d2H.mean), method = "lm", se = TRUE, color = MOJcolor.minor)+
  geom_abline(slope = 8, intercept = 10, color = "black", linetype = 2)+
  
  MStheme_isos+
  xlab(expression(delta^"18"*"O (\u2030)"))+
  ylab(expression(delta^"2"*"H (\u2030)"))+
  scale_x_continuous(limits = c(-14,12), expand = c(0, 0), breaks=seq(-10,10,5)) +
  scale_y_continuous(limits = c(-95,25), expand = c(0, 0), breaks=seq(-80,20,20)) 
  #scale_x_continuous(limits = c(-25,15), expand = c(0, 0), breaks=seq(-25,15,5)) +
  #scale_y_continuous(limits = c(-175,34), expand = c(0, 0), breaks=seq(-170,30,30)) 
MOJp1

MOJp2 <- ggplot()+
  geom_point(data = MOJprecip, aes(x = d18O.mean, y= d_excess.mean),fill = MOJcolor, size = 3, shape = 22)+
  geom_point(data = MOJmetw[MOJmetw$Water.Type != "Precipitation",], aes(x = d18O.mean, y= d_excess.mean),fill = MOJcolor, size = 3, shape = 23)+
  geom_point(data = MOJsw, aes(x = d18O.mean, y= d_excess.mean, fill = -bottom.depth, shape = siteID.2), size = 3)+
  geom_point(aes(x = mean.weighted.precip.MOJ.d18O, y= mean.weighted.precip.MOJ.dexcess),fill = "yellow", size = 5, shape = 22)+
  
  scale_fill_gradient(low = MOJcolor, high = MOJcolor.minor)+  
  scale_shape_manual(values = c(21,24, 25))+
  #lines
  geom_smooth(data = MOJprecip,  aes(x = d18O.mean, y= d_excess.mean), method = "lm", se = TRUE, color = MOJcolor)+
  geom_smooth(data = MOJsw,  aes(x = d18O.mean, y= d_excess.mean), method = "lm", se = TRUE, color = MOJcolor.minor)+
  MStheme_isos+
  xlab(expression(delta^"18"*"O (\u2030)"))+
  ylab(expression(italic(d)*"-excess"*" (\u2030)"))+
  scale_x_continuous(limits = c(-14,12), expand = c(0, 0), breaks=seq(-10,10,5)) +
  scale_y_continuous(limits = c(-122,15), expand = c(0, 0), breaks=seq(-120,15,20)) 
  
  #scale_x_continuous(limits = c(-25,15), expand = c(0, 0), breaks=seq(-25,15,5)) +
  #scale_y_continuous(limits = c(-122,26), expand = c(0, 0), breaks=seq(-120,20,20)) 
MOJp2

MOJp3 <- ggplot()+
  
  geom_errorbar(data = MOJsw, aes(x = d18O.mean, ymax = D17O_pmg + D17O_err, ymin = D17O_pmg - D17O_err), width = 0)+
  geom_errorbar(data = MOJmetw, aes(x = d18O.mean, ymax = D17O_pmg + D17O_err, ymin = D17O_pmg - D17O_err), width = 0)+
  
  geom_point(data = MOJprecip, aes(x = d18O.mean, y= D17O_pmg),fill = MOJcolor, size = 3, shape = 22)+
  geom_point(data = MOJmetw[MOJmetw$Water.Type != "Precipitation",], aes(x = d18O.mean, y= D17O_pmg),fill = MOJcolor, size = 3, shape = 23)+
  geom_point(data = MOJsw, aes(x = d18O.mean, y= D17O_pmg, fill = -bottom.depth, shape = siteID.2), size = 3)+
  scale_fill_gradient(low = MOJcolor, high = MOJcolor.minor)+ 
  scale_shape_manual(values = c(21,24, 25))+
  geom_point(aes(x = mean.weighted.precip.MOJ.d18O, y= mean(MOJmetw$D17O_pmg, na.rm = TRUE)),fill = "yellow", size = 5, shape = 22)+
  
  #lines
  geom_smooth(data = MOJprecip,  aes(x = d18O.mean, y= D17O_pmg), method = "lm", se = TRUE, color = MOJcolor)+
  geom_smooth(data = MOJsw,  aes(x = d18O.mean, y= D17O_pmg), method = "lm", se = TRUE, color = MOJcolor.minor)+
  MStheme_isos+
  xlab(expression(delta^"18"*"O (\u2030)"))+
  ylab(expression(Delta^"'17"*"O (per meg)")) +
  scale_x_continuous(limits = c(-14,12), expand = c(0, 0), breaks=seq(-10,10,5)) +
  scale_y_continuous(limits = c(-155,50), expand = c(0,0), breaks = seq(-150, 40, 30))  

  #scale_x_continuous(limits = c(-25,15), expand = c(0, 0), breaks=seq(-25,15,5)) +
  #scale_y_continuous(limits = c(-155,52), expand = c(0,0), breaks = seq(-150, 50, 50))  
MOJp3

#MOJsumm <- ggarrange(MOJp1, MOJp2, MOJp3, ncol = 3, common.legend = TRUE, legend = "right")
MOJsumm <- ggarrange(MOJp1, MOJp2, MOJp3, ncol = 3, common.legend = TRUE, legend = "none")

MOJsumm

##ESGR three plots#
ESGRp1 <- ggplot()+
  geom_point(data = ESGRprecip, aes(x = d18O.mean, y= d2H.mean),fill = ESGRcolor, size = 3, shape = 22)+
  geom_point(data = ESGRmetw[ESGRmetw$Water.Type != "Precipitation",], aes(x = d18O.mean, y= d2H.mean),fill = ESGRcolor, size = 3, shape = 23)+
  
  geom_point(data = ESGRsw, aes(x = d18O.mean, y= d2H.mean, fill = -bottom.depth), size = 3, shape= 21)+
  scale_fill_gradient(low = ESGRcolor, high = ESGRcolor.minor)+  
  #lines
  geom_point(aes(x = mean.weighted.precip.ESGR.d18O, y= mean.weighted.precip.ESGR.d2H),fill = "yellow", size = 5, shape = 22)+
  
  geom_smooth(data = ESGRprecip,  aes(x = d18O.mean, y= d2H.mean), method = "lm", se = TRUE, color = ESGRcolor)+
  geom_smooth(data = ESGRsw,  aes(x = d18O.mean, y= d2H.mean), method = "lm", se = TRUE, color = ESGRcolor.minor)+
  geom_abline(slope = 8, intercept = 10, color = "black", linetype = 2)+
  
  
  MStheme_isos+
  xlab(expression(delta^"18"*"O (\u2030)"))+
  ylab(expression(delta^"2"*"H (\u2030)"))+
  #geom_label(data = ESGRsw, aes(x = d18O.mean, y= d2H.mean, label = Identifier_1 ), nudge = 2, size = 1.5)+
  scale_x_continuous(limits = c(-25,1), expand = c(0, 0), breaks=seq(-25,0,5)) +
  scale_y_continuous(limits = c(-175,15), expand = c(0, 0), breaks=seq(-170,15,30)) 
  #scale_x_continuous(limits = c(-25,15), expand = c(0, 0), breaks=seq(-25,15,5)) +
  #scale_y_continuous(limits = c(-175,34), expand = c(0, 0), breaks=seq(-170,30,30)) 

ESGRp1


ESGRp2 <- ggplot()+
  geom_point(data = ESGRprecip, aes(x = d18O.mean, y= d_excess.mean),fill = ESGRcolor, size = 3, shape = 22)+
  geom_point(data = ESGRmetw[ESGRmetw$Water.Type != "Precipitation",], aes(x = d18O.mean, y= d_excess.mean),fill = ESGRcolor, size = 3, shape = 23)+
  
  geom_point(data = ESGRsw, aes(x = d18O.mean, y= d_excess.mean, fill = -bottom.depth), size = 3, shape= 21  )+
  scale_fill_gradient(low = ESGRcolor, high = ESGRcolor.minor)+  
  geom_point(aes(x = mean.weighted.precip.ESGR.d18O, y= mean.weighted.precip.ESGR.dexcess),fill = "yellow", size = 5, shape = 22)+
  
  #lines
  geom_smooth(data = ESGRprecip,  aes(x = d18O.mean, y= d_excess.mean), method = "lm", se = TRUE, color = ESGRcolor)+
  geom_smooth(data = ESGRsw,  aes(x = d18O.mean, y= d_excess.mean), method = "lm", se = TRUE, color = ESGRcolor.minor)+
  MStheme_isos+
  xlab(expression(delta^"18"*"O (\u2030)"))+
  ylab(expression(italic(d)*"-excess"*" (\u2030)"))+
  
  scale_x_continuous(limits = c(-25,1), expand = c(0, 0), breaks=seq(-25,0,5)) +
  scale_y_continuous(limits = c(-25,30), expand = c(0, 0), breaks=seq(-20,30,10)) 
  
  #scale_x_continuous(limits = c(-25,15), expand = c(0, 0), breaks=seq(-25,15,5)) +
  #scale_y_continuous(limits = c(-122,26), expand = c(0, 0), breaks=seq(-120,20,20)) 
ESGRp2

ESGRp3 <- ggplot()+
  
  geom_errorbar(data = ESGRsw, aes(x = d18O.mean, ymax = D17O_pmg + D17O_err, ymin = D17O_pmg - D17O_err), width = 0)+
  geom_errorbar(data = ESGRmetw, aes(x = d18O.mean, ymax = D17O_pmg + D17O_err, ymin = D17O_pmg - D17O_err), width = 0)+
  
  geom_point(data = ESGRprecip, aes(x = d18O.mean, y= D17O_pmg),fill = ESGRcolor, size = 3, shape = 22)+
  geom_point(data = ESGRmetw[ESGRmetw$Water.Type != "Precipitation",], aes(x = d18O.mean, y= D17O_pmg),fill = ESGRcolor, size = 3, shape = 23)+
  geom_point(data = ESGRsw, aes(x = d18O.mean, y= D17O_pmg, fill = -bottom.depth), size = 3, shape= 21)+
  scale_fill_gradient(low = ESGRcolor, high = ESGRcolor.minor)+
  
  geom_point(aes(x = mean.weighted.precip.ESGR.d18O, y= mean(ESGRmetw$D17O_pmg, na.rm = TRUE)),fill = "yellow", size = 5, shape = 22)+
  
  #lines
  geom_smooth(data = ESGRprecip,  aes(x = d18O.mean, y= D17O_pmg), method = "lm", se = TRUE, color = ESGRcolor)+
  geom_smooth(data = ESGRsw,  aes(x = d18O.mean, y= D17O_pmg), method = "lm", se = TRUE, color = ESGRcolor.minor)+
  MStheme_isos+
  xlab(expression(delta^"18"*"O (\u2030)"))+
  ylab(expression(Delta^"'17"*"O (per meg)")) +
  
  scale_x_continuous(limits = c(-25,1), expand = c(0, 0), breaks=seq(-25,0,5)) +
  scale_y_continuous(limits = c(-50,70), expand = c(0, 0), breaks=seq(-50,70,20)) 
  #scale_x_continuous(limits = c(-25,15), expand = c(0, 0), breaks=seq(-25,15,5)) +
  #scale_y_continuous(limits = c(-155,52), expand = c(0,0), breaks = seq(-150, 50, 50))  
ESGRp3

#ESGRsumm <- ggarrange(ESGRp1, ESGRp2, ESGRp3, ncol = 3, common.legend = TRUE, legend = "right")
ESGRsumm <- ggarrange(ESGRp1, ESGRp2, ESGRp3, ncol = 3, common.legend = TRUE, legend = "none")

multipanelsumm <- ggarrange(REYsumm, JORsumm, MOJsumm, ESGRsumm, nrow = 4)

multipanelsumm

ggsave(filename = "isos_summ_rescaled.pdf", plot = multipanelsumm, device = cairo_pdf, height=8.25,width=10.75, path = path.to.figs.refined)
#not working, not sure why not, will extract in illustrator
#iso_summ_leg <- get_legend(multipanelsumm)
#as_ggplot(iso_summ_leg)

#########ALL SITES TOGETHER for final panel###########
plines <- ggplot()+
  #points
  geom_point(data = REYmetw, aes(x = d18O.mean, y= d2H.mean),fill = REYcolor, size = 3, shape = 22)+
  geom_point(data = REYsw, aes(x = d18O.mean, y= d2H.mean),fill = REYcolor.minor, size = 3, shape= 21  )+
  
  geom_point(data = MOJmetw, aes(x = d18O.mean, y= d2H.mean), shape = 22, fill = MOJcolor, size = 3 )+
  geom_point(data = MOJsw, aes(x = d18O.mean, y= d2H.mean), fill = MOJcolor.minor, size = 3, shape = 21 )+
  
  geom_point(data = JORmetw, aes(x = d18O.mean, y= d2H.mean), fill = JORcolor, size = 3, shape = 22 )+
  geom_point(data = JORsw, aes(x = d18O.mean, y= d2H.mean), fill = JORcolor.minor, size = 3, shape = 21 )+
  
  geom_point(data = ESGRmetw, aes(x = d18O.mean, y= d2H.mean), fill = ESGRcolor, size = 3, shape = 22 )+
  geom_point(data = ESGRsw, aes(x = d18O.mean, y= d2H.mean), fill = ESGRcolor.minor, size = 3, shape = 21 )+
  
  #lines
  geom_smooth(data = REYmetw,  aes(x = d18O.mean, y= d2H.mean), method = "lm", se = TRUE, color = REYcolor)+
  geom_smooth(data = REYsw,  aes(x = d18O.mean, y= d2H.mean), method = "lm", se = TRUE, color = REYcolor.minor)+
  
  geom_smooth(data = MOJmetw,  aes(x = d18O.mean, y= d2H.mean), method = "lm", se = TRUE, color = MOJcolor)+
  geom_smooth(data = MOJsw,  aes(x = d18O.mean, y= d2H.mean), method = "lm", se = TRUE, color = MOJcolor.minor)+
  
  geom_smooth(data = JORmetw,  aes(x = d18O.mean, y= d2H.mean), method = "lm", se = TRUE, color = JORcolor)+
  geom_smooth(data = JORsw,  aes(x = d18O.mean, y= d2H.mean), method = "lm", se = TRUE, color = JORcolor.minor)+
  geom_smooth(data = ESGRmetw,  aes(x = d18O.mean, y= d2H.mean), method = "lm", se = TRUE, color = ESGRcolor)+
  geom_smooth(data = ESGRsw,  aes(x = d18O.mean, y= d2H.mean), method = "lm", se = TRUE, color = ESGRcolor.minor)+

  #add in a MWL for comparison
  geom_abline(intercept = 10, slope = 8, linetype = 2)+

  xlab(expression(delta^"18"*"O (\u2030)"))+
  ylab(expression(delta^"2"*"H (\u2030)"))+
  MStheme_isos+ 
  scale_x_continuous(limits = c(-25,15), expand = c(0, 0), breaks=seq(-25,15,5)) +
  scale_y_continuous(limits = c(-175,34), expand = c(0, 0), breaks=seq(-170,30,30)) 

plot(plines)

ggsave(filename = "plines_esgr.pdf", plot = plines, device = cairo_pdf, height=3,width=3, path = path.to.figs) 



plines2 <- ggplot()+
  #points
  geom_point(data = REYmetw, aes(x = d18O.mean, y= d_excess.mean),fill = REYcolor, size = 3, shape = 22)+
  geom_point(data = REYsw, aes(x = d18O.mean, y= d_excess.mean),fill = REYcolor.minor, size = 3, shape= 21  )+
  
  geom_point(data = MOJmetw, aes(x = d18O.mean, y= d_excess.mean), shape = 22, fill = MOJcolor, size = 3 )+
  geom_point(data = MOJsw, aes(x = d18O.mean, y= d_excess.mean), fill = MOJcolor.minor, size = 3, shape = 21 )+
  
  geom_point(data = JORmetw, aes(x = d18O.mean, y= d_excess.mean), fill = JORcolor, size = 3, shape = 22 )+
  geom_point(data = JORsw, aes(x = d18O.mean, y= d_excess.mean), fill = JORcolor.minor, size = 3, shape = 21 )+
  
  geom_point(data = ESGRmetw, aes(x = d18O.mean, y= d_excess.mean),fill = ESGRcolor, size = 3, shape = 22)+
  geom_point(data = ESGRsw, aes(x = d18O.mean, y= d_excess.mean), fill = ESGRcolor.minor, size = 3, shape= 21)+
  
  #lines
  geom_smooth(data = REYmetw,  aes(x = d18O.mean, y= d_excess.mean), method = "lm", se = TRUE, color = REYcolor)+
  geom_smooth(data = REYsw,  aes(x = d18O.mean, y= d_excess.mean), method = "lm", se = TRUE, color = REYcolor.minor)+
  
  geom_smooth(data = MOJmetw,  aes(x = d18O.mean, y= d_excess.mean), method = "lm", se = TRUE, color = MOJcolor)+
  geom_smooth(data = MOJsw,  aes(x = d18O.mean, y= d_excess.mean), method = "lm", se = TRUE, color = MOJcolor.minor)+
  
  geom_smooth(data = JORmetw,  aes(x = d18O.mean, y= d_excess.mean), method = "lm", se = TRUE, color = JORcolor)+
  geom_smooth(data = JORsw,  aes(x = d18O.mean, y= d_excess.mean), method = "lm", se = TRUE, color = JORcolor.minor)+
  
  geom_smooth(data = ESGRmetw,  aes(x = d18O.mean, y= d_excess.mean), method = "lm", se = TRUE, color = ESGRcolor)+
  geom_smooth(data = ESGRsw,  aes(x = d18O.mean, y= d_excess.mean), method = "lm", se = TRUE, color = ESGRcolor.minor)+
  
  #add in a MWL for comparison
  #geom_abline(intercept = 10, slope = 8)+
  
  xlab(expression(delta^"18"*"O (\u2030)"))+
  ylab(expression(delta^"2"*"H (\u2030)"))+
  ylab(expression(italic(d)*"-excess"*" (\u2030)"))+
  scale_x_continuous(limits = c(-25,15), expand = c(0, 0), breaks=seq(-25,15,5)) +
  scale_y_continuous(limits = c(-122,26), expand = c(0, 0), breaks=seq(-120,20,20)) +
  MStheme_isos

plot(plines2)


plines3 <- ggplot()+
  
  #lines
  geom_smooth(data = REYprecip,  aes(x = d18O.mean, y= D17O_pmg), method = "lm", se = TRUE, color = REYcolor)+
  geom_smooth(data = REYsw,  aes(x = d18O.mean, y= D17O_pmg), method = "lm", se = TRUE, color = REYcolor.minor)+
  
  geom_smooth(data = MOJprecip,  aes(x = d18O.mean, y= D17O_pmg), method = "lm", se = TRUE, color = MOJcolor)+
  geom_smooth(data = MOJsw,  aes(x = d18O.mean, y= D17O_pmg), method = "lm", se = TRUE, color = MOJcolor.minor)+
  
  geom_smooth(data = JORprecip,  aes(x = d18O.mean, y= D17O_pmg), method = "lm", se = TRUE, color = JORcolor)+
  geom_smooth(data = JORsw,  aes(x = d18O.mean, y= D17O_pmg), method = "lm", se = TRUE, color = JORcolor.minor)+
  
  geom_smooth(data = ESGRprecip, aes(x = d18O.mean, y= D17O_pmg),method= "lm", se = TRUE, color = ESGRcolor)+
  geom_smooth(data = ESGRsw, aes(x = d18O.mean, y= D17O_pmg), method= "lm", se = TRUE, color = ESGRcolor)+
  
  #points
  geom_point(data = REYprecip, aes(x = d18O.mean, y= D17O_pmg),fill = REYcolor, size = 3, shape = 22)+
  geom_point(data = REYsw, aes(x = d18O.mean, y= D17O_pmg),fill = REYcolor.minor, size = 3, shape= 21  )+
  
  geom_point(data = MOJprecip, aes(x = d18O.mean, y= D17O_pmg), shape = 22, fill = MOJcolor, size = 3 )+
  geom_point(data = MOJsw, aes(x = d18O.mean, y= D17O_pmg), fill = MOJcolor.minor, size = 3, shape = 21 )+
  
  geom_point(data = JORprecip, aes(x = d18O.mean, y= D17O_pmg), fill = JORcolor, size = 3, shape = 22 )+
  geom_point(data = JORsw, aes(x = d18O.mean, y= D17O_pmg), fill = JORcolor.minor, size = 3, shape = 21 )+
  
 geom_point(data = ESGRprecip, aes(x = d18O.mean, y= D17O_pmg),fill = ESGRcolor, size = 3, shape = 22)+
  geom_point(data = ESGRsw, aes(x = d18O.mean, y= D17O_pmg), fill = ESGRcolor.minor, size = 3, shape= 21)+
  
  #add in a MWL for comparison
  #geom_abline(intercept = 10, slope = 8)+
  
  xlab(expression(delta^"18"*"O (\u2030)"))+
  ylab(expression(Delta^"'17"*"O (per meg)")) +
  scale_x_continuous(limits = c(-25,15), expand = c(0, 0), breaks=seq(-25,15,5)) +
  scale_y_continuous(limits = c(-155,52), expand = c(0,0), breaks = seq(-150, 50, 50))  +
  MStheme_isos

plot(plines3)

ggsave(filename = "plines3_yesesgr.pdf", plot = plines3, device = cairo_pdf, height=3,width=3, path = path.to.figs) 


alltogethersumm <- ggarrange(plines, plines2, plines3, nrow = 1, ncol = 3, legend = "right")
alltogethersumm
ggsave(filename = "alltogethersumm.pdf", plot = alltogethersumm, device = cairo_pdf, height=3,width=10, path = path.to.figs) 

multipanelsumm_5panels <- ggarrange(REYsumm, JORsumm, MOJsumm, ESGRsumm,alltogethersumm, nrow = 5) 

path.to.figs.postR1 <- "~/Documents/Soil Water 17O Manuscript/code_figures/MS_figures_postR1/" 
ggsave(filename = "isos_summ_5panels.pdf", plot = multipanelsumm_5panels, width=8.4,height=10.5, device = cairo_pdf, path = path.to.figs.postR1)


