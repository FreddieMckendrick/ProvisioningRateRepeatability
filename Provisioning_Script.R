#Provisioning Rate Repeatability Estimates 

# 1 - Load packages that might be needed ####

library(MCMCglmm)
library(dplyr)
library(ggplot2)
library(dotwhisker)
library(tidybayes)
library(forcats)
library(rptR)
library(readr)

# 2 -  Load Dataset and tidy ####
#data produced from query ""
Prov.rate <- read.csv("ProvisioningRateRpt.csv", header =T)

Prov.rate$BirdID <- as.factor(Prov.rate$BirdID) 
Prov.rate$Status <- as.factor(Prov.rate$Status)
Prov.rate$Sex <- as.factor(Prov.rate$Sex)
Prov.rate$BGID <- as.factor(Prov.rate$BreedGroupID)
Prov.rate$NWID <- as.factor(Prov.rate$NestWatchID)
Prov.rate$MateID <- as.factor(Prov.rate$Mate) 
Prov.rate$NestStatus <- as.factor(Prov.rate$NestStatus)
Prov.rate$WatchDate <- as.Date(Prov.rate$WatchDate, "%d/%m/%y")
Prov.rate$Month<- as.factor(Prov.rate$Month)
Prov.rate$Year<- as.factor(Prov.rate$Year)
Prov.rate$HelperNo <- as.factor(Prov.rate$HelperNumber)
Prov.rate$Groupsize <- as.factor(Prov.rate$GroupSize)
Prov.rate$Age <- as.integer(Prov.rate$Age)
#age squared??? 
Prov.rate$TerritoryID <- as.factor(Prov.rate$TerritoryID) 
Prov.rate$TQ <- as.numeric(Prov.rate$TerritoryQuality)
Prov.rate$ClutchSize <- as.integer(Prov.rate$ClutchSize)
Prov.rate$BroodSize <- as.integer(Prov.rate$BroodSize)
Prov.rate$NoChicks <- as.integer(Prov.rate$NumberofChicks)
Prov.rate$Obs <- as.factor(Prov.rate$Observer)
Prov.rate$FPID <- as.factor(Prov.rate$FieldPeriodID)


#Provisioning rate data
Prov.rate$ProvFeeds <- as.integer(Prov.rate$ProvisioningFeeds)
Prov.rate$ProvFeedsPerHour <- Prov.rate$ProvFeeds/(Prov.rate$ObsDuration/60)
#apply a couple of standard transformation just in case
Prov.rate$logFeedsPerHour <- log(Prov.rate$ProvFeedsHour)
Prov.rate$sqrtFeedsPerHour <- sqrt(Prov.rate$ProvFeedsHour)

#check layout and structure of data
#if variables are not right go back and re run code 
head(Prov.rate)
str(Prov.rate) 

#check variability in provisioning rate 
hist(Prov.rate$ProvFeedsPerHour)
hist(Prov.rate$logFeedsPerHour)
hist(Prov.rate$sqrtFeedsPerHour)

# 3  - Visualise data ####

#looking at provisioning rate and the number of helpers 
ggplot(data = Prov.rate, aes(x = HelperNo, y = ProvFeedsPerHour, colour= Sex)) + geom_boxplot() + geom_point() + theme_bw()

#looking at provisioning rate and group size
ggplot(data = Prov.rate, aes(x = GroupSize, y = ProvFeedsPerHour, colour= Sex)) + geom_boxplot() + geom_point() + theme_bw()

#looking at provisioning rate and age 
ggplot(data = Prov.rate, aes(x = Age, y = ProvFeedsPerHour)) + geom_point() + theme_bw() + facet_wrap(~Sex)


# 4 - Define Priors  ####

#Inverse gamma prior = non-informative
#Predict that the posterior is not different from 0 

prior1<-list(R=list(V=1, nu=0.002), 
             G = list(G1 = list(V = 1, nu = 0.002), 
                      G2 = list(V = 1, nu = 0.002)))

#Close to inverse wishart (informative)

prior1<-list(R=list(V=1, nu=0.2), 
             G = list(G1 = list(V = 1, nu = 0.2), 
                      G2 = list(V = 1, nu = 0.2)))

#Expanded prior
prior3<- list(R = list(V = 1, nu=0.002), G = list(G1 = list(V = 1,nu= 1,alpha.mu=0,alpha.V=1000), 
                                                  G2 = list(V = 1,nu= 1,alpha.mu=0,alpha.V=1000)))

#Half Lotte, half cauchy
prior4<- list(R = list(V = 1, nu=0.002), G = list(G1 = list(V = 1,nu= 0.002,alpha.mu=0,alpha.V=1000), 
                                                  G2 = list(V=1, nu=1, alpha.mu = 0, alpha.V = 25^2)))