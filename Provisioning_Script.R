#Provisioning Rate Repeatability Estimates 

rm(list=ls())
# 1 - Load packages that might be needed ####

library(MCMCglmm)
library(dplyr)
library(ggplot2)
library(dotwhisker)
library(tidybayes)
library(forcats)
library(beepr)
library(rptR)
library(lme4)

# 2 -  Load Dataset and tidy ####
#data produced from query ""
# In excel extract month and year from watch date into seperate columns
Prov.rate <- read.csv("FM_ProvisioningRateData_20200106.csv", header = T)

Prov.rate$BirdID <- as.factor(Prov.rate$BirdID) 
Prov.rate$Status <- as.factor(Prov.rate$Status)
Prov.rate$Sex <- as.factor(Prov.rate$SexEstimate) # 1 = male
Prov.rate$BGID <- as.factor(Prov.rate$BreedGroupID)
Prov.rate$NWID <- as.factor(Prov.rate$NestWatchID)
Prov.rate$NestID <- as.factor(Prov.rate$NestID)
#Prov.rate$MateID <- as.factor(Prov.rate$Mate) 
Prov.rate$WatchDate <- as.Date(Prov.rate$WatchDate, "%d/%m/%y")
Prov.rate$Month<- as.factor(Prov.rate$Month)
Prov.rate$Year<- as.factor(Prov.rate$Year)
Prov.rate$HelperNo <- as.factor(Prov.rate$HelperNumber)
Prov.rate$Groupsize <- as.factor(Prov.rate$GroupSize)
Prov.rate$Age <- as.integer(Prov.rate$Age)
#age squared??? 
Prov.rate$TerritoryID <- as.factor(Prov.rate$TerritoryID) 
Prov.rate$TQ <- as.numeric(Prov.rate$TQ)
Prov.rate$TQcorrected <- as.numeric(Prov.rate$TQcorrected)
Prov.rate$ClutchSize <- as.integer(Prov.rate$ClutchSize)
Prov.rate$BroodSize <- as.integer(Prov.rate$BroodSize)
Prov.rate$NoChicks <- as.integer(Prov.rate$NoFledglings)
Prov.rate$Obs <- as.factor(Prov.rate$Observer)
Prov.rate$FPID <- as.factor(Prov.rate$FieldPeriodID)

#mean centre variables 
#Prov.rate$Age.centred <- Prov.rate$Age-mean(Prov.rate$Age)

#Provisioning rate data
Prov.rate$ProvFeeds <- as.integer(Prov.rate$ProvisioningVisit)
Prov.rate$ProvFeedsPerHour <- Prov.rate$ProvFeeds/(Prov.rate$ObsDuration/60) # Standardise Provsioning feeds to account for varying lenghts of nest watches
Prov.rate$ProvFeedsPerHour <- round(Prov.rate$ProvFeedsPerHour) #round provisioning feeds to the nearest full number
Prov.rate$ProvFeedsPerHour <- as.integer(Prov.rate$ProvFeedsPerHour)

#apply a couple of standard transformation just in case
Prov.rate$logFeedsPerHour <- log(Prov.rate$ProvFeedsPerHour)
Prov.rate$sqrtFeedsPerHour <- sqrt(Prov.rate$ProvFeedsPerHour)

#Remove NAs
Prov.rate <- Prov.rate[complete.cases(Prov.rate),]
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
ggplot(data = Prov.rate, aes(x = HelperNo, y = ProvFeedsPerHour, colour= Sex)) + geom_boxplot() + geom_point(position=position_dodge(width=0.75)) + theme_bw()

#looking at provisioning rate and group size
ggplot(data = Prov.rate, aes(x = GroupSize, y = ProvFeedsPerHour, group = GroupSize)) + geom_boxplot() + geom_point() + theme_bw()

#looking at provisioning rate and age 
ggplot(data = Prov.rate, aes(x = Age, y = ProvFeedsPerHour)) + geom_point() + theme_bw() + facet_wrap(~Sex)


# 4 - Define Priors  ####

#Inverse gamma prior = non-informative
#Predict that the posterior is not different from 0 

prior1<-list(R=list(V=1, nu=0.002), 
             G = list(G1 = list(V = 1, nu = 0.002), 
                      G2 = list(V = 1, nu = 0.002)))

#Close to inverse wishart (informative)

prior2<-list(R=list(V=1, nu=0.2), 
             G = list(G1 = list(V = 1, nu = 0.2), 
                      G2 = list(V = 1, nu = 0.2)))

#Expanded prior
prior3<- list(R = list(V = 1, nu=0.002), G = list(G1 = list(V = 1,nu= 1,alpha.mu=0,alpha.V=1000), 
                                                  G2 = list(V = 1,nu= 1,alpha.mu=0,alpha.V=1000)))

#Half Lotte, half cauchy
prior4<- list(R = list(V = 1, nu=0.002), G = list(G1 = list(V = 1,nu= 0.002,alpha.mu=0,alpha.V=1000), 
                                                  G2 = list(V=1, nu=1, alpha.mu = 0, alpha.V = 25^2)))

# 5 - Running First Models ####

model1 <-MCMCglmm(ProvFeedsPerHour ~ HelperNo + Sex + GroupSize + BroodSize + TQ + Age,
                 random=~BirdID + NWID, nitt=200300, burnin=3000, thin=100, prior = prior2,
                 verbose=TRUE, 
                 family="poisson", 
                 data=Prov.rate)
beep(sound=4)


# 6 - Check diagnostics ####

summary(model1)

plot(model1$Sol)#Fixed effects (variance components)
plot(model1$VCV)#Random effects (variance components)


#Check for convergences
#Passed? (all should pass)
heidel.diag(model1$Sol)
heidel.diag(model1$VCV)


#<2 z-score, so non-significant is good (tests similarity between first and last 10% of iterations)
geweke.diag(model1$Sol) # brood size is greater than 2??
geweke.diag(model1$VCV)


#Test the independence of the fixed effects/random effects
#if autocorrelation is high you might need more thinning??
autocorr(model1$Sol)
autocorr(model1$VCV)

#Should be at least 1000, with 10,000 as goal
#maximum amount of iterations is 2,000 
#Would like to be as close to that as possible
round(sort(effectiveSize(model1$Sol)))
round(sort(effectiveSize(model1$VCV)))


# 7 - Repeatability analysis ####


#Repeatability
#Between individual variance Divided by total phenotypic variance (between individual variance + within individual variance).
rep.Prov.rate <- (model1$VCV[, "BirdID"])/
  (model1$VCV[, "BirdID"]+
     model1$VCV[,"NWID"] +
     model1$VCV[,"units"]+log(1/exp(model1$Sol[,"(Intercept)"])+1))

#posterior.mode gives the repeatability score for provisioning rate
posterior.mode(rep.Prov.rate) 
# was 0.0918 first time ran with NWID:NestID
#0.0871 with just NWID as other random effect
HPDinterval(rep.Prov.rate) 
# 0.046-0.147 for NWID:NestID
# 0.0519-0.153 for NWID as only other random


# Non_bayesian Models ####

#Straight forward standard repeatability estimate using rptR and controlling for nest watch duration
#model does not converge
# gives repeatability score of 0.078, CI 0-0.15
rpt(ProvFeedsPerHour ~ (1|BirdID), grname = "BirdID", data = Prov.rate, datatype = "Poisson", nboot = 3000, npermut = 0)

hist(Prov.rate$LogAge)
hist(Prov.rate$SqrtAge)
Prov.rate$LogAge <- log(Prov.rate$Age)
Prov.rate$SqrtAge <- sqrt(Prov.rate$Age)

hist(Prov.rate$LogTQ)
hist(Prov.rate$SqrtTQ)
Prov.rate$LogTQ <- log(Prov.rate$TQ)
Prov.rate$SqrtTQ <- sqrt(Prov.rate$TQ)

model2 <- glmer(ProvFeedsPerHour ~ HelperNo + GroupSize + Sex + SqrtAge + BroodSize + SqrtTQ + (1|BirdID) + (1|NWID), family = poisson, data = Prov.rate)
summary(model2)

