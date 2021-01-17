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
Prov.rate <- read.csv("FM_ProvisioningRate_20210114.csv", header = T)

Prov.rate$BirdID <- as.factor(Prov.rate$BirdID) 
Prov.rate$Sex <- as.factor(Prov.rate$SexEstimate) # 1 = male
Prov.rate$BGID <- as.factor(Prov.rate$BreedGroupID)
Prov.rate$NWID <- as.factor(Prov.rate$NestWatchID)
Prov.rate$NestID <- as.factor(Prov.rate$NestID)
Prov.rate$WatchDate <- as.Date(Prov.rate$WatchDate, "%d/%m/%y")
Prov.rate$Month<- as.factor(Prov.rate$Month)
Prov.rate$Year<- as.factor(Prov.rate$Year)
Prov.rate$HelperNo <- as.factor(Prov.rate$HelperNumber) # as integer or factor???
Prov.rate$Groupsize <- as.factor(Prov.rate$GroupSize)
Prov.rate$TerritoryID <- as.factor(Prov.rate$TerritoryID) 
Prov.rate$NoChicks <- as.integer(Prov.rate$NoFledglings)
Prov.rate$FPID <- as.factor(Prov.rate$FieldPeriodID)

#mean centre Age
Prov.rate$Age.centred <- Prov.rate$Age-mean(Prov.rate$Age)

#Provisioning rate data
Prov.rate$ProvFeeds <- as.integer(Prov.rate$ProvisioningVisit)
# Standardise Provsioning feeds to account for varying lenghts of nest watches
Prov.rate$ProvFeedsPerHour <- Prov.rate$ProvFeeds/(Prov.rate$ObsDuration/60) 
#round provisioning feeds to the nearest integer
Prov.rate$ProvFeedsPerHour <- round(Prov.rate$ProvFeedsPerHour) 
Prov.rate$ProvFeedsPerHour <- as.integer(Prov.rate$ProvFeedsPerHour)


#Remove NAs
Prov.rate <- Prov.rate[complete.cases(Prov.rate),]

#check layout and structure of data
head(Prov.rate)
str(Prov.rate) 

#check variability in provisioning rate 
#Does not seem to be any outliers
hist(Prov.rate$ProvFeedsPerHour)

# 3  - Visualise data ####

#This isn't really needed, just wanted to look at the relationship between Provisioning rate and fixed effects 

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


#Expanded prior from HE/Tara
#prior2<- list(R = list(V = 1, nu=0.002), G = list(G1 = list(V = 1,nu= 1,alpha.mu=0,alpha.V=1000), G2 = list(V = 1,nu= 1,alpha.mu=0,alpha.V=1000)))


# 5 - Running First Models ####

model1 <-MCMCglmm(ProvFeedsPerHour ~ HelperNo + Sex + GroupSize + BroodSize + Age,
                 random=~BirdID + NestID:NWID, nitt=200300, burnin=3000, thin=100,
                 verbose=TRUE, 
                 family="poisson", 
                 data=Prov.rate)
beep(sound=4)


# 6 - Check diagnostics ####

summary(model1)
#HelperNo, Brood Size and Age all having an effect on provisioning rate 

plot(model1$Sol)#Fixed effects (variance components) 
plot(model1$VCV)#Random effects (variance components)
#All Traces seem like that the model has worked well 
#When ran with data from just 2010-2020 the model did not converge well with BirdID
#When using data from 1995 model converges

#Check for convergences
#HelperNo3 and GroupSize failed on Halfwidth
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
     model1$VCV[,"NestID:NWID"] +
     model1$VCV[,"units"]+log(1/exp(model1$Sol[,"(Intercept)"])+1))

#posterior.mode gives the repeatability score for provisioning rate
posterior.mode(rep.Prov.rate) 
#When I ran it with set prior I got an estimate of 0.119

HPDinterval(rep.Prov.rate) 
# HPD Interval of 0.0685-0.184


# Look at males and Females individually to see if repeatability varies between sexes ####

#Subset data into two new data frames with males and females 
Male.PR <- subset(Prov.rate, Sex == 1)
Female.PR <- subset(Prov.rate, Sex ==0)

#Remove any NAs
Male.PR <- Male.PR[complete.cases(Male.PR),]
Female.PR <- Female.PR[complete.cases(Female.PR),]

#Run same model as before 
#Males
model2 <- MCMCglmm(ProvFeedsPerHour ~ HelperNo + GroupSize + BroodSize + Age,
                     random=~BirdID + NestID:NWID, nitt=200300, burnin=3000, thin=100,
                     verbose=TRUE, 
                     family="poisson", 
                     data=Male.PR)
#Females
model3 <- MCMCglmm(ProvFeedsPerHour ~ HelperNo + GroupSize + BroodSize + Age,
                   random=~BirdID + NestID:NWID, nitt=200300, burnin=3000, thin=100,
                   verbose=TRUE, 
                   family="poisson", 
                   data=Female.PR)



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

glmm1 <- glmer(ProvFeedsPerHour ~(1|BirdID) + (1|NWID), family = poisson, data = Prov.rate)
summary(model2)

