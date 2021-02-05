#Provisioning Rate Repeatability Estimates 

rm(list=ls())
# 1 - Load packages that might be needed ####

library(tidyverse)
library(Matrix)
library(coda)
library(ape)
library(MCMCglmm)
library(dotwhisker)
library(tidybayes)
library(beepr)
library(rptR)
library(lme4)
library(insight)
library(parallel)
library(readxl)
library(lubridate)

# 2 -  Load Dataset and tidy ####
#data produced from query ""
# In excel extract month and year from watch date into separate columns
Prov.rate <- read_csv("Data/FM_ProvisioningRate_20210114.csv")

Prov.rate$BirdID <- as.factor(Prov.rate$BirdID) 
Prov.rate$Status <- as.factor(Prov.rate$Status)
Prov.rate$Sex <- as.factor(Prov.rate$SexEstimate) # 1 = male
Prov.rate$BGID <- as.factor(Prov.rate$BreedGroupID)
Prov.rate$NWID <- as.factor(Prov.rate$NestWatchID)
Prov.rate$NestID <- as.factor(Prov.rate$NestID)
Prov.rate$WatchDate <- dmy(Prov.rate$WatchDate)
Prov.rate$Month<- as.factor(Prov.rate$Month)
Prov.rate$Year<- as.factor(Prov.rate$Year)
Prov.rate$HelperNo <- as.integer(Prov.rate$HelperNumber) # as integer or factor???
Prov.rate$GroupSize <- as.integer(Prov.rate$GroupSize)
Prov.rate$TerritoryID <- as.factor(Prov.rate$TerritoryID) 
Prov.rate$NoChicks <- as.integer(Prov.rate$NoFledglings)
Prov.rate$FPID <- as.factor(Prov.rate$FieldPeriodID)

#mean center Age
Prov.rate$Age_centred <- Prov.rate$Age-mean(Prov.rate$Age)
#in HE 2017 Age was also divided by 2 standard deviations
Age_sd <- sd(Prov.rate$Age)
Age_Two_sd <- 2*Age_sd
Prov.rate$Age_centred <- Prov.rate$Age_centred / Age_Two_sd

#Provisioning rate data
Prov.rate$ProvFeeds <- as.integer(Prov.rate$ProvisioningVisit)
#Standardize Provisioning feeds to account for varying lengths of nest watches
Prov.rate$ProvFeedsPerHour <- Prov.rate$ProvFeeds/(Prov.rate$ObsDuration/60) 
#round provisioning feeds to the nearest integer
Prov.rate$ProvFeedsPerHour <- round(Prov.rate$ProvFeedsPerHour) 
Prov.rate$ProvFeedsPerHour <- as.integer(Prov.rate$ProvFeedsPerHour)

#Remove NAs
#35 NAs (will identify them in Excel I believe they are in fledgling column)
NA_Values <- Prov.rate %>%  filter_all(any_vars(is.na(.)))
Prov.rate %>% 
  select_if(function(x) any(is.na(x))) %>% 
  summarise_each(funs(sum(is.na(.)))) -> extra_NA

#Remove 5 records
#Records that have no Brood size are NWID 2480,2481 from Nest 4778
#Record that has no observer = NWID 659
Prov.rate <- Prov.rate %>% 
  filter( NestID != 4778) %>% 
  filter(NWID != 659)

#check layout and structure of data
head(Prov.rate)
str(Prov.rate) 

#check variability in provisioning rate 
#Does not seem to be any outliers
hist(Prov.rate$ProvFeedsPerHour)

#Add NSA individuals to group size --------

#Add data of individuals with only one breed group
#change Field period and BirdID to factors
NSA_OneGroup <- read_csv("Data/FM_NSA_OneBG.csv")
NSA_OneGroup
NSA_OneGroup$FieldPeriodID <- as.factor(NSA_OneGroup$FieldPeriodID)
NSA_OneGroup <- NSA_OneGroup %>% 
  dplyr::rename(FPID = FieldPeriodID)
NSA_OneGroup$BirdID <- as.factor(NSA_OneGroup$BirdID)

#Add data of breed groups of individual birds
NSA_BGID <- read_csv("Data/FM_NSA.csv")
NSA_BGID
NSA_BGID$FieldPeriodID <- as.factor(NSA_BGID$FieldPeriodID)
NSA_BGID <- NSA_BGID %>% 
  dplyr::rename(FPID = FieldPeriodID)
NSA_BGID$BirdID <- as.factor(NSA_BGID$BirdID)
NSA_BGID <- NSA_BGID %>% 
  dplyr::select(-(NSA_Birds))

#Join tables together so that you have the breed group and territory numbers for birds with one breed group
NSA_OneGroup <- NSA_OneGroup %>% left_join(NSA_BGID)
NSA_OneGroup <- NSA_OneGroup %>% 
  dplyr::rename(BGID = BreedGroupID)
NSA_OneGroup$BGID <- as.factor(NSA_OneGroup$BGID)

# Identify if there is breed group with more than one NSA in the same year 
#eight groups with 2 and two groups with three NSA
NSA_OneGroup %>%
  dplyr::count(FPID, BGID, TerritoryNumber) %>% 
  dplyr::filter(n>1)

#Count the number of NSA individuals per breed group per Field Period
NumberOfBirds <- NSA_OneGroup %>%
  dplyr::count(FPID, BGID, TerritoryNumber)

#Set FPID and BGID as factor so match Prov.rate data 
NumberOfBirds$BGID <- as.factor(NumberOfBirds$BGID) 
NumberOfBirds$FPID <- as.factor(NumberOfBirds$FPID)

#Join Number of NSA individuals to Prov.rate
#Joined by matching FPID, BGID and Territory Number
Prov.rate <- Prov.rate %>% left_join(NumberOfBirds)

#Change name to recognisable variable
Prov.rate <- Prov.rate %>% 
  dplyr::rename(NSA.Individuals = n)

#Change Na to 0
Prov.rate <- Prov.rate %>% 
  dplyr::mutate_at(dplyr::vars(NSA.Individuals), ~tidyr::replace_na(.,0))

#Add NSA individuals into existing group size column
Prov.rate <- Prov.rate %>% 
  dplyr::mutate(GroupSize = GroupSize + NSA.Individuals)

# Add Territory Quality  ----------------------------------

#Load in datasets
TerritoryQuality_calc <- read_excel("Data/sys_TerritoryQuality.xlsx")
FPS_SummerIndex <- read_excel("Data/AllFPS_SummerIndex.xlsx")

# Territory quality: fixing missing datapoints 

#Is it a summer or winter season? (Calculate for all FPs whether offspring were born or not), then calculate midpoint (used to determine chronology)
FPS_SummerIndex$Season <- NA
FPS_SummerIndex$Season[FPS_SummerIndex$SummerIndex>=0.5] <- "Summer"
FPS_SummerIndex$Season[FPS_SummerIndex$SummerIndex<0.5] <- "Winter"

FPS_SummerIndex$midpoint <- as.Date((FPS_SummerIndex$PeriodStart) + 
                                      ((difftime(FPS_SummerIndex$PeriodEnd, FPS_SummerIndex$PeriodStart , units = c("days")))/2))

FPS_SummerIndex <- subset(FPS_SummerIndex, select = c("TerritoryID", "FieldPeriodID", "Season", "midpoint"))

#Bind information about TQ for all FP where data is available
TerritoryQuality_calc <- unique(subset(TerritoryQuality_calc, select = c("TerritoryID", "FieldPeriodID", "TQcorrected"), Island=="CN"))

#Link the dataset we have territory quality info for to the complete list of FPs and territories (even if there's no TQ data)
TerritoryQuality_calc <- unique(merge(TerritoryQuality_calc, FPS_SummerIndex, by.x=c("TerritoryID", "FieldPeriodID"),
                                      by.y=c("TerritoryID", "FieldPeriodID"), 
                                      all.x = TRUE,
                                      all.y = TRUE))

#Arrange 
TerritoryQuality_calc <- TerritoryQuality_calc %>%
  arrange(TerritoryID, midpoint) %>%
  group_by(TerritoryID) 

colnames(TerritoryQuality_calc)[colnames(TerritoryQuality_calc) == 'TQcorrected'] <- 'TQ'
TerritoryQuality_calc$prevTQ <- TerritoryQuality_calc$TQ
TerritoryQuality_calc$nextTQ <- TerritoryQuality_calc$TQ

#Break down into summer and winter FPs as we will need to calculate the estimated TQ based on matching seasons 
TQ_winter <- subset(TerritoryQuality_calc, Season=="Winter")
TQ_summer <- subset(TerritoryQuality_calc, Season=="Summer")

#####Sara's code 
#If there is an NA for a particular FP, work out the mean between prevTQ and nextTQ, and then replace any NAs in the TQ column with the meanTQ

#Print the value of the most recent field period before current NA
TQ_winter <- TQ_winter %>% 
  group_by(TerritoryID) %>% 
  tidyr::fill(prevTQ, .direction='updown')

#Print the value of the most recent field period after current NA
TQ_winter <- TQ_winter %>% 
  group_by(TerritoryID) %>% 
  tidyr::fill(nextTQ, .direction='downup')

#Calculate mean of prev and next TQs 
TQ_winter$TQ <- (TQ_winter$prevTQ + TQ_winter$nextTQ)/2

#Repeat for summer seasons

#Print the value of the most recent field period before current NA
TQ_summer <- TQ_summer %>% 
  group_by(TerritoryID) %>% 
  tidyr::fill(prevTQ, .direction='updown')

#Print the value of the most recent field period after current NA
TQ_summer <- TQ_summer %>% 
  group_by(TerritoryID) %>% 
  tidyr::fill(nextTQ, .direction='downup')

#Calculate mean of prev and next TQs 
TQ_summer$TQ <- (TQ_summer$prevTQ + TQ_summer$nextTQ)/2

#Bind summer and winter seasons back together
TQ <- rbind(TQ_summer, TQ_winter)

#Refine TQ dataset 
TQ <- subset(TQ, select = c("TerritoryID", "FieldPeriodID", "TQ"))


#Rename variables and set as factors so matches Provisioning data
TQ <- TQ %>% 
  dplyr::rename(FPID = FieldPeriodID)
TQ$FPID <- as.factor(TQ$FPID)
TQ$TerritoryID <- as.factor(TQ$TerritoryID)

#Join TQ to provisioning data
Prov.rate <- Prov.rate %>% 
  left_join(TQ)

#mean centre Territory Quality
Prov.rate$TQ_centred <- Prov.rate$TQ-mean(Prov.rate$TQ)
#in HE 2017 Age was also divided by 2 standard deviations
TQ_sd <- sd(Prov.rate$TQ)
TQ_Two_sd <- 2*TQ_sd
Prov.rate$TQ_centred <- Prov.rate$TQ_centred / TQ_Two_sd

#Set Prov.rate as a dataframe now it has finished being tidied 
Prov.rate <- as.data.frame(Prov.rate)


# 3  - Visualise data ####

#This isn't really needed, just wanted to look at the relationship between Provisioning rate and fixed effects 

#looking at provisioning rate and the number of helpers 
ggplot(data = Prov.rate, aes(x = HelperNo, y = ProvFeedsPerHour, colour= Sex)) + geom_boxplot() + geom_point(position=position_dodge(width=0.75)) + theme_bw()

#looking at provisioning rate and group size
ggplot(data = Prov.rate, aes(x = GroupSize, y = ProvFeedsPerHour, group = GroupSize)) + geom_boxplot() + geom_point() + theme_bw()

#looking at provisioning rate and age 
ggplot(data = Prov.rate, aes(x = Age, y = ProvFeedsPerHour)) + geom_point() + theme_bw() + facet_wrap(~Sex)


# 4 - Define Priors  ####

#Inverse-gamma prior
#non-informative prior that should be same as default
#Predicts that the posterior distribution is not different from 0 
#Three levels of G for three random effects? 
prior1<-list(R=list(V=1, nu=0.002), 
             G = list(G1 = list(V = 1, nu = 0.002), 
                      G2 = list(V = 1, nu = 0.002),
                      G3 = list(V = 1, nu = 0.002),
                      G4 = list(V = 1, nu = 0.002)))

#parameter expanded prior 
#covariance matrix set as 1000 
prior2 <- list(R = list(V = 1, nu = 0.002),
               G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                        G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                        G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                        G4 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

#same as prior2 but with lower nu
prior3 <- list(R = list(V = 1, nu = 0.002),
               G = list(G1 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1000),
                        G2 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1000),
                        G3 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1000),
                        G4 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1000)))

# 5 - Running First Models ####

#Run with default prior 
Prov_rate_MCMC_model1 <-MCMCglmm(ProvFeedsPerHour ~ HelperNo + Sex + GroupSize + BroodSize + Age_centred + TQ_centred, random = ~BirdID + NestID + NestID:NWID + Observer, nitt=650000, burnin=50000, thin=300, prior = prior2, verbose = TRUE, family = "poisson",  data = Prov.rate)
beep(sound=4)

#run with prior3
#doesn't really alter variances much at all
Prov_rate_MCMC_model2 <-MCMCglmm(ProvFeedsPerHour ~ HelperNo + Sex + GroupSize + BroodSize + Age_centred + TQ_centred, random = ~BirdID + NestID + NestID:NWID + Observer, nitt=650000, burnin=50000, thin=300, prior = prior3, verbose = TRUE, family = "poisson",  data = Prov.rate)
beep(sound=4)


#Produced different post.means and variances so it looks like NWID nested within side is effecting the results 

# 6 - Check diagnostics ####

summary(Prov_rate_MCMC_model1)
#HelperNo, Brood Size and Age all having an effect on provisioning rate 

plot(Prov_rate_MCMC_model1$Sol) #Fixed effects (variance components) 
plot(Prov_rate_MCMC_model1$VCV) #Random effects (variance components)
#All Traces seem like that the model has worked well 
#When ran with data from just 2010-2020 the model did not converge well with BirdID
#When using data from 1995 model converges

#Heidel/Welch convergence diagnostic 
#based on relative accuracy for the estimate of the mean
#If either fails then the model run should be extended 
heidel.diag(Prov_rate_MCMC_model1$Sol) # All passed
heidel.diag(Prov_rate_MCMC_model1$VCV) # All passed

#Geweke convergence disgnostic
#<2 z-score, so non-significant is good (tests similarity between first 10% and last 50% of iterations)
geweke.diag(Prov_rate_MCMC_model1$Sol) #All less than 2
geweke.diag(Prov_rate_MCMC_model1$VCV) #NestID:NWID greater than 2

#Test the independence of the fixed effects/random effects
#if autocorrelation is high then more thinning is needed before calculation of summary statistics 
#none above 0.1 so no more thinning is needed
autocorr(Prov_rate_MCMC_model1$Sol)
autocorr(Prov_rate_MCMC_model1$VCV)

#Gives the effective sampled size for each column
#maximum amount of iterations is 2,000 
#Would like to be as close to that as possible
round(sort(effectiveSize(Prov_rate_MCMC_model1$Sol)))
round(sort(effectiveSize(Prov_rate_MCMC_model1$VCV)))
#majority are 2000 with all at least 1000

# 7 - Repeatability analysis ####

#Repeatability
#Between individual variance Divided by total phenotypic variance (between individual variance + within individual variance).
rep.Prov.rate <- (Prov_rate_MCMC_model1$VCV[, "BirdID"])/ (Prov_rate_MCMC_model1$VCV[, "BirdID"]+ Prov_rate_MCMC_model1$VCV[,"NestID"] + Prov_rate_MCMC_model1$VCV[,"NestID:NWID"] + Prov_rate_MCMC_model1$VCV[,"Observer"] + Prov_rate_MCMC_model1$VCV[,"units"] + log(1/exp(Prov_rate_MCMC_model1$Sol[,"(Intercept)"])+1))

#posterior.mode gives the repeatability score for provisioning rate
posterior.mode(rep.Prov.rate) 
#When I ran it with prior 2 I got an estimate of 0.09032

HPDinterval(rep.Prov.rate) 
# HPD Interval of 0.04779-0.14101


# 8 - exploring data ------------------------------------

percentagehelpers <- subset(Prov.rate, HelperNumber >=1)
# 725 records where helpers were present 
Helper.breed.groups <- length(unique(percentagehelpers$BGID))
#256 breed groups have a helper out of 510
Percentage.BG.Helpers <- (256/510)*100
#50.39% of groups have a helper
table(percentagehelpers$HelperNo)

OneHelper <- subset(percentagehelpers, HelperNumber == 1)
OneHelperBG <- length(unique(OneHelper$BGID)) # 208/256
TwoHelper <- subset(percentagehelpers, HelperNumber == 2)
TwoHelperBG <- length(unique(TwoHelper$BGID)) #45/256
ThreeHelper <- subset(percentagehelpers, HelperNumber == 3)
ThreeHelperBG <- length(unique(ThreeHelper$BGID)) #3/256

#how many individuals have provisioning data 
Individuals.prov <- length(unique(Prov.rate$BirdID))
#how many nest watches
NWatches<- length(unique(Prov.rate$NWID))

# Gelman-Rubin Criterion estimate ####

#checks whether multiple chains converge on the same posterior distribution 
#not expected to happen by chance but only when data constraining to certain values 
GR.model <- mclapply(1:4, 
                     function(i) 
                       {MCMCglmm(ProvFeedsPerHour ~ HelperNo + Sex + GroupSize + BroodSize + Age_centred + TQ_centred, random=~BirdID + NestID + NestID:NWID + Observer, prior = prior2, nitt=650000, burnin=50000, thin=300, verbose=TRUE, family="poisson",  data=Prov.rate)}, mc.cores=4)

GR.model <- lapply(GR.model, function(m) m$Sol)
GR.model <- do.call(mcmc.list, GR.model)


#plot scale reduction factors
#closer the factor to 1 the better the convergence 
par(mfrow=c(4,2), mar = c(2,2,1,2))
gelman.plot(GR.model, auto.layout = F)
par(mfrow=c(4,2), mar=c(2,2,1,2))
gelman.plot(GR.model, auto.layout=F) # All converge at one after about half iterations 
gelman.diag(GR.model) # All are one

par(mfrow=c(8,2), mar=c(2, 1, 1, 1))
plot(GR.model, ask=F, auto.layout=F)
#traces all seem to follow similar trajectory

#Plotting fixed effects ####

#Set posterior mode and confidence intervals to own dataframe
FixedEffects <- data.frame(posterior.mode(Prov_rate_MCMC_model1$Sol), HPDinterval(Prov_rate_MCMC_model1$Sol))
FixedEffects <- tibble::rownames_to_column(FixedEffects, "term")
FixedEffects <- rename(FixedEffects, estimate = posterior.mode.Prov_rate_MCMC_model1.Sol.)

FixedEffects$term[FixedEffects$term == "HelperNo"] <- "Number of Helpers*"
FixedEffects$term[FixedEffects$term == "Sex1"] <- "Male"
FixedEffects$term[FixedEffects$term == "GroupSize"] <- "Group Size"
FixedEffects$term[FixedEffects$term == "BroodSize"] <- "Brood Size*"
FixedEffects$term[FixedEffects$term == "TQ_centred"] <- "Territory Quality"
FixedEffects$term[FixedEffects$term == "Age_centred"] <- "Age*"

FixedEffects$estimate <- as.numeric(FixedEffects$estimate)
FixedEffects$lower <- as.numeric(FixedEffects$lower)
FixedEffects$upper <- as.numeric(FixedEffects$upper)

FixedEffects$term <- factor(FixedEffects$term, level = c('(Intercept)', 'Number of Helpers*', 'Male', 'Group Size', 'Brood Size*', 'Age*', 'Territory Quality'))

library(wesanderson)
zissou1 <- wes_palette("Zissou1", 7, type = "continuous")
Fox1 <- wes_palette("FantasticFox1", 7, type = "continuous")
FixedEffects %>%
  ggplot(aes(y = term, x = estimate, xmin = lower, xmax = upper, colour = term)) +
  theme_bw() +
  geom_pointinterval(size = 2) +
  geom_vline(xintercept=c(0), linetype="dotted", size = 1) + scale_colour_manual(values = Fox1) + theme(legend.position= "none") + xlab("Provisioning Rate Parameter Estimates") + ylab("")
                   
# Look at males and Females individually to see if repeatability varies between sexes ####
#Note these models do not run currently due to priors not being strong enough
#Subset data into two new data frames with males and females 
Male.PR <- subset(Prov.rate, Sex == 1)
Female.PR <- subset(Prov.rate, Sex == 0)

#Run same model as before 
#Males
Male_model_Prov_rate <- MCMCglmm(ProvFeedsPerHour ~ HelperNo + GroupSize + BroodSize + Age, random=~BirdID + NestID:NWID + Observer, prior = prior2, nitt=430000, burnin=30000, thin=200, verbose=TRUE, family="poisson",  data= Male.PR)

summary(Male_model_Prov_rate)
plot(Male_model_Prov_rate$Sol)
plot(Male_model_Prov_rate$VCV)

#quick repeatability
Male.rep.Prov.rate <- (Male_model_Prov_rate$VCV[, "BirdID"])/ (Male_model_Prov_rate$VCV[, "BirdID"]+ Male_model_Prov_rate$VCV[,"NestID"] + Male_model_Prov_rate$VCV[,"NestID:NWID"] + Male_model_Prov_rate$VCV[,"Observer"] + Male_model_Prov_rate$VCV[,"units"]+log(1/exp(Male_model_Prov_rate$Sol[,"(Intercept)"])+1))

posterior.mode(Male.rep.Prov.rate) #0.101
HPDinterval(Male.rep.Prov.rate)


#Females
Female_model_Prov_rate <- MCMCglmm(ProvFeedsPerHour ~ HelperNo + GroupSize + BroodSize + Age, random=~BirdID + NestID:NWID + Observer, prior = prior2, nitt=430000, burnin=30000, thin=200, verbose=TRUE, family="poisson",  data= Female.PR)

summary(Female_model_Prov_rate)
plot(Female_model_Prov_rate$Sol)
plot(Female_model_Prov_rate$VCV)

#quick repeatability
Female.rep.Prov.rate <- (Female_model_Prov_rate$VCV[, "BirdID"])/ (Female_model_Prov_rate$VCV[, "BirdID"]+ Female_model_Prov_rate$VCV[,"NestID"] + Female_model_Prov_rate$VCV[,"NestID:NWID"] + Female_model_Prov_rate$VCV[,"Observer"] + Female_model_Prov_rate$VCV[,"units"]+log(1/exp(Female_model_Prov_rate$Sol[,"(Intercept)"])+1))

posterior.mode(Female.rep.Prov.rate)
HPDinterval(Female.rep.Prov.rate)

# Non Bayesian Models ####

#Straight forward standard repeatability estimate using rptR and controlling for nest watch duration
#model does not converge
# gives repeatability score of 0.078, CI 0-0.15
rptPoisson(ProvFeedsPerHour ~ HelperNo + Sex + GroupSize + BroodSize + Age_centred + TQ_centred + (1|BirdID) + (1|NestID) + (1|NestID:NWID) + (1|Observer), grname = "BirdID", data = Prov.rate, nboot = 1000, npermut = 1000)


#repeatability using glmer 
#same fixed and random effects as MCMC model
Prov_Rate_glmer_model1 <- glmer(ProvFeedsPerHour ~ HelperNo + Sex + GroupSize + BroodSize + Age_centred + TQ_centred + (1|BirdID) + (1|NestID) + (1|NestID:NWID) + (1|Observer), family = poisson, data = Prov.rate)

summary(Prov_Rate_glmer_model1) # does not converge 

#checking singularity
#If some of the constrained parameters of the random effects (theta) are on the boundary)
#if singular than should be close to one
#but this is not (0.1607)
tt <- getME(Prov_Rate_glmer_model1, "theta")
ll <- getME(Prov_Rate_glmer_model1, "lower")
min(tt[ll==0])

#Check gradient calculations 
derivs1 <- Prov_Rate_glmer_model1@optinfo$derivs
sc_grad1 <- with(derivs1,solve(Hessian,gradient))
max(abs(sc_grad1))
max(pmin(abs(sc_grad1),abs(derivs1$gradient)))    

#increase the number of iterations 
ss <- getME(Prov_Rate_glmer_model1, c("theta", "fixef"))
Prov_Rate_glmer_model2 <- update(Prov_Rate_glmer_model1, start = ss, control = glmerControl(optCtrl = list(maxfun=2e4)))  
    
summary(Prov_Rate_glmer_model2)


#calculate repeatability from glmer
#extract varaince components
variance_components <- get_variance(Prov_Rate_glmer_model2)
#extract estimate for the intercept 
#(intercept) is the first in the vector

glmer_fixed <- as.numeric(fixef(Prov_Rate_glmer_model2))
FE_names <- c("Intercept", "HelperNo", "Sex", "GroupSize", "BroodSize", "Age", "TQ")
glmer_fixed_estimates <- tibble(glmer_fixed, FE_names) %>% 
  filter(FE_names == "Intercept")

#calculate variance from birdID divided by total phenotypic variance 
glmer_repeatability <- (variance_components$var.intercept["BirdID"])/(variance_components$var.random + variance_components$var.residual + log(1/exp(glmer_fixed_estimates$glmer_fixed)+1))

glmer_repeatability #very similar to that produced by MCMC glmm


