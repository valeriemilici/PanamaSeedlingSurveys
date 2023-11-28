#### Seed Mass Model ###

### Initialize Workspace -------------------------------------------------------
rm(list = ls())

library(tidyverse) #because... always
library(lme4) #models
library(DHARMa) #model diagnostics
library(glmmTMB)

dat <- read.csv("Data/ModDataTraits.csv")

### Model

seedmass.mod <- glmer(Symptomatic/N.obs ~ seedmass.scaled  +
                       as.factor(Census) +
                       (1|Transect/Site) + (1|Sp),
                     data = dat,
                     weights = N.obs,
                     family = binomial,
                     control = glmerControl(optimizer = "bobyqa"))

summary(seedmass.mod)
testDispersion(seedmass.mod)
simulationOutput <- simulateResiduals(fittedModel = seedmass.mod, plot = F)
plot(simulationOutput)
testZeroInflation(simulationOutput)
#seedmass alone does not explain anything. 


