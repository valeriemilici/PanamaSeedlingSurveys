#### Seed Dispersal Mechanism Model ###

### Initialize Workspace -------------------------------------------------------
rm(list = ls())

library(tidyverse) #because... always
library(lme4) #models
library(lmerTest)
library(DHARMa) #model diagnostics
library(performance) #VIF and model diagnostics

dat <- read.csv("Data/ModDataTraits.csv")

### How many observations of each dispersal mechanism do I have? --------------

wind <- data_frame(dat$wind, dat$N.obs)

wind <- wind %>% filter(dat$wind > 0)

sum(wind$`dat$N.obs`)

sum(dat$N.obs)


#hardly any obs of explosive (36) or water (1) dispersal. remove from models

### Models --------------------------------------------------------------------

dispall.mod <- glmer(Symptomatic/N.obs ~  bat  + nonvolant_mammal +
                       small_bird + large_bird + wind  + scale(MAP) +
                       as.factor(Census) +
                       (1|Transect/Site) + (1|Sp),
                     data = dat,
                     weights = N.obs,
                     family = binomial,
                     control = glmerControl(optimizer = "bobyqa"))

summary(dispall.mod)

anova(dispall.mod)
testDispersion(dispall.mod)
simulationOutput <- simulateResiduals(fittedModel = dispall.mod, plot = F)
plot(simulationOutput)
#deviation signification but the lines look good given the number of data points. 
#nonvolant_mammal dispersal is associated with an increase in obs. symptoms

