#### MAP*SeedMass Model ###

### Initialize Workspace -------------------------------------------------------
rm(list = ls())

library(tidyverse) #because... always
library(lme4) #models
library(DHARMa) #model diagnostics

dat <- read.csv("Data/ModDataTraits.csv")

### Model

MAPxMass.mod <- glmer(Symptomatic/N.obs ~ seedmass.scaled*MAP.scaled  +
                        as.factor(Census) +
                        (1|Site/Transect) + (1|Sp),
                      data = dat,
                      weights = N.obs,
                      family = binomial,
                      control = glmerControl(optimizer = "bobyqa"))

summary(MAPxMass.mod)
#seedmass and MAP do not interact