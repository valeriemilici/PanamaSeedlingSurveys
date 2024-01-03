#### General Effects Model ####

### Initialize workspace ------------------------------------------------------
rm(list = ls())

library(tidyverse) #because... always
library(lme4) #models
library(DHARMa) #model diagnostics
library(performance) #VIF and model diagnostics
library(lmerTest)

dat <- read.csv("Data/ModDataTraits.csv")

### Model ---------------------------------------------------------------------
dat$Census <- factor(dat$Census,
                     levels = c( "8", "15", "6", "16", "7"))
contrasts(dat$Census)<- contr.sum(5)


CensoMAP_mod <- glmer(Symptomatic/N.obs ~ 
                        MAP.scaled*as.factor(Census)+ 
                        (1|Site/Transect) + (1|Sp),
                      family = binomial,
                      weights = N.obs, 
                      data = dat,  
                      control = glmerControl(optimizer = "bobyqa"))

summary(CensoMAP_mod)

car::Anova(CensoMAP_mod)
saveRDS(CensoMAP_mod, file = "modeloutput/CensoMAP_mod.RDS")