#### Adult Density Model ###

### Initialize Workspace -------------------------------------------------------
rm(list = ls())

library(tidyverse) #because... always
library(lme4) #models
library(DHARMa) #model diagnostics

dat <- read.csv("Data/ModDataTraits.csv")

dat$Census <- as.factor(dat$Census)

contrasts(dat$Census) <- contr.sum(5)

dat1 <- filter(dat, SumBA > 0)

### Model

sumBA.mod <- glmer(Symptomatic/N.obs ~ scale(SumBA)  +
                        Census +
                     scale(MAP) +
                        (1|Site/Transect) + (1|Sp),
                      data = dat,
                      weights = N.obs,
                      family = binomial,
                      control = glmerControl(optimizer = "bobyqa"))

summary(sumBA.mod)

testDispersion(sumBA.mod)
simulationOutput <- simulateResiduals(fittedModel = sumBA.mod, plot = F)
plot(simulationOutput)
testZeroInflation(simulationOutput) #this model looks really good.
#As adult abundance increases, their seedlings are more likely to become sick.

saveRDS(sumBA.mod, "modeloutput/sumBA_mod") #good
