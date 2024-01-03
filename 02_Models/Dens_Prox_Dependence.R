#### CNDDxMAP Model ###

### Initialize Workspace -------------------------------------------------------
rm(list = ls())

library(tidyverse) #because... always
library(lme4) #models
library(DHARMa) #model diagnostics
library(performance) #VIF and model diagnostics
library(pbkrtest) #LRT in dist model
library(gamm4) #nonlinear model testing
library(afex) #automatically uses all model convergence optimizers

dat <- read.csv("Data/NCDD_Data.csv")
#Data found in Zenodo Repository 10.5281/zenodo.10456096. 

### Prepare the Data ----------------------------------
sp.prox.var <- dat %>% group_by(Sp) %>%
  summarize(spvar = var(prox)) %>%
  filter(spvar > 0)
sp.prox.var <- pull(sp.prox.var, Sp)

sp.dens.var <- dat %>% group_by(Sp) %>%
  summarize(spvar = var(N.obs)) %>%
  filter(spvar > 0)
sp.dens.var <- pull(sp.dens.var, Sp)

#filter data to fit selection criteria for both density and proximity
#to use in analysis that includes both variables.
sp.both.var <- dat %>% group_by(Sp) %>%
  summarise(spvar = var(N.obs) & var(prox)) %>%
filter(spvar > 0)

sp.both.var <- pull(sp.both.var, Sp)


dat.dens <- dat %>% filter(Sp %in% sp.dens.var)

dat.dist <- dat %>% filter( Sp %in% sp.prox.var)

dat.both <- dat %>% filter(Sp %in% sp.both.var)

dat.dens$Census <- as.factor(dat.dens$Census)
dat.both$Census <- as.factor(dat.both$Census)

contrasts(dat.both$Census) <- contr.sum(5)
contrasts(dat.dens$Census)<- contr.sum(5)

### Basic CDD Models ---------------------------------------------------------

# Are density and proximity correlated?

cor.test(dat$N.obs, dat$prox, method = "pearson") 
#no, not really! (correlation = 0.07)
# this correlation is significant, but it is slight and unlikely to lead to
# variance inflation/ identifiability issues

## Density + Proximity * MAP -------------------------

#Note: This model tests whether both density and proximity can be used in the
# same model, because they are not correlated. Tests show that the model has
# a degenerate Hessian when both are included, which causes issues in the
# estimation of SE. Density and Proximity should be modeled separately.

both_mod <- glmer(Symptomatic/N.obs ~ 
                    scale(MAP) *( scale(log1p(het.dens)) +
                                    scale(log(N.obs)) +
                                    scale(log1p(prox.h)) +
                                    scale(log1p(prox)) )+
                    as.factor(Census)+ 
                    (1|Site/Transect) + (log(N.obs)||Sp) + (log1p(prox)||Sp),
                  family = binomial,
                  weights = N.obs, 
                  data = dat.both, 
                  control = glmerControl(optimizer = "bobyqa"))

summary(both_mod)
#model fails to converge. 

allFit(both_mod) #will an alternative optimizer get it to converge?
#alternative optimizers are ok, but model still fails to converge and has
#degenerate Hessian. 

## Density x MAP -------------------------------
densMAP_mod <- glmer(cbind(Symptomatic, N.obs - Symptomatic) ~ 
                       scale(MAP) * (scale(log1p(het.dens)) +
                      scale(log(N.obs))) + as.factor(Census)+ 
                        (1|Site/Transect) + (log(N.obs)||Sp),
                     family = binomial,
                     data = dat.dens, 
                     control = glmerControl(optimizer = "bobyqa"))

simres <- simulateResiduals(densMAP_mod)
plotResiduals(simres, form = dat.dens$N.obs)
summary(densMAP_mod)


saveRDS(densMAP_mod, "modeloutput/densMAP_mod") #good

## Proximity x MAP ------------------------------
dat.dist$Census <- as.factor(dat.dist$Census)
contrasts(dat.dist$Census) <- contr.sum(5)

distMAP_mod <- glmer(Symptomatic/N.obs ~ 
                       (scale(MAP)  * (scale(log1p(prox.h)) + scale(log1p(prox)))) +
                       as.factor(Census) + 
                       (1|Site/Transect) + ((log1p(prox))||Sp),
                     family = binomial,
                     weight = N.obs, 
                     data = dat.dist, 
                     control = glmerControl(optimizer="bobyqa"))

summary(distMAP_mod)

saveRDS(file = "modeloutput/distMAP.RDS",object =  distMAP_mod)

# Density Random Slope LRT ----------------------
densMAP_mod_i <- glmer(Symptomatic/N.obs ~ 
                       scale(MAP) * (scale(log1p(het.dens)) +
                                       scale(log(N.obs))) + as.factor(Census)+ 
                       (1|Site/Transect) + (1|Sp),
                     family = binomial,
                     weights = N.obs, 
                     data = dat.dens, 
                     control = glmerControl(optimizer = "bobyqa"))

anova(densMAP_mod, densMAP_mod_i)


saveRDS(file = "modeloutput/densMAP.RDS",object =  densMAP_mod)
# evidence for CNDD, and more disease where it's drier, but no interaction.

# Proximity Random Slope LRT -----
distMAP_mod_i <- glmer(Symptomatic/N.obs ~ 
                       (scale(MAP)  * (scale(log1p(prox.h)) + scale(log1p(prox)))) +
                       as.factor(Census) + 
                       (1|Site/Transect) + ((1|Sp)),
                     family = binomial,
                     weight = N.obs, 
                     data = dat.dist, 
                     control = glmerControl(optimizer="bobyqa"))

anova(distMAP_mod, distMAP_mod_i)

## Nonlinear responses --------------------
# Effect of Density GAMM ----------------------
densMAP_gam <- gamm4(cbind(Symptomatic, N.obs - Symptomatic) ~ 
                       scale(log1p(het.dens)) +
                       scale(log1p(het.dens)):scale(MAP) +
                       s(scale(MAP), scale(log(N.obs)), bs ="tp", k = 20) +
                       as.factor(Census), 
                     random = ~ (1|Site/Transect) + (log(N.obs)||Sp),
                     family = binomial,
                     data = dat.dens,
                     REML = F)

summary(densMAP_gam$gam)
#EDF is 1.18, p = 0.12. The effect of density is
#essentially linear.
summary(densMAP_gam$mer)


compare_performance(densMAP_gam$mer, densMAP_mod, rank = T)
#the linear model has more support than the GAM. We should interpret the data 
#using the linear model.
#lme4 model performance score = 75%, gamm model performance score = 25%
AIC(densMAP_gam$mer, densMAP_mod)


# Effect of Proximity GAMM ------------------------
distMAP_gam <- gamm4(cbind(Symptomatic, N.obs - Symptomatic) ~ 
                      scale(log1p(prox.h))+
                      scale(log1p(prox.h)):scale(MAP) +
                      s(scale(MAP), scale(log1p(prox)), bs = "tp", k = 20) +
                      as.factor(Census), 
                      random = ~ (1|Site/Transect) + ((log1p(prox))||Sp),
                     family = binomial,
                     data = dat.dist, 
                     REML = F)

summary(distMAP_gam$gam)
#EDF is 1, p = 0.002, which means that the effect of proximity is linear

summary(distMAP_gam$mer)

compare_performance(distMAP_mod, distMAP_gam$mer, rank = T)
AIC(distMAP_gam$mer, distMAP_mod)
#both tests indicate that the lme4 model is better than the gamm.
#dAIC favoring lme4 model = 1.3
#lme4 performance score = 75%, gam performance score = 25%
