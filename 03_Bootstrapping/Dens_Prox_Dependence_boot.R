## Initialize Workspace -------
rm(list = ls())
library(parallel) #allows bootstrapping in parallel
library(stats) #for quantile extraction

#Create cluster within computer to bootstrap in parallel
cl <- makeCluster(detectCores()) 
clusterEvalQ(cl, library(lme4))

#Read in the model
densMAP_mod <- readRDS("modeloutput/densMAP.RDS")
distMAP_mod <- readRDS("modeloutput/distMAP.RDS")

### Density and Proximity Dependence Models -------

pred.dens.fun <- function(.) {
  preddens <- expand.grid(MAP = 2311,
                         Census = as.factor(c(8,15)),
                         het.dens = 5,
                         focal.dens = seq(1,50,5))
  predict(., newdata = preddens, re.form = ~0)
}

dens.boot <- bootMer(densMAP_mod, nsim = 1000, FUN = pred.dens.fun,
                     parallel="snow", ncpus = detectCores(), 
                     cl=cl)
save(dens.boot, file = "BootOutput/dens.boot")


pred.dist.fun <- function(.) {
  preddist <- expand.grid(MAP = c(2311),
                          Census = as.factor(c(8,15)),
                          prox.h = 0.05,
                          prox = seq(0.05,1,0.05))
  predict(., newdata = preddist, re.form = ~0)
}

dist.boot <- bootMer(distMAP_mod, nsim = 1000, FUN = pred.dist.fun,
                     parallel="snow", ncpus = detectCores(), 
                     cl=cl)
save(dist.boot, file = "BootOutput/dist.boot")


stopCluster(cl = cl) 
