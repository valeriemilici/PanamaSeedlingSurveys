## Initialize Workspace ---------
rm(list = ls())
library(parallel) #allows bootstrapping in parallel
library(stats) #for quantile extraction

#Create cluster within computer to bootstrap in parallel
cl <- makeCluster(detectCores()) 
clusterEvalQ(cl, library(lme4))


# Adult Abundance (SumBA) Model ------

pred.BA.fun <- function(.){
  predBA <- expand.grid(SumBA = seq(1500, 12000000, 100000),
                        Census = as.factor(8),
                        MAP = 2311)
  
  predict(., newdata = predBA, re.form = ~0)
}

BA.mod <- readRDS("modeloutput/sumBA_mod")
BA.boot <- bootMer(BA.mod, nsim = 1000, FUN = pred.BA.fun,
                   parallel = "snow", ncpus = detectCores(),
                   cl = cl)

save(BA.boot, file = "BootOutput/BA.boot")

stopCluster(cl = cl) 