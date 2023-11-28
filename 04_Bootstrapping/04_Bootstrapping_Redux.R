## Bootstrapping Code 
rm(list = ls())
library(parallel) #allows bootstrapping in parallel
library(stats) #for quantile extraction

#Create cluster within computer to bootstrap in parallel
cl <- makeCluster(detectCores()) 
clusterEvalQ(cl, library(lme4))

#########################
######## NCDD * MAP #####
########################
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

################################
######### Dispersal Model ######
################################

pred.disp.fun <- function(.) {
  preddisp <- expand.grid(Census = as.factor(15),
                          bat = c(0.5, -0.5/3),
                          nonvolant_mammal = c(0.5, -0.5/3),
                          small_bird = c(0.5, -0.5/3),
                          large_bird = c(0.5, -0.5/3))
  predict(., newdata = preddisp, re.form = ~0)
}

disp.boot <- bootMer(dispall.mod, nsim = 1000, FUN = pred.disp.fun,
                     parallel="snow", ncpus = detectCores(), 
                     cl=cl)
save(disp.boot, file = "BootOutput/disp.boot")
 

############################
######### SumBA ###########
###########################

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

### Density * SumBA -----------------------------------------------------------

densBA.fun <- function(.) {
  densBA <- expand.grid(MAP = 2311,
                          Census = as.factor(8),
                          het.dens = 5,
                          N.obs = c(1,4),
                        sumBA = seq(1500, 12000000, 100000))
  predict(., newdata = densBA, re.form = ~0)
}

densBA.boot <- bootMer(BAdens.mod, nsim = 1000, FUN = densBA.fun,
                     parallel="snow", ncpus = detectCores(), 
                     cl=cl)
save(densBA.boot, file = "BootOutput/densBA.boot")



### Distance * SumBA -----------------------------------------------------------
distBA.fun <- function(.) {
  distBA <- expand.grid(MAP = 2311,
                        Census = as.factor(8),
                        prox.h = 2.917301,
                        prox = c(0.02487738,0.08802606, 0.44538877),
                        sumBA = seq(1500, 12000000, 100000))
  predict(., newdata = distBA, re.form = ~0)
}
distBA.boot <- bootMer(distSumBA_mod, nsim = 1000, FUN = distBA.fun,
                       parallel="snow", ncpus = detectCores(), 
                     cl=cl)
save(distBA.boot, file = "BootOutput/distBA.boot")




stopCluster(cl = cl) 
