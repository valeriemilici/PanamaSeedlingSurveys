#Intialize Workspace --------
library(parallel) #allows bootstrapping in parallel
library(stats) #for quantile extraction
library(lme4)
#Create cluster within computer to bootstrap in parallel
cl <- makeCluster(detectCores()) 
clusterEvalQ(cl, library(lme4))

# read in the model
CensoMAP_mod <- readRDS("modeloutput/CensoMAP_mod.RDS")


#get the numbers for MAP.scaled
#dat <- read.csv("Data/ModDataTraits.csv")
#quantile(unique(dat$MAP.scaled), probs = seq(0.05,0.95,0.08))
#range(dat$MAP.scaled)

## Step 1: Bootstrap the model for each census along precip gradient

#Create the bootstrapping function 
pred.fun <- function(.) {
  preddat <- expand.grid(MAP.scaled = seq(-1.168,1.30, len =50),
                         Census =as.factor(c(6,7,8,15,16)))
  predict(., newdata = preddat, re.form = ~0)
}

#Bootstrap the model
CensoMAP.boot <- bootMer(CensoMAP_mod, nsim = 1000, FUN = pred.fun,
                     parallel="snow", ncpus = detectCores(), 
                     cl=cl)
#Save the output
save(CensoMAP.boot, file = "BootOutput/CensoMAP2.boot")

## Step 2: Bootstrap intercept only model to visualize main effect of precip
pred.MAP.fun <- function(CensoMAP_mod) {
model.matrix(~MAP.scaled,
data = expand.grid(MAP.scaled =seq(-1.168,1.30, len =50))) %*%
  fixef(CensoMAP_mod)[-c(3:10)]
}

sims <- bootMer(CensoMAP_mod, FUN = pred.MAP.fun,
                nsim = 1000, parallel = "snow", ncpus = detectCores(),cl = cl)
save(sims, file = "BootOutput/MAP1.boot")

#Remove clustering
stopCluster(cl = cl) 
