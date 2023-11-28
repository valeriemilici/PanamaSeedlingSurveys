### Plot the relationship between census and MAP on disease incidence

rm(list = ls())

library(tidyverse) #because it's always needed
library(splines) #allow for bends in the line
library(ggplot2) # because I'm plotting
library(scales) #plogis figure y-axis
library(functional) #for plogis transformation
library(stats) #for quantile extraction
library(patchwork) #for tiling graphs
library(lubridate) # for date objects

#Read in model and data --------------------------------------------------------
CensoMAP_mod <- readRDS("modeloutput/CensoMAP_mod.RDS")
dat <- read.csv("Data/ModDataTraits.csv")

#Make predictions --------------------------------------------------------------
preddat <- expand.grid(MAP.scaled = seq(-1.168,1.30, len =50), 
                      Census = as.factor(c(6,7,8,15,16))) 
preddat$preds <- predict(CensoMAP_mod, newdata = preddat, re.form = ~0) 

#merging bootstrapped confidence intervals
load("BootOutput/CensoMAP2.boot")
preddat<- data.frame(preddat, confint(CensoMAP.boot))

names(preddat)[4:5] <- c("lwr", "upr") 

#Giving Census a month
Census <- c(6,7,8,15,16)

fig.dates <- data.frame(Census,
                        CensusDate = factor(c(
                          "June", "July", "August", "March", "April"), levels = c(
                            "June", "July", "August", "March", "April") ),
                        Season = factor(c("Wet Season", "Wet Season","Wet Season", "Dry Season", "Dry Season"),
                                        levels = c("Wet Season", "Dry Season")  ))#bind it together

#make sure fixed and random effect lines have census date column
preddat <- merge(preddat, fig.dates, by = "Census") 

### Liza's Suggestion ----------------------------------------------------------
MAP <- unique(dat$MAP)

mean(MAP) #2335.6
sd(MAP) #396.3

preddat$MAP <- (preddat$MAP.scaled*sd(MAP)) + mean(MAP)

preddat <- preddat[,-c(1,2,7)]

## Panel for MAP only

load("BootOutput/MAP1.boot")

MAPpreds <- data.frame(t(apply(sims$t, 2, quantile, c(0.5,0.025, 0.975))))
MAPpreds$MAP.scaled<-  seq(-1.168,1.30, len =50)


MAPpreds$MAP <- MAPpreds$MAP.scaled*sd(MAP) + mean(MAP)
colnames(MAPpreds) <- c("preds", "lwr", "upr", "MAP.scaled", "MAP")

MAPpreds <- MAPpreds[,-4]

MAPpreds$CensusDate <- "Mean Effect"

MAPpreds <- MAPpreds[,c(1:3,5,4)]

preddat2<- rbind(preddat, MAPpreds)

preddat2$CensusDate <- factor(preddat2$CensusDate, levels = c("Mean Effect", "June",
                                                              "July", "August",
                                                              "March", "April"))
#Plot it together
CensoMAP <- ggplot() +
  geom_ribbon(preddat2, mapping = aes(x = MAP, y = plogis(preds), ymin = plogis(lwr), ymax = plogis(upr)), 
              alpha = 0.3) +
  geom_line(preddat2, mapping = aes(MAP, plogis(preds))) +
  labs(y = "Disease Incidence", x = "Mean Annual Precipitation (mm)") +
  facet_wrap(~CensusDate, scales = "free") +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        strip.text.x = element_text(size = 10))

CensoMAP

ggsave(filename = "Figures/CensoMAP1.png", plot = CensoMAP, height = 4, width = 6,
       units = "in")
