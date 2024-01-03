### Initialize Workspace -------------------------------------------------------
rm(list = ls())
library(tidyverse) #because it's always needed
library(splines) #allow for bends in the line
library(ggplot2) # because I'm plotting
library(scales) #plogis figure y-axis
library(functional) #for plogis transformation
library(stats) #for quantile extraction
library(patchwork) #for tiling graphs
library(lubridate) # for date objects

#Load Model Output
BA.mod <- readRDS("modeloutput/sumBA_mod")
#Load Boot Output
load("BootOutput/BA.boot")
#Load Data
dat <- read.csv("Data/ModDataTraits.csv")

### Plot Figure ----------------------------------------------------------------
#create new data for predictions
sumBA <- expand.grid(SumBA = seq(1500, 12000000, 100000),
                     Census = as.factor(8),
                     MAP = quantile(unique(dat$MAP), probs = 0.5))
#add model predictions
sumBA$preds <- predict(BA.mod, newdata = sumBA, re.form = ~0) 

#convert to data scale
sumBA <- mutate(sumBA, probs = plogis(preds))

#merging bootstrapped confidence intervals
sumBA<- data.frame(sumBA, confint(BA.boot))

names(sumBA)[6:7] <- c("lwr", "upr") 

#create the figure
BAfig <- ggplot(sumBA) + geom_ribbon(aes(y = probs, x = SumBA,
                                          ymin = plogis(lwr),
                                          ymax = plogis(upr)),
                                      alpha = 0.3,
                                      show.legend = T) +
  geom_line(mapping = aes(SumBA, probs),
            show.legend = T, size = 1.2) +
  labs(x = expression(paste("Conspecific Adult Basal Area ",(ha^-1))),
       y = "Disease Incidence") +
  theme_classic() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size= 12))

BAfig

ggsave(plot = BAfig, filename = "Figures/BAfig.png",
       height = 3, width = 4, units = "in")
