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
densMAP_mod <- readRDS("modeloutput/densMAP.RDS")
distMAP_mod <- readRDS("modeloutput/distMAP.RDS")
#Load Bootstrap Output
load("BootOutput/dens.boot")
load("BootOutput/dist.boot")
#Load Data
dat <- read.csv("Data/mod_data_2.csv")

### Prepare the Data ----------------------------------
sp.prox.var <- dat %>% group_by(Sp) %>%
  summarize(spvar = var(prox)) %>%
  filter(spvar > 0)
sp.prox.var <- pull(sp.prox.var, Sp)

sp.dens.var <- dat %>% group_by(Sp) %>%
  summarize(spvar = var(N.obs)) %>%
  filter(spvar > 0)
sp.dens.var <- pull(sp.dens.var, Sp)

dat.dens <- dat %>% filter(Sp %in% sp.dens.var)

dat.dist <- dat %>% filter( Sp %in% sp.prox.var)
### Density only figure --------------------------------------------------------
dens <- expand.grid(MAP = quantile(unique(dat.dens$MAP), probs = 0.5), 
                    Census = as.factor(8),
                    het.dens = median(dat.dens$het.dens),
                    N.obs = seq(1, 45, 5))

dens$preds <- predict(densMAP_mod, newdata = dens, re.form = ~0) 

dens <- mutate(dens, probs = plogis(preds))

#merging bootstrapped confidence intervals
dens<- data.frame(dens, confint(dens.boot))

names(dens)[7:8] <- c("lwr", "upr") 


## plot fig
densfig <- ggplot(dens) + geom_ribbon(aes(y = probs, x = N.obs,
                                                   ymin = plogis(lwr),
                                                   ymax = plogis(upr)),
                                               alpha = 0.3,
                                               show.legend = T) +
  geom_line(mapping = aes(N.obs, probs),
            show.legend = T, size = 1.2) +
  labs(x = expression(paste("Conspecific Seedling Density ",(m^-2))),
       y = "Disease Incidence") +
  theme_classic() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size= 12))

densfig

ggsave(plot = densfig, filename = "Figures/densfig.jpeg",
      height = 3, width = 4, units = "in")

### Distance only figure ------------------------------------------------------- 

dist <- expand.grid(MAP = quantile(unique(dat.dist$MAP), probs = 0.5), 
                          Census = as.factor(8),
                          prox.h = median(dat.dist$prox.h),
                          prox = quantile(dat.dist$prox,
                                          probs = seq(0.1,0.9,0.1)))



dist$preds <- predict(distMAP_mod, newdata = dist, re.form = ~0) 

dist <- mutate(dist, probs = plogis(preds))

#merging bootstrapped confidence intervals
dist<- data.frame(dist, confint(dist.boot))

names(dist)[7:8] <- c("lwr", "upr") 

dist <- mutate(dist, distance = 1/prox)

## plot fig
distfig <- ggplot(dist) + geom_ribbon(aes(y = probs, x = distance,
                                          ymin = plogis(lwr),
                                          ymax = plogis(upr)),
                                      alpha = 0.3,
                                      show.legend = T) +
  geom_line(mapping = aes(distance, probs),
            show.legend = T, size = 1.2) +
  labs(x = "Nearest Conspecific Adult (m)",
       y = "") +
  ylim(c(0,0.08)) +
  theme_classic() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size= 12),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

distfig

ggsave(plot = distfig, filename = "Figures/densfig.jpeg",
       height = 3, width = 4, units = "in")

## Join Together (Figure 3) ---------------------------------------------------

distdensfig <- densfig + distfig + 
  plot_annotation(tag_levels = "A",
                  theme = theme(plot.tag = element_text(size = 10)))
                                                   
distdensfig
ggsave("Figures/distdensfig.png", plot = distdensfig, 
       height = 4, width = 7, units = "in")


