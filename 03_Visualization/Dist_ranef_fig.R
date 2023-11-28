rm(list = ls())

library(tidyverse) #because it's always needed
library(ggplot2) # because I'm plotting
library(stats) #for quantile extraction
library(patchwork) #for tiling graphs
library(stringr) #add leading 0s
library(lme4)

dist <- readRDS(file = "modeloutput/distMAP.RDS")
dat <- read.csv("Data/mod_data_2.csv")
###Species effects don't vary for density model

# Add a column with the sp. name and abundance in dataset (Common for both panels)
abund <- dat %>% dplyr::select(Sp, species_abund, genus, species,abund, sumBA) %>%
 slice_sample(n = 1, by = Sp)

abund$abund <- str_pad(abund$abund, 3, pad = "0")
### Distance Random Effect Blup Plot -------------------------------------------

#Step 1: Extract Random Effects and calculate lwr and upr
prox_re <- data.frame(ranef(dist)) #extract random effects from model

prox.sp_re <-filter(prox_re, grpvar == "Sp" & term == "log1p(prox)") #filter for the random slope of prox

prox.sp_re <- prox.sp_re[-(1:2)] #remove unneccesary columns

prox.sp_re <- prox.sp_re %>% mutate( lwr = condval - condsd) %>%
  mutate(upr = condval + condsd)%>%
  dplyr::rename(Sp = grp) #add columns for upper and lower CI values


#Step 2: Join together
prox.sp_re <- inner_join(prox.sp_re, abund, by = "Sp") #Join

#Step 3: Create columns with full species name
prox.sp_re <- prox.sp_re %>%
  mutate(abund.p = paste0("(", abund, ")")) %>%
  replace(is.na(.), 0)

prox.sp_re$sp.name <- paste(prox.sp_re$genus, prox.sp_re$species, sep = " ")
#entecy is missing?

prox.sp_re[15,12] <- "Enterolobium cyclocarpum"

prox.sp_re$obs.dom <- paste(prox.sp_re$sp.name, prox.sp_re$abund.p, sep = " ")
#ok now it looks good

### Make the figure (original) -------------------------------------------------
dist.re <- ggplot(prox.sp_re, aes(condval, reorder(obs.dom, condval)))+
  geom_point() +
  geom_errorbarh(mapping = aes(xmin = lwr, xmax = upr)) +
  geom_vline(xintercept = 0) + 
  scale_y_discrete(guide = guide_axis(check.overlap = T)) +
  labs(x = "Effect of Proximity",
       y = "Species") +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        #legend.title = element_text(size = 14),
       # legend.text = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.y = element_text(face = "italic"))

dist.re

ggsave("Figures/dist_ranef.jpeg", plot = dist.re,
       height = 8, width = 5, units = "in")

### Order by species dominance (the adults) ------------------------------------

dist.re2 <- ggplot(prox.sp_re, aes(condval, reorder(obs.dom, as.numeric(sumBA))))+
  geom_point() +
  geom_errorbarh(mapping = aes(xmin = lwr, xmax = upr)) +
  geom_vline(xintercept = 0) + 
  scale_y_discrete(guide = guide_axis(check.overlap = T)) +
  labs(x = "Effect of Proximity",
       y = "Species") +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        #legend.title = element_text(size = 14),
       # legend.text = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.y = element_text(face = "italic"))

dist.re2

ggsave("Figures/dist_ranef2.png", plot = dist.re2,
       height = 8, width = 5, units = "in")

### Order by species dominance (seedling observations) ------------------------
prox.sp_re <- arrange(prox.sp_re, desc(abund)) %>%
  mutate(abund.p = paste0("(", abund, ")"))

prox.sp_re$sp.name <- paste(prox.sp_re$genus, prox.sp_re$species, sep = " ")
#entecy is missing?

prox.sp_re[15,10] <- "Enterolobium cyclocarpum"



prox.sp_re$obs.dom <- paste(prox.sp_re$sp.name, prox.sp_re$abund.p, sep = " ")
#ok now it looks good

dist.re3 <- ggplot(prox.sp_re, aes(condval, reorder(obs.dom, as.numeric(abund))))+
  geom_point() +
  geom_errorbarh(mapping = aes(xmin = lwr, xmax = upr)) +
  geom_vline(xintercept = 0) + 
  scale_y_discrete(guide = guide_axis(check.overlap = T)) +
  labs(x = "Effect of Proximity",
       y = "Species") +
  theme_classic() +
  theme(axis.title = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.text.y = element_text(face = "italic"))

dist.re3

ggsave("Figures/dist_ranef3.jpeg", plot = dist.re3,
       height = 8, width = 5, units = "in")




