## Creates dataframe with species traits for analyses that don't include
## conspecific or heterospecific density --------------------------------------

## Initialize Workspace -------------------------------------------------------

rm(list = ls())

library(tidyverse)
library(readxl)
library(stringr)
library(spatstat)

dat <- read.csv("Data/Traits_Data_Raw.csv") 
#census data, includes SumBA, ready for merge with dispersal and sm data

dispersal <- read.csv("Data/dispersal_SJW.csv") #dispersal data
#See Wright et al., 2016 (Ecology;
# https://onlinelibrary.wiley.com/doi/10.1002/ecy.1519) for dispersal data 

sm <- read.csv("Data/Fruit_Seed_masses_20210111.csv") #seed mass data
#See Wright et al., 2010 (Ecology;
# https://doi.org/10.1890/09-2335.1) for seed mass data


### Add in dispersal regime (Wright et al., 2016) -----------------------------

names(dispersal)[2] <- "Sp"

dispersal$Sp <- tolower(dispersal$Sp)

disp <- select(dispersal, c(2, 7:13))

#Join together

datdisp<- left_join(dat, disp, by = "Sp")

## Manual set dispersal types to sum to zero contrasts
datdisp$bat <- ifelse(datdisp$bat == "TRUE", 0.5, -0.5/4)
datdisp$small_bird <- ifelse(datdisp$small_bird == "TRUE", 0.5, -0.5/4)
datdisp$large_bird <- ifelse(datdisp$large_bird == "TRUE", 0.5, -0.5/4)
datdisp$nonvolant_mammal <- ifelse(datdisp$nonvolant_mammal == "TRUE", 0.5, -0.5/4)
datdisp$wind <- ifelse(datdisp$wind == "TRUE", 0.5, -0.5/4)

### Add in SM data (Wright et al. 2010) ---------------------------------------

drymass <- select(sm, c(14,22))
drymass$sp6 <- tolower(drymass$sp6)
names(drymass)[2] <- "Sp"
drymass.clean <- na.omit(drymass)

meanmass <- drymass.clean %>%
  group_by(Sp) %>% 
  summarise("seedmass" = mean(seed_dry))

datmass <- left_join(datdisp, meanmass, by = "Sp")

#Seedmass and sumBA
datScaled <- datmass %>%
  mutate(seedmass.scaled = ((seedmass - mean(seedmass, na.rm = T))/
                              sd(seedmass, na.rm = T)))

### You are done. Save the output for the models! -----------------------------
write.csv(datScaled, file = "Data/ModDataTraits.csv")

