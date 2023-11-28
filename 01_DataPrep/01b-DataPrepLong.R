## Creates dataframe with species traits for analyses that don't include
## conspecific or heterospecific density --------------------------------------

## Initialize Workspace -------------------------------------------------------

rm(list = ls())

library(tidyverse)
library(readxl)
library(stringr)
library(spatstat)

datst <-read.csv("Data/mod_data_long.csv") #census data 

spnames <- read_excel("Data/speciescode_Key.xlsx") #species name info

dispersal <- read.csv("Data/dispersal_SJW.csv") #dispersal data

sm <- read.csv("Data/Fruit_Seed_masses_20210111.csv") #seed mass data

adult <- read.csv("Data/siteadults.csv") #adults at each site data
### Add formal species names to data ------------------------------------------

spp <- select(spnames, c(1,2,3,5))

#separate column 1 into two columns

gs <- data.frame(str_split_fixed(spp$Species, " ", 2))

names(gs) <- c("genus", "species")

#merge together

spp <- cbind(spp, gs)

spp <- select(spp, c(5,6,1,2,3,4))

names(spp)[6] <- "Sp"

#Join data and species information together

datsp <- left_join(datst, spp, by = "Sp")

### Add in dispersal regime (Wright et al., 2016) -----------------------------

names(dispersal)[2] <- "Sp"

dispersal$Sp <- tolower(dispersal$Sp)

disp <- select(dispersal, c(2, 7:13))

#Join together

datdisp<- left_join(datsp, disp, by = "Sp")

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

### Add in total adult basal area at each site --------------------------------
sumBA <- adult %>%
  rename(Sp = Mnemonic) %>%
  mutate(BA = pi*((DBH/2)^2)) %>%
  group_by(Site, Sp) %>%
  summarise("SumBA" = sum(BA)) #for the record, this worked on the 1st try!

datBA <- left_join(datmass, sumBA, by = c("Site", "Sp"))

### Scale the continuous variables --------------------------------------------

#MAP (must scale them as unique numbers)
MAP <- data.frame(MAP = as.numeric(unique(datBA$MAP)),
                 Site = unique(datBA$Site))
MAP$MAP.scaled <- ((MAP$MAP - mean(MAP$MAP))/ sd(MAP$MAP))
MAP <- MAP[,-1]

datBA <- left_join(datBA, MAP, by = "Site")

#Seedmass and sumBA
datScaled <- datBA %>%
  mutate(seedmass.scaled = ((seedmass - mean(seedmass, na.rm = T))/sd(seedmass, na.rm = T)),
         SumBA.scaled = ((SumBA - mean(SumBA, na.rm = T))/sd(SumBA, na.rm = T)))

# Adding nearest adult column for dispersal information

adult$site.sp <- paste(adult$Site, adult$Mnemonic, sep = "_")

# Nearest adult
closest_adult<- vector()
datScaled$closest_adult<- NA
a<- vector()

for (i in 1:nrow(datScaled)) {
  
  a<- adult[which(as.character(adult$site.sp) == as.character(datScaled$site.sp[i])),]
  
  closest_adult <- min(crossdist.default(datScaled$x[i], datScaled$y[i], a$PX, a$PY))
  
  datScaled$closest_adult[i]<- closest_adult
}



### You are done. Save the output for the models! -----------------------------
write.csv(datScaled, file = "Data/ModDataTraits.csv")

