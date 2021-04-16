# 16 April 2021
# Samantha Rumschlag
# Generate dataset to be used in macroinvertebrate biodiversity analyses

# load package
library(StreamData)
library(tidyverse)

# read in Devin generated ecoregion table (includes EPA US Level III Ecoregions
# and EPA NARS Level III Ecoregions)
ecoregions <- read.csv("EcoRegionGroups.csv",
                       header = T,
                       stringsAsFactors = FALSE) 

# read in site info - this will be used to join ecoregions contains 
# EPA US Level III Ecoregions for each site
site.info <- read.csv("20201217.0749.SiteInfo.csv", 
                      header = T,
                      stringsAsFactors = FALSE,
                      colClasses = c("SiteNumber" = "character")) %>%
  select(SiteNumber, Ecoregion_US_L3CODE) %>%
  #join NARS
  left_join(ecoregions, by = "Ecoregion_US_L3CODE") %>%
  select(SiteNumber, Ecoregion_NARS)

# set seed to generate the same dataset every time
# WE NEED TO FIX THIS BECAUSE DATASET IS DIFFERENT EVERY TIME...
set.seed(101421)

# generate initial dataset
dat <- getInvertData(
  dataType = "occur",
  taxonLevel = "Genus",
  taxonFix = "remove",
  program = "National Water Quality Assessment",
  lifestage = FALSE,
  abunMeasure = "abundance",
  rarefy = TRUE) 

dat_alpha_all <- dat %>% 
  # exclude sites in Alaska and Hawaii
  filter(!(StateFIPSCode %in% c("02", "15"))) %>%
  # exclude sites collected after 2017
  filter(CollectionYear < 2017) %>%
  # create column for site-level richness
  rowwise() %>%
  mutate(alpha = sum(c_across(30:ncol(dat))>0)) %>% 
  ungroup() %>%
  # make a year column - 1993 = 1
  mutate(Year = CollectionYear - min(CollectionYear) + 1) %>%
  # make a "continuous" year column e.g. 1993 Jan 1 = 1 + 1/365
  # this could be used in richness analyses
  mutate(YearCont = Year + (CollectionDayOfYear / 365)) %>%
  # create column that is number of times a site was sampled through years
  group_by(SiteNumber) %>%
  mutate(n = length(unique(Year))) %>%
  ungroup() %>%
  # add ecoregion data
  left_join(site.info, by = "SiteNumber") %>%
  # remove EcoRegion - Year combos where only 1 site was measured 
  mutate(EcoRegion_Year = paste(Ecoregion_NARS, Year, sep = "_")) 

dat_alpha_all <- dat_alpha_all %>%
  filter(EcoRegion_Year %in% (dat_alpha_all %>%
            group_by(EcoRegion_Year) %>%
            summarize(n = n()) %>%
            ungroup() %>%
            filter(n > 1))$EcoRegion_Year) 

# remove EcoRegions that are only represented in one year across the dataset
# count numer of observations within EcoRegion - Year combos
# all ecoregions are measured in 8 or more years
#
# dat_alpha_all %>%
#  group_by(Ecoregion_NARS) %>%
#  summarize(tal = length(unique(Year))) %>%
#  ungroup()

#
# calculate diversity metrics with all sites, including those that were only 
# sampled one or twice
#

# calculate alpha.bar (mean richness among sites within an ecoregion and year)
dat_alpha.bar <- dat_alpha_all %>%
  # for each Ecoregion_NARS & Year combo, take the mean alpha value
  group_by(Ecoregion_NARS, Year) %>%
  summarize(alpha.bar = mean(alpha)) %>%
  ungroup()

# calculate gamma (total number of genera within an ecoregion and year)
dat_gamma <- dat_alpha_all %>%
  # pivot datset long across the genera columns
  pivot_longer(cols = 30:ncol(dat), 
               names_to = "genera", 
               values_to = "presence") %>%
  # remove genera not present at sites
  filter(presence == 1) %>%
  # for each Ecoregion_NARS & Year combo, count number of genera
  group_by(Ecoregion_NARS, Year) %>%
  summarize(gamma = length(unique(genera))) %>%
  ungroup()

#join datasets and calculate multiplicative & proportional beta
dat_div_all <- dat_alpha.bar %>%
  left_join(dat_gamma, by=c("Ecoregion_NARS", "Year")) %>%
  mutate(mult.beta = gamma/alpha.bar,
         prop.beta = 1 - (alpha.bar / gamma))
  
#
# calculate diversity metrics for only sites that are samples 3 or more times
#

# calcualte diversity metrics with sites that are sample 3 or more times
dat_alpha_3more <- dat_alpha_all %>%
  filter(n >= 3) %>%
  # remove EcoRegion - Year combos where only 1 site was measured 
  mutate(EcoRegion_Year = paste(Ecoregion_NARS, Year, sep = "_")) 

dat_alpha_3more <- dat_alpha_3more %>%
  filter(EcoRegion_Year %in% (dat_alpha_3more %>%
                                group_by(EcoRegion_Year) %>%
                                summarize(n = n()) %>%
                                ungroup() %>%
                                filter(n > 1))$EcoRegion_Year)  

# remove EcoRegions that are only represented in one year across the dataset
# count numer of observations within EcoRegion - Year combos
# all ecoregions are measured in 6 or more years
#
# dat_alpha_3more %>%
#   group_by(Ecoregion_NARS) %>%
#   summarize(tal = length(unique(Year))) %>%
#   ungroup()

# calculate alpha.bar (mean richness among sites within an ecoregion and year)
dat_alpha.bar_3more <- dat_alpha_3more %>%
  # for each Ecoregion_NARS & Year combo, take the mean alpha value
  group_by(Ecoregion_NARS, Year) %>%
  summarize(alpha.bar = mean(alpha)) %>%
  ungroup()

# calculate gamma (total number of genera within an ecoregion and year)
dat_gamma_3more <- dat_alpha_3more %>%
  # pivot datset long across the genera columns
  pivot_longer(cols = 30:ncol(dat), 
               names_to = "genera", 
               values_to = "presence") %>%
  # remove genera not present at sites
  filter(presence == 1) %>%
  # for each Ecoregion_NARS & Year combo, count number of genera
  group_by(Ecoregion_NARS, Year) %>%
  summarize(gamma = length(unique(genera))) %>%
  ungroup()

#join datasets and calculate multiplicative & proportional beta
dat_div_3more <- dat_alpha.bar_3more %>%
  left_join(dat_gamma_3more, by=c("Ecoregion_NARS", "Year")) %>%
  mutate(mult.beta = gamma/alpha.bar,
         prop.beta = 1 - (alpha.bar / gamma))

#read out the final datasets for analyses
write.csv(dat_alpha_all, "dat_alpha_all.csv")
write.csv(dat_div_all, "dat_div_all.csv")

write.csv(dat_alpha_3more, "dat_alpha_3more.csv")
write.csv(dat_div_3more, "dat_div_3more.csv")

