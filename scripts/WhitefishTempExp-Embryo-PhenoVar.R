#### CLEAR THE ENVIRONMENT FIRST -----------------------------------------------------------------

#rm(list = ls(all.names = TRUE))


#### LOAD PACKAGES  ------------------------------------------------------------------------------

library(tidyverse)
library(readxl)
library(data.table)
library(fullfact)
library(parallel)


#### LOAD INCUBATION TEMPERATURE DATA ------------------------------------------------------------

ADD <- read.csv("data/Artedi-Temperature-ADD-2020.csv", header = TRUE) %>% 
  dplyr::select(population, temperature, ADD) %>% 
  group_by(population, temperature) %>% 
  mutate(dpf = 1:n())


#### LOAD HATCHING DATA --------------------------------------------------------------------------

hatch.USA <- read_excel("data/Coregonine-Temperature-Experiment-NA-Hatch.xlsx", sheet = "2020HatchingData") %>% 
  filter(is.na(notes) | notes != "empty well") %>% 
  filter(block != "A" | population != "superior") %>% 
  mutate(eye = as.numeric(eye),
         hatch = as.numeric(hatch)) %>% 
  filter(!is.na(eye), !is.na(hatch)) %>% 
  left_join(ADD) %>% 
  dplyr::select(population, species, family, male, female, block, temperature, eye, hatch, dpf, ADD)

hatch.Finland <- read_excel("data/Coregonine-Temperature-Experiment-FI-Hatch.xlsx", sheet = "2019HatchingData") %>% 
  mutate(premature = 0) %>% 
  dplyr::select(population, species, family, male, female, block, temperature, eye, hatch, dpf, ADD)

## Combine all populations and years
hatch <- bind_rows(hatch.USA, hatch.Finland) %>% 
  mutate(population = factor(population, levels = c("konnevesi", "superior", "ontario"), ordered = TRUE),
         temperature = factor(temperature, ordered = TRUE, 
                              levels = c(2, 2.2, 4.0, 4.4, 6.9, 8, 8.9),
                              labels = c("2.0", "2.2", "4.0", "4.4", "6.9", "8.0", "8.9")),
         female = factor(female, levels = seq(1, 12, 1),
                         labels = c("F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9", "F10", "F11", "F12")),
         male = factor(male, levels = seq(1, 16, 1),
                       labels = c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9", "M10", "M11", "M12", "M13", "M14", "M15", "M16")),
         # Create a variable with population and species combined
         population = factor(interaction(population, species), ordered = TRUE,
                             levels = c("konnevesi.albula", "konnevesi.lavaretus", "superior.artedi", "ontario.artedi"),
                             labels = c("LK-Vendace", "LK-Whitefish", "LS-Cisco", "LO-Cisco")),
         group = interaction(population, temperature)) %>% 
  rename(sire = male, dam = female)

## Clean up environment
rm(hatch.USA, hatch.Finland, ADD)


#### FILTER TO EACH TRAITS' DATASET --------------------------------------------------------------

## filter to only eyed embryos
hatch.survival <- hatch %>% filter(eye != 0)

## filter to only hatched embryos
hatch.dpf <- hatch %>% filter(!is.na(dpf), hatch == 1)

## filter to only hatched embryos
hatch.ADD <- hatch %>% filter(!is.na(ADD), hatch == 1)


#### STATISTICAL ANALYSIS - GENERATE OBSERVED VARIANCES ------------------------------------------

## Embryo Survival
phenoVar.survival.obs <- do.call(rbind, lapply(as.character(unique(hatch.survival$group)), function(grp) {
  ## Filter to only a single group
  data.group <- hatch.survival %>% filter(group == grp) %>% 
      select(family, dam, sire, block, hatch)
  
  obs.survival <- observGlmer(observ = data.group, dam = "dam", sire = "sire", response = "hatch",
                               fam_link = binomial(logit))
  
  obs.survival.df <- data.frame(group = substr(grp, 1, nchar(grp)-4),
                                temperature = as.numeric(substr(grp, nchar(grp)-2, nchar(grp))),
                                dam.var = obs.survival$random[3,2],
                                dam.p = obs.survival$random[3,7],
                                dam.perc = obs.survival$random[3,3],
                                sire.var = obs.survival$random[2,2],
                                sire.p = obs.survival$random[2,7],
                                sire.perc = obs.survival$random[2,3],
                                dam.sire.var = obs.survival$random[1,2],
                                dam.sire.p = obs.survival$random[1,7],
                                dam.sire.perc = obs.survival$random[1,3],
                                residual.var = obs.survival$other[1,2],
                                residual.perc = obs.survival$other[1,3]) %>% 
    mutate_if(is.numeric, round, 4)
})) %>% mutate(trait = "survival")

## DPF
phenoVar.dpf.obs <- do.call(rbind, lapply(as.character(unique(hatch.dpf$group)), function(grp) {
  ## Filter to only a single group
  data.group <- hatch.dpf %>% filter(group == grp) %>% 
      select(family, dam, sire, block, dpf)
    
  obs.dpf <- observLmer(observ = data.group, dam = "dam", sire = "sire", response = "dpf")
  
  obs.dpf.df <- data.frame(group = substr(grp, 1, nchar(grp)-4),
                           temperature = as.numeric(substr(grp, nchar(grp)-2, nchar(grp))),
                           dam.var = obs.dpf$random[3,2],
                           dam.p = obs.dpf$random[3,7],
                           dam.perc = obs.dpf$random[3,3],
                           sire.var = obs.dpf$random[2,2],
                           sire.p = obs.dpf$random[2,7],
                           sire.perc = obs.dpf$random[2,3],
                           dam.sire.var = obs.dpf$random[1,2],
                           dam.sire.p = obs.dpf$random[1,7],
                           dam.sire.perc = obs.dpf$random[1,3],
                           residual.var = obs.dpf$other[1,2],
                           residual.perc = obs.dpf$other[1,3]) %>% 
    mutate_if(is.numeric, round, 4)
  })) %>% mutate(trait = "dpf")

## ADD
phenoVar.ADD.obs <- do.call(rbind, lapply(as.character(unique(hatch.ADD$group)), function(grp) {
  ## Filter to only a single group
  data.group <- hatch.ADD %>% filter(group == grp) %>% 
      select(family, dam, sire, block, ADD)
  
  obs.ADD <- observLmer(observ = data.group, dam = "dam", sire = "sire", response = "ADD")
  
  obs.ADD.df <- data.frame(group = substr(grp, 1, nchar(grp)-4),
                           temperature = as.numeric(substr(grp, nchar(grp)-2, nchar(grp))),
                           dam.var = obs.ADD$random[3,2],
                           dam.p = obs.ADD$random[3,7],
                           dam.perc = obs.ADD$random[3,3],
                           sire.var = obs.ADD$random[2,2],
                           sire.p = obs.ADD$random[2,7],
                           sire.perc = obs.ADD$random[2,3],
                           dam.sire.var = obs.ADD$random[1,2],
                           dam.sire.p = obs.ADD$random[1,7],
                           dam.sire.perc = obs.ADD$random[1,3],
                           residual.var = obs.ADD$other[1,2],
                           residual.perc = obs.ADD$other[1,3]) %>% 
    mutate_if(is.numeric, round, 4)  
  })) %>% mutate(trait = "ADD")


#### CREATE TEMPERATURE TREATMENT DATAFRAME ------------------------------------------------------

temp <- data.frame(group = c("LK-Whitefish", "LK-Whitefish", "LK-Whitefish", "LK-Whitefish",
                             "LK-Vendace", "LK-Vendace", "LK-Vendace", "LK-Vendace",
                             "LS-Cisco", "LS-Cisco", "LS-Cisco", "LS-Cisco",
                             "LO-Cisco", "LO-Cisco", "LO-Cisco", "LO-Cisco"),
                   temperature = c(rep(c(2.2, 4.0, 6.9, 8.0),2), rep(c(2.0, 4.4, 6.9, 8.9),2)),
                   temp.treatment = factor(rep(c("Coldest", "Cold", "Warm", "Warmest"), 4), 
                                           ordered = TRUE, levels = c("Coldest", "Cold", "Warm", "Warmest")))


#### CALCUALTE MEAN VARIANCE ACROSS TEMPERATURES -------------------------------------------------

phenoVar.embryo.mean <- bind_rows(phenoVar.survival.obs, phenoVar.dpf.obs, phenoVar.ADD.obs) %>% 
  group_by(group, trait) %>% 
  summarize(dam.perc.mean = mean(dam.perc),
            sire.perc.mean = mean(sire.perc),
            dam.sire.perc.mean = mean(dam.sire.perc),
            residual.perc.mean = mean(residual.perc)) %>% 
  pivot_longer(3:6, names_to = "component", values_to = "variance") %>% 
  filter(group != "LK-Whitefish") %>% 
  mutate(component = factor(component, ordered = TRUE,
                            levels = c("dam.perc.mean", "sire.perc.mean", "dam.sire.perc.mean", "residual.perc.mean"),
                            labels = c("Dam", "Sire", "Dam.Sire", "Error")),
         group = factor(group, ordered = TRUE, 
                        levels = c("LK-Vendace", "LS-Cisco", "LO-Cisco")))


#### CALCUALTE ERROR ACROSS TEMPERATURES ---------------------------------------------------------

phenoVar.embryo.error <- bind_rows(phenoVar.survival.obs, phenoVar.dpf.obs, phenoVar.ADD.obs) %>% 
  group_by(group, trait) %>% 
  summarize(dam.perc.se = sd(dam.perc)/sqrt(n()),
            sire.perc.se = sd(sire.perc)/sqrt(n()),
            dam.sire.perc.se = sd(dam.sire.perc)/sqrt(n()),
            residual.perc.se = sd(residual.perc)/sqrt(n())) %>% 
  pivot_longer(3:6, names_to = "component", values_to = "error") %>% 
  filter(group != "LK-Whitefish") %>% 
  mutate(component = factor(component, ordered = TRUE,
                            levels = c("dam.perc.se", "sire.perc.se", "dam.sire.perc.se", "residual.perc.se"),
                            labels = c("Dam", "Sire", "Dam.Sire", "Error")),
         group = factor(group, ordered = TRUE, 
                        levels = c("LK-Vendace", "LS-Cisco", "LO-Cisco")))


#### JOIN MEAN AND ERROR -------------------------------------------------------------------------

phenoVar.embryo.all <- left_join(phenoVar.embryo.mean, phenoVar.embryo.error)


#### COMBINE ALL TRAITS --------------------------------------------------------------------------

#phenoVar.embryo.all <- bind_rows(phenoVar.survival.obs, phenoVar.dpf.obs, phenoVar.ADD.obs) %>% 
#  left_join(temp) %>% 
#  filter(group != "LK-Whitefish") %>% 
#  select(group, temp.treatment, trait, dam.perc, sire.perc, dam.sire.perc, residual.perc) %>% 
#  pivot_longer(4:7, names_to = "component", values_to = "variance") %>% 
#  mutate(component = factor(component, ordered = TRUE,
#                            levels = c("dam.perc", "sire.perc", "dam.sire.perc", "residual.perc"),
#                            labels = c("Dam", "Sire", "Dam.Sire", "Error")),
#         component.trt = factor(interaction(component, temp.treatment), ordered = TRUE,
#                                levels = c("Dam.Coldest", "Dam.Cold", "Dam.Warm", "Dam.Warmest",
#                                           "Sire.Coldest", "Sire.Cold", "Sire.Warm", "Sire.Warmest",
#                                           "Dam.Sire.Coldest", "Dam.Sire.Cold", "Dam.Sire.Warm", "Dam.Sire.Warmest",
#                                           "Error.Coldest", "Error.Cold", "Error.Warm", "Error.Warmest")),
#         group = factor(group, ordered = TRUE, 
#                        levels = c("LK-Vendace", "LS-Cisco", "LO-Cisco"))) %>% 
#  left_join(phenoVar.embryo.boot)


#### CALCULATE CORRELATIONS ----------------------------------------------------------------------

## Embryo Survival
phenoVar.survival.cor <- phenoVar.survival.obs %>% 
  filter(group != "LK-Whitefish") %>% 
  group_by(group) %>% 
  summarize(dam.cor = cor(dam.perc, temperature),
            sire.cor = cor(sire.perc, temperature),
            dam.sire.cor = cor(dam.sire.perc, temperature),
            error.cor = cor(residual.perc, temperature)) %>% 
  mutate_if(is.numeric, round, 2) %>% 
  mutate(dam.cor.2 = ifelse(dam.cor >= 0.7, "POSITIVE", ifelse(dam.cor <= -0.7, "NEGATIVE", "NC")),
         sire.cor.2 = ifelse(sire.cor >= 0.7, "POSITIVE", ifelse(sire.cor <= -0.7, "NEGATIVE", "NC")),
         dam.sire.cor.2 = ifelse(dam.sire.cor >= 0.7, "POSITIVE", ifelse(dam.sire.cor <= -0.7, "NEGATIVE", "NC")),
         error.cor.2 = ifelse(error.cor >= 0.7, "POSITIVE", ifelse(error.cor <= -0.7, "NEGATIVE", "NC")))

## DPF
phenoVar.dpf.cor <- phenoVar.dpf.obs %>% 
  filter(group != "LK-Whitefish") %>% 
  group_by(group) %>% 
  summarize(dam.cor = cor(dam.perc, temperature),
            sire.cor = cor(sire.perc, temperature),
            dam.sire.cor = cor(dam.sire.perc, temperature),
            error.cor = cor(residual.perc, temperature)) %>% 
  mutate_if(is.numeric, round, 2) %>% 
  mutate(dam.cor.2 = ifelse(dam.cor >= 0.7, "POSITIVE", ifelse(dam.cor <= -0.7, "NEGATIVE", "NC")),
         sire.cor.2 = ifelse(sire.cor >= 0.7, "POSITIVE", ifelse(sire.cor <= -0.7, "NEGATIVE", "NC")),
         dam.sire.cor.2 = ifelse(dam.sire.cor >= 0.7, "POSITIVE", ifelse(dam.sire.cor <= -0.7, "NEGATIVE", "NC")),
         error.cor.2 = ifelse(error.cor >= 0.7, "POSITIVE", ifelse(error.cor <= -0.7, "NEGATIVE", "NC")))

## ADD
phenoVar.ADD.cor <- phenoVar.ADD.obs %>% 
  filter(group != "LK-Whitefish") %>% 
  group_by(group) %>% 
  summarize(dam.cor = cor(dam.perc, temperature),
            sire.cor = cor(sire.perc, temperature),
            dam.sire.cor = cor(dam.sire.perc, temperature),
            error.cor = cor(residual.perc, temperature)) %>% 
  mutate_if(is.numeric, round, 2) %>% 
  mutate(dam.cor.2 = ifelse(dam.cor >= 0.7, "POSITIVE", ifelse(dam.cor <= -0.7, "NEGATIVE", "NC")),
         sire.cor.2 = ifelse(sire.cor >= 0.7, "POSITIVE", ifelse(sire.cor <= -0.7, "NEGATIVE", "NC")),
         dam.sire.cor.2 = ifelse(dam.sire.cor >= 0.7, "POSITIVE", ifelse(dam.sire.cor <= -0.7, "NEGATIVE", "NC")),
         error.cor.2 = ifelse(error.cor >= 0.7, "POSITIVE", ifelse(error.cor <= -0.7, "NEGATIVE", "NC")))

