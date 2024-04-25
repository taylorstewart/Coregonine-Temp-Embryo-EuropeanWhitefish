#### CLEAR THE ENVIRONMENT FIRST -----------------------------------------------------------------

rm(list = ls(all.names = TRUE))


#### LOAD PACKAGES  ------------------------------------------------------------------------------

library(tidyverse)
library(readxl)
library(data.table)
library(fullfact)
library(parallel)


#### LOAD HATCHING DATA --------------------------------------------------------------------------

hatch.FR <- read_excel("data/Coregonine-Temperature-Experiment-EuropeFrance-Hatch.xlsx", sheet = "hatching") %>% 
  select(population, species_form, family, male, female, female_tl_mm, female_fm_g, male_tl_mm, male_fm_g, block, no, temperature, plate, eye, hatch, dpf, ADD, include.incubation) %>% 
  mutate(female = factor(female),
         male = factor(male),
         family = factor(family),
         block = factor(block),
         # Create a variable with population and species combined
         population = factor(interaction(population, species_form)))

hatch.FI <- read_excel("data/Coregonine-Temperature-Experiment-EuropeFinland-Hatch.xlsx", sheet = "hatching") %>% 
  select(population, species_form = species, family, male, female, female_tl_mm, female_fm_g, male_tl_mm, male_fm_g, block, no, temperature, plate, eye, hatch, dpf, ADD, include.incubation) %>% 
  filter(temperature %in% c(6.9, 8), species_form == "lavaretus") %>% 
  mutate(female = factor(female),
         male = factor(male),
         family = factor(family),
         block = factor(block),
         # Create a variable with population and species combined
         population = factor("konnevesi.littoral"))

hatch <- bind_rows(hatch.FR, hatch.FI) %>% 
  mutate(temperature = factor(temperature, ordered = TRUE, levels = c(6.9, 7, 8, 9), labels = c(7, 7, 9, 9)),
         population = factor(population, ordered = TRUE, levels = c("konnevesi.littoral", "constance.littoral", "constance.pelagic", "leman.littoral", "bourget.littoral")),
         group = interaction(population, temperature)) %>% 
  filter(include.incubation == "y") %>% 
  rename(sire = male, dam = female)

## Clean up environment
rm(hatch.FR, hatch.FI)


#### FILTER TO EACH TRAITS' DATASET --------------------------------------------------------------

## filter to only eyed embryos
hatch.survival <- hatch %>% filter(eye != 0)

## filter to only hatched embryos
hatch.survived <- hatch %>% filter(!is.na(dpf), hatch == 1)


#### STATISTICAL ANALYSIS - GENERATE OBSERVED VARIANCES ------------------------------------------

## Embryo Survival
phenoVar.survival.obs <- do.call(rbind, lapply(as.character(unique(hatch.survival$group)), function(grp) {
  ## Filter to only a single group
  data.group <- hatch.survival %>% filter(group == grp) %>% 
      select(family, dam, sire, block, hatch)
  
  obs.survival <- observGlmer(observ = data.group, dam = "dam", sire = "sire", response = "hatch",
                              fam_link = binomial(link = "logit"))
  
  obs.survival.df <- data.frame(group = substr(grp, 1, nchar(grp)-2),
                                temperature = as.numeric(substr(grp, nchar(grp), nchar(grp))),
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
phenoVar.dpf.obs <- do.call(rbind, lapply(as.character(unique(hatch.survived$group)), function(grp) {
  ## Filter to only a single group
  data.group <- hatch.survived %>% filter(group == grp) %>% 
      select(family, dam, sire, block, dpf)
    
  obs.dpf <- observLmer(observ = data.group, dam = "dam", sire = "sire", response = "dpf")
  
  obs.dpf.df <- data.frame(group = substr(grp, 1, nchar(grp)-2),
                           temperature = as.numeric(substr(grp, nchar(grp), nchar(grp))),
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
phenoVar.ADD.obs <- do.call(rbind, lapply(as.character(unique(hatch.survived$group)), function(grp) {
  ## Filter to only a single group
  data.group <- hatch.survived %>% filter(group == grp) %>% 
      select(family, dam, sire, block, ADD)
  
  obs.ADD <- observLmer(observ = data.group, dam = "dam", sire = "sire", response = "ADD")
  
  obs.ADD.df <- data.frame(group = substr(grp, 1, nchar(grp)-2),
                           temperature = as.numeric(substr(grp, nchar(grp), nchar(grp))),
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

temp <- data.frame(group = c("konnevesi.littoral", "konnevesi.littoral",
                             "constance.pelagic", "constance.pelagic",
                             "constance.littoral", "constance.littoral",
                             "leman.littoral", "leman.littoral",
                             "bourget.littoral", "bourget.littoral"),
                   temperature = factor(c(7, 8, rep(c(7, 9), 4)), ordered = TRUE, levels = c(7, 8, 9)))


#### CALCUALTE MEAN VARIANCE ACROSS TEMPERATURES -------------------------------------------------

phenoVar.embryo.mean <- bind_rows(phenoVar.survival.obs, phenoVar.dpf.obs, phenoVar.ADD.obs) %>% 
  group_by(group, trait) %>% 
  summarize(dam.perc.mean = mean(dam.perc),
            sire.perc.mean = mean(sire.perc),
            dam.sire.perc.mean = mean(dam.sire.perc),
            residual.perc.mean = mean(residual.perc)) %>% 
  pivot_longer(3:6, names_to = "component", values_to = "variance") %>% 
  mutate(component = factor(component, ordered = TRUE,
                            levels = c("dam.perc.mean", "sire.perc.mean", "dam.sire.perc.mean", "residual.perc.mean"),
                            labels = c("Dam", "Sire", "Dam.Sire", "Error")),
         group = factor(group, ordered = TRUE, levels = c("konnevesi.littoral", "constance.pelagic", "constance.littoral", "leman.littoral", "bourget.littoral"),
                        labels = c("S. Konnevesi", "Constance\n(Pelagic)", "Constance\n(Littoral)", "Geneva", "Bourget")))


#### CALCUALTE ERROR ACROSS TEMPERATURES ---------------------------------------------------------

phenoVar.embryo.error <- bind_rows(phenoVar.survival.obs, phenoVar.dpf.obs, phenoVar.ADD.obs) %>% 
  group_by(group, trait) %>% 
  summarize(dam.perc.se = sd(dam.perc)/sqrt(n()),
            sire.perc.se = sd(sire.perc)/sqrt(n()),
            dam.sire.perc.se = sd(dam.sire.perc)/sqrt(n()),
            residual.perc.se = sd(residual.perc)/sqrt(n())) %>% 
  pivot_longer(3:6, names_to = "component", values_to = "error") %>% 
  mutate(component = factor(component, ordered = TRUE,
                            levels = c("dam.perc.se", "sire.perc.se", "dam.sire.perc.se", "residual.perc.se"),
                            labels = c("Dam", "Sire", "Dam.Sire", "Error")),
         group = factor(group, ordered = TRUE, levels = c("konnevesi.littoral", "constance.pelagic", "constance.littoral", "leman.littoral", "bourget.littoral"),
                        labels = c("S. Konnevesi", "Constance\n(Pelagic)", "Constance\n(Littoral)", "Geneva", "Bourget")))


#### JOIN MEAN AND ERROR -------------------------------------------------------------------------

phenoVar.embryo.all <- left_join(phenoVar.embryo.mean, phenoVar.embryo.error) %>% 
  mutate(trait = factor(trait, ordered = TRUE, levels = c("survival", "dpf", "ADD"),
                        labels = c("Embryo Survival", "Days Post-fertilization", "Accumulated Degrees-days")))


#### VISUALIZATION -------------------------------------------------------------------------------

## Create base plot
ggplot(phenoVar.embryo.all, aes(x = group, y = variance, group = component, fill = component)) + 
  geom_bar(stat = "identity", size = 0.2, position = position_dodge(0.9), color = "black") +
  geom_errorbar(aes(ymin = ifelse(variance - error < 0, 0, variance - error), 
                    ymax = ifelse(variance + error > 100, 99.5, variance + error)), 
                position = position_dodge(0.9), size = 0.5, width = 0.4, color = "gray15", show.legend = FALSE) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0.5)) +
  scale_fill_manual(values = c("#7bccc4", "#f0f9e8", "#bae4bc", "#2b8cbe"),
                    labels = c("Female  ", "Male  ", "Female x Male  ", "Residual Error")) +
  labs(y = "% of Total Phenotypic Variation", x = "Population") +
  theme_bw() +
  theme(axis.title.x = element_text(color = "Black", size = 20, margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 20, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 17),
        axis.text.y = element_text(size = 17),
        axis.ticks.length = unit(2, 'mm'),
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1, 'cm'),
        legend.position = "top",
        strip.text = element_text(size = 18),
        strip.background = element_rect(color = "transparent", fill = "white"),
        panel.spacing = unit(1.5, "lines"),
        plot.margin = unit(c(5, 5, 5, 5), 'mm')) +
  facet_wrap(~trait)

## Save figure
ggsave("figures/Fig2-PhenoVar.png", width = 14, height = 7, dpi = 600)
