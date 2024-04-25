#### CLEAR THE ENVIRONMENT FIRST ---------------------------------------------

rm(list = ls(all.names = TRUE))


#### LOAD PACKAGES -----------------------------------------------------------

library(tidyverse)
library(readxl)
library(lme4)
library(lmerTest)
library(afex)
library(buildmer)
library(ggplot2)
library(gridExtra)
library(grid)
library(cowplot)
library(emmeans)


#### LOAD HATCHING DATA ------------------------------------------------------

hatch.FR <- read_excel("data/Coregonine-Temperature-Experiment-EuropeFrance-Hatch.xlsx", sheet = "hatching") %>% 
  select(population, species_form, family, male, female, female_tl_mm, female_fm_g, male_tl_mm, male_fm_g, block, no, temperature, plate, eye, hatch, dpf, ADD, include.incubation) %>% 
  mutate(female = factor(female),
         male = factor(male),
         family = factor(family),
         block = factor(block),
         # Create a variable with population and species combined
         group = factor(interaction(population, species_form)))

hatch.FI <- read_excel("data/Coregonine-Temperature-Experiment-EuropeFinland-Hatch.xlsx", sheet = "hatching") %>% 
  select(population, species_form = species, family, male, female, female_tl_mm, female_fm_g, male_tl_mm, male_fm_g, block, no, temperature, plate, eye, hatch, dpf, ADD, include.incubation) %>% 
  filter(temperature %in% c(6.9, 8), species_form == "lavaretus") %>% 
  mutate(female = factor(female),
         male = factor(male),
         family = factor(family),
         block = factor(block),
         # Create a variable with population and species combined
         group = factor("konnevesi.littoral"))

hatch <- bind_rows(hatch.FR, hatch.FI) %>% 
  mutate(temperature = factor(temperature, ordered = TRUE, levels = c(6.9, 7, 8, 9), labels = c(7, 7, 8, 9)),
         group = factor(group, ordered = TRUE, levels = c("konnevesi.littoral", "constance.littoral", "constance.pelagic", "leman.littoral", "bourget.littoral")),
         trans.dpf = dpf^(1/4),
         trans.ADD = ADD^(1/4)) %>% 
  filter(include.incubation == "y")



str(hatch)

prop_eye <- hatch %>% group_by(group, temperature) %>% 
  summarize(n = n(),
            eye = sum(eye)) %>% ungroup() %>% 
  mutate(prop = eye/n)


